"""Unsupervised representation learning on the k-space, then a label-free check.

Every earlier model was supervised. This asks a different question: does the
k-space contain gaze structure that a classifier failed to *extract*, or does it
contain none at all? Labels are used only to EVALUATE what was learned, never to
fit it.

Two methods, both on the per-direction-cluster detrended residual. That matters:
the dominant variance in raw readouts is the spoke direction (a supervised probe
predicts it at 0.9865), so anything unsupervised run on raw data rediscovers the
trajectory and nothing else.

  kmeans      cluster window embeddings into k=4, match clusters to gaze bins by
              optimal assignment, report balanced accuracy.
  autoencoder FiLM/TCN encoder -> bottleneck -> decoder, trained only to
              reconstruct its input. Freeze, then fit a linear probe on the
              frozen embeddings (the standard linear-evaluation protocol).

Both are scored against a circular-shift null: the clustering / representation
is unchanged, only the evaluation labels move. That isolates real
label-structure agreement from chance agreement between two blocky sequences.
"""
import argparse
import glob
import time

import numpy as np
import scipy.io as sio
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NR, NB, NV = 79894, 160, 8
NSEG, NSHOT, NOFF = 44, 1872, 14
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
SHIFTS = (11000, 23000, 37000, 51000)
DIRS = np.load(OUT / 'dirs.npy')


def gaze_labels(sub, kind='filtered'):
    lab = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/{kind}/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        lab[m[NOFF:, 1:].ravel().astype(bool)] = r
    return lab


def detrended(sub):
    p = OUT / 'detrended' / f'sub-{sub:03d}.npy'
    if not p.exists():
        raise FileNotFoundError(f'{p} - build it with train_gaze.py first')
    return np.load(p, mmap_mode='r')


def windows(lab, W, stride, purity=0.7):
    """Window start indices and their majority label, or -1 if impure."""
    oh = np.stack([(lab == k) for k in range(4)], 1).astype(np.int32)
    cl = np.concatenate([np.zeros((1, 4), np.int64), np.cumsum(oh, 0)])
    starts, ys = [], []
    for a in range(0, NR - W + 1, stride):
        c = cl[a + W] - cl[a]
        n = c.sum()
        starts.append(a)
        ys.append(int(np.argmax(c)) if n >= 0.3 * W and c.max() >= purity * n else -1)
    return np.array(starts), np.array(ys)


def matched_bacc(pred, true, k=4):
    """Balanced accuracy after the best cluster->class assignment."""
    m = true >= 0
    pred, true = pred[m], true[m]
    C = np.zeros((k, 4))
    for p, t in zip(pred, true):
        C[p, t] += 1
    r, c = linear_sum_assignment(-C)
    mapping = dict(zip(r, c))
    yh = np.array([mapping.get(p, -1) for p in pred])
    hit = np.zeros(4); tot = np.zeros(4)
    for j in range(4):
        s = true == j
        tot[j] = s.sum(); hit[j] = (yh[s] == j).sum()
    ok = tot > 0
    return float(np.mean(hit[ok] / tot[ok])) if ok.any() else np.nan


# ------------------------------------------------------------------ autoencoder

class AE(nn.Module):
    """Per-readout conv autoencoder, FiLM-conditioned on the spoke direction, with
    a temporal trunk so the code can carry structure across readouts."""

    def __init__(self, cin=2 * NV, d=32):
        super().__init__()
        self.dirmlp = nn.Sequential(nn.Linear(3, 32), nn.GELU(), nn.Linear(32, 32), nn.GELU())
        ch = [cin, 16, 32, 64]
        self.convs = nn.ModuleList([nn.Conv1d(ch[i], ch[i + 1], 7, 2, 3) for i in range(3)])
        self.norms = nn.ModuleList([nn.GroupNorm(4, ch[i + 1]) for i in range(3)])
        self.films = nn.ModuleList([nn.Linear(32, 2 * ch[i + 1]) for i in range(3)])
        self.to_code = nn.Linear(64 * (NB // 8), d)
        self.tcn = nn.Sequential(
            nn.Conv1d(d, d, 3, padding=2, dilation=2), nn.GroupNorm(4, d), nn.GELU(),
            nn.Conv1d(d, d, 3, padding=8, dilation=8), nn.GroupNorm(4, d), nn.GELU())
        self.from_code = nn.Linear(d, 64 * (NB // 8))
        self.deconv = nn.Sequential(
            nn.ConvTranspose1d(64, 32, 4, 2, 1), nn.GELU(),
            nn.ConvTranspose1d(32, 16, 4, 2, 1), nn.GELU(),
            nn.ConvTranspose1d(16, cin, 4, 2, 1))

    def encode(self, x, dirs):
        B, W, C, N = x.shape
        e = self.dirmlp(dirs.reshape(B * W, 3))
        h = x.reshape(B * W, C, N)
        for conv, norm, film in zip(self.convs, self.norms, self.films):
            h = norm(conv(h))
            g, b = film(e).chunk(2, -1)
            h = F.gelu(h * (1 + g[..., None]) + b[..., None])
        z = self.to_code(h.flatten(1)).reshape(B, W, -1).permute(0, 2, 1)
        return self.tcn(z) + z                       # [B, d, W]

    def forward(self, x, dirs):
        B, W = x.shape[0], x.shape[1]
        z = self.encode(x, dirs)
        h = self.from_code(z.permute(0, 2, 1).reshape(B * W, -1))
        return self.deconv(h.reshape(B * W, 64, NB // 8)).reshape(B, W, -1, NB)


class WinData(torch.utils.data.Dataset):
    def __init__(self, subs, W, stride):
        self.items, self.X = [], {}
        for s in subs:
            self.X[s] = detrended(s)
            for a in range(0, NR - W + 1, stride):
                self.items.append((s, a))
        self.W = W

    def __len__(self):
        return len(self.items)

    def __getitem__(self, i):
        s, a = self.items[i]
        x = np.asarray(self.X[s][a:a + self.W])
        x = np.concatenate([x.real, x.imag], 2)
        return (torch.from_numpy(x.astype(np.float32)).permute(0, 2, 1),
                torch.from_numpy(DIRS[a:a + self.W].astype(np.float32)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--subs', type=int, nargs='*', default=list(KEEP))
    ap.add_argument('--window', type=int, default=256)
    ap.add_argument('--epochs', type=int, default=25)
    ap.add_argument('--code', type=int, default=32)
    a = ap.parse_args()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    W, stride = a.window, a.window

    # ---------------- method 1: k-means on PCA window embeddings --------------
    print('=== k-means on detrended window embeddings (k=4) ===')
    print(f'  {"sub":8s} {"windows":>8s} {"true":>8s} {"null":>8s} {"z":>7s}')
    kz = []
    for s in a.subs:
        X = np.asarray(detrended(s)).reshape(NR, -1)
        Xr = np.column_stack([X.real, X.imag]).astype(np.float32)
        st, y = windows(gaze_labels(s), W, stride)
        emb = np.stack([Xr[i:i + W].mean(0) for i in st])
        emb = StandardScaler().fit_transform(emb)
        u, sv, vt = np.linalg.svd(emb - emb.mean(0), full_matrices=False)
        emb = (emb @ vt[:16].T)
        km = KMeans(4, n_init=10, random_state=0).fit_predict(emb)
        true = matched_bacc(km, y)
        nul = []
        for sh in SHIFTS:
            _, ys = windows(np.roll(gaze_labels(s), sh), W, stride)
            nul.append(matched_bacc(km, ys))
        nul = np.array(nul, float)
        z = (true - nul.mean()) / (nul.std(ddof=1) + 1e-9)
        kz.append(z)
        print(f'  sub-{s:03d} {len(st):8d} {true:8.4f} {nul.mean():8.4f} {z:+7.2f}', flush=True)
    print(f'  median z {np.median(kz):+.2f}   (chance balanced accuracy 0.25)\n')

    # ---------------- method 2: autoencoder + frozen linear probe -------------
    print(f'=== autoencoder ({a.code}-d code), then linear probe on frozen features ===')
    ds = WinData(a.subs, W, stride)
    dl = torch.utils.data.DataLoader(ds, batch_size=8, shuffle=True, num_workers=4)
    net = AE(d=a.code).to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=3e-4, weight_decay=1e-4)
    for ep in range(a.epochs):
        net.train(); tot = n = 0; t0 = time.time()
        for x, d in dl:
            x, d = x.to(dev), d.to(dev)
            loss = F.mse_loss(net(x, d), x)
            opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.item(); n += 1
        print(f'  ep {ep:3d} recon MSE {tot/n:.5f}  [{time.time()-t0:.0f}s]', flush=True)

    net.eval()
    print(f'\n  {"sub":8s} {"probe true":>11s} {"probe null":>11s} {"z":>7s}')
    pz = []
    for s in a.subs:
        X = detrended(s)
        st, y = windows(gaze_labels(s), W, stride)
        Z = []
        with torch.no_grad():
            for i in range(0, len(st), 16):
                xs, ds = [], []
                for aa in st[i:i + 16]:
                    x = np.asarray(X[aa:aa + W])
                    xs.append(np.concatenate([x.real, x.imag], 2).astype(np.float32))
                    ds.append(DIRS[aa:aa + W].astype(np.float32))
                xb = torch.from_numpy(np.stack(xs)).permute(0, 1, 3, 2).to(dev)
                db = torch.from_numpy(np.stack(ds)).to(dev)
                Z.append(net.encode(xb, db).mean(-1).cpu().numpy())
        Z = StandardScaler().fit_transform(np.concatenate(Z))
        def probe(lbl):
            m = lbl >= 0
            if m.sum() < 40 or len(np.unique(lbl[m])) < 3:
                return np.nan
            cut = int(0.7 * len(lbl))
            tr = m & (np.arange(len(lbl)) < cut)
            te = m & (np.arange(len(lbl)) >= cut)
            if tr.sum() < 20 or te.sum() < 20:
                return np.nan
            clf = LogisticRegression(max_iter=2000, C=0.1).fit(Z[tr], lbl[tr])
            p = clf.predict(Z[te])
            hit = np.zeros(4); tot = np.zeros(4)
            for j in range(4):
                sel = lbl[te] == j
                tot[j] = sel.sum(); hit[j] = (p[sel] == j).sum()
            ok = tot > 0
            return float(np.mean(hit[ok] / tot[ok]))
        true = probe(y)
        nul = np.array([probe(windows(np.roll(gaze_labels(s), sh), W, stride)[1])
                        for sh in SHIFTS], float)
        z = (true - np.nanmean(nul)) / (np.nanstd(nul, ddof=1) + 1e-9)
        pz.append(z)
        print(f'  sub-{s:03d} {true:11.4f} {np.nanmean(nul):11.4f} {z:+7.2f}', flush=True)
    print(f'  median z {np.nanmedian(pz):+.2f}   (chance 0.25)')


if __name__ == '__main__':
    main()
