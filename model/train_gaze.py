"""Gaze decoding from ROI-PCA-8 k-space: sliding windows -> per-readout classes.

Input  : W consecutive readouts, each 160 k-space samples x 8 virtual channels
         (complex -> 16 real planes), plus the spoke direction of that readout.
Output : per-readout class probabilities over the 4 gaze directions, so a whole
         scan comes back as an array of readouts labelled up/down/left/right.
         Extra classes (blink, saccade) drop in by widening the head.

Design choices are forced by measurements in DESIGN.md, not taste:

  direction conditioning  a readout's profile depends mostly on its spoke
                          direction, so the encoder is FiLM-modulated by it;
                          without this the dominant input variance is the
                          trajectory rotating, not the eye moving.
  per-cluster detrend     slow drift let earlier probes decode *time*; the same
                          detrend that fixed them is applied here, cached.
  within-subject splits   the ROI-PCA transform is per subject, so virtual
                          channel 3 of sub-001 is unrelated to sub-004's. A
                          shared trunk gets a small per-subject input adapter.
  balanced + shift null   raw accuracy is unreadable here (chance is the
                          window-majority rate). Always report balanced
                          accuracy against a circular-shift null.
"""
import argparse
import glob
import time
from pathlib import Path

import numpy as np
import scipy.io as sio
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
CACHE = OUT / 'detrended'
NR, NB, NV = 79894, 160, 8
NSEG, NSHOT, NOFF = 44, 1872, 14
NC, MAVG = 600, 11
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)

DIRS = np.load(OUT / 'dirs.npy')
CL = cKDTree(DIRS[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(DIRS)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]


# --------------------------------------------------------------------------- data

def detrended(sub):
    """Per-cluster local detrend, cached. Removes drift slower than ~55 s."""
    CACHE.mkdir(exist_ok=True)
    p = CACHE / f'sub-{sub:03d}.npy'
    if p.exists():
        return np.load(p, mmap_mode='r')
    X = np.load(OUT / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy').reshape(NR, -1)
    X = (X / (np.abs(X).mean() + 1e-12)).astype(np.complex64)
    R = np.empty_like(X)
    ker = np.ones(MAVG) / MAVG
    for idx in CIDX:
        blk = X[idx]
        if len(idx) >= MAVG:
            pad = np.concatenate([blk[:MAVG // 2][::-1], blk, blk[-(MAVG // 2):][::-1]])
            loc = np.apply_along_axis(lambda v: np.convolve(v, ker, 'valid'), 0, pad)
            R[idx] = blk - loc
        else:
            R[idx] = blk - blk.mean(0)
    np.save(p, R.reshape(NR, NB, NV))
    return np.load(p, mmap_mode='r')


def gaze_labels(sub, kind='filtered'):
    lab = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/{kind}/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        lab[m[NOFF:, 1:].ravel().astype(bool)] = r
    return lab


class Windows(torch.utils.data.Dataset):
    """Sliding windows of readouts from one subject, restricted to a time range."""

    def __init__(self, sub, lo, hi, W=384, stride=64, kind='filtered', sidx=0):
        self.X = detrended(sub)
        self.y = gaze_labels(sub, kind)
        self.W, self.sub, self.sidx = W, sub, sidx
        self.starts = np.arange(lo, max(lo, hi - W + 1), stride)

    def __len__(self):
        return len(self.starts)

    def __getitem__(self, i):
        a = int(self.starts[i]); b = a + self.W
        x = np.asarray(self.X[a:b])                      # [W, NB, NV] complex64
        x = np.concatenate([x.real, x.imag], axis=2)     # [W, NB, 2*NV]
        return (torch.from_numpy(x.astype(np.float32)).permute(0, 2, 1),   # [W, C, NB]
                torch.from_numpy(DIRS[a:b].astype(np.float32)),            # [W, 3]
                torch.from_numpy(self.y[a:b].astype(np.int64)),
                self.sidx)


# -------------------------------------------------------------------------- model

class ReadoutEncoder(nn.Module):
    """1D conv over the k-space sample axis, FiLM-modulated by the spoke direction."""

    def __init__(self, cin=2 * NV, d=128, ndir=32):
        super().__init__()
        self.dir = nn.Sequential(nn.Linear(3, 64), nn.GELU(), nn.Linear(64, ndir), nn.GELU())
        ch = [cin, 32, 64, 128]
        self.convs = nn.ModuleList(
            [nn.Conv1d(ch[i], ch[i + 1], 7, stride=2, padding=3) for i in range(3)])
        self.norms = nn.ModuleList([nn.GroupNorm(8, ch[i + 1]) for i in range(3)])
        self.films = nn.ModuleList([nn.Linear(ndir, 2 * ch[i + 1]) for i in range(3)])
        self.head = nn.Linear(128, d)

    def forward(self, x, dirs):
        # x [N, C, NB], dirs [N, 3]
        e = self.dir(dirs)
        for conv, norm, film in zip(self.convs, self.norms, self.films):
            x = norm(conv(x))
            g, b = film(e).chunk(2, -1)
            x = F.gelu(x * (1 + g[..., None]) + b[..., None])
        return self.head(x.mean(-1))


class TemporalTrunk(nn.Module):
    """Dilated causal-free TCN over the readout axis."""

    def __init__(self, d=128, nlayer=5):
        super().__init__()
        self.blocks = nn.ModuleList()
        for i in range(nlayer):
            dil = 2 ** i
            self.blocks.append(nn.Sequential(
                nn.Conv1d(d, d, 3, padding=dil, dilation=dil),
                nn.GroupNorm(8, d), nn.GELU(),
                nn.Conv1d(d, d, 1), nn.GroupNorm(8, d), nn.GELU()))

    def forward(self, h):                      # [B, d, W]
        for blk in self.blocks:
            h = h + blk(h)
        return h


class GazeNet(nn.Module):
    def __init__(self, nsub, d=128, ncls=4):
        super().__init__()
        # per-subject channel adapter: the ROI-PCA basis is subject-specific
        self.adapt = nn.Parameter(torch.stack([torch.eye(2 * NV) for _ in range(nsub)]))
        self.enc = ReadoutEncoder(d=d)
        self.trunk = TemporalTrunk(d)
        self.head = nn.Conv1d(d, ncls, 1)

    def forward(self, x, dirs, sidx):
        B, W, C, NBs = x.shape
        A = self.adapt[sidx]                              # [B, C, C]
        x = torch.einsum('bwcn,bdc->bwdn', x, A)
        h = self.enc(x.reshape(B * W, C, NBs), dirs.reshape(B * W, 3)).reshape(B, W, -1)
        h = self.trunk(h.permute(0, 2, 1))
        return self.head(h).permute(0, 2, 1)              # [B, W, ncls]


# ----------------------------------------------------------------------- training

def balanced_acc(pred, true):
    hit = np.zeros(4); tot = np.zeros(4)
    for k in range(4):
        m = true == k
        tot[k] = m.sum(); hit[k] = (pred[m] == k).sum()
    ok = tot > 0
    return float(np.mean(hit[ok] / tot[ok])) if ok.any() else float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--subs', type=int, nargs='*', default=list(KEEP))
    ap.add_argument('--window', type=int, default=384)
    ap.add_argument('--epochs', type=int, default=20)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--lr', type=float, default=3e-4)
    ap.add_argument('--kind', default='filtered')
    ap.add_argument('--smoke', action='store_true')
    a = ap.parse_args()

    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    sidx = {s: i for i, s in enumerate(a.subs)}
    # within-subject temporal split: first 70% train, last 30% test
    tr = torch.utils.data.ConcatDataset(
        [Windows(s, 0, int(0.7 * NR), a.window, 64, a.kind, sidx[s]) for s in a.subs])
    te = torch.utils.data.ConcatDataset(
        [Windows(s, int(0.7 * NR), NR, a.window, a.window, a.kind, sidx[s]) for s in a.subs])
    print(f'train windows {len(tr)}, test windows {len(te)}, device {dev}')

    net = GazeNet(len(a.subs)).to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=a.lr, weight_decay=1e-4)
    nep = 1 if a.smoke else a.epochs

    def loader(ds, shuf, order):
        return torch.utils.data.DataLoader(ds, batch_size=a.batch, shuffle=shuf,
                                           num_workers=4, drop_last=False)

    for ep in range(nep):
        net.train(); t0 = time.time(); tot = n = 0
        for bi, (x, d, y, si) in enumerate(loader(tr, True, None)):
            x, d, y, si = x.to(dev), d.to(dev), y.to(dev), si.to(dev)
            logit = net(x, d, si)
            loss = F.cross_entropy(logit.reshape(-1, 4), y.reshape(-1), ignore_index=-1)
            opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.item(); n += 1
            if a.smoke and bi >= 3:
                break
        net.eval(); P = []; T = []
        with torch.no_grad():
            for bi, (x, d, y, si) in enumerate(loader(te, False, None)):
                p = net(x.to(dev), d.to(dev), si.to(dev)).argmax(-1).cpu().numpy().ravel()
                t = y.numpy().ravel()
                m = t >= 0
                P.append(p[m]); T.append(t[m])
                if a.smoke and bi >= 3:
                    break
        P = np.concatenate(P); T = np.concatenate(T)
        print(f'ep {ep:3d}  loss {tot/max(n,1):.4f}  test balanced acc '
              f'{balanced_acc(P, T):.4f}  (chance 0.25)  [{time.time()-t0:.0f}s]', flush=True)
        if a.smoke:
            break
    print('done')


if __name__ == '__main__':
    main()
