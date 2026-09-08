"""Shrunken gaze model: fixes the three faults of the first run.

  1. capacity   - encoder 16/32/64 (was 32/64/128), d=64 (was 128), dropout.
  2. receptive  - dilations are chosen so the TCN actually spans its window.
                  The first run used 1,2,4,8,16 = 63 readouts against a 384
                  readout window, so the trunk saw 16% of its own input.
  3. protocol   - real 3-way split (60/20/20 in time), model selected on val,
                  test touched once. The first run printed test every epoch.

Windows are strided by W//2 rather than 64, so they overlap 2x instead of 6x.
"""
import argparse, glob, time
import numpy as np
import scipy.io as sio
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
NR, NB = 79894, 160
NV = 8
NSEG, NSHOT, NOFF = 44, 1872, 14
NC, MAVG = 600, 11
NCLS = 4
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
DIRS = np.load(OUT / 'dirs.npy')
CL = cKDTree(DIRS[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(DIRS)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]


SRC = {'ROI-PCA_8': (OUT, 'sub-{s:03d}_band160-320_ROI-PCA_8.npy', 8),
       'SVD_8':  (Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/svd'),
                  'sub-{s:03d}_band160-320_SVD_8.npy', 8),
       'SVD_16': (Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/svd'),
                  'sub-{s:03d}_band160-320_SVD_16.npy', 16),
       'raw52': (Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband_raw52'),
                 'sub-{s:03d}_band160-320_raw52.npy', 52)}
VARIANT = 'ROI-PCA_8'


def detrended(sub):
    """Per-cluster local detrend, cached per variant."""
    root, pat, _ = SRC[VARIANT]
    cache = OUT / f'detrended_{VARIANT}'
    cache.mkdir(exist_ok=True)
    p = cache / f'sub-{sub:03d}.npy'
    if VARIANT == 'ROI-PCA_8' and (OUT / 'detrended' / f'sub-{sub:03d}.npy').exists():
        return np.load(OUT / 'detrended' / f'sub-{sub:03d}.npy', mmap_mode='r')
    if p.exists():
        return np.load(p, mmap_mode='r')
    X = np.load(root / pat.format(s=sub)).reshape(NR, -1)
    X = (X / (np.abs(X).mean() + 1e-12)).astype(np.complex64)
    R = np.empty_like(X)
    ker = np.ones(MAVG) / MAVG
    for idx in CIDX:
        blk = X[idx]
        if len(idx) >= MAVG:
            pad = np.concatenate([blk[:MAVG // 2][::-1], blk, blk[-(MAVG // 2):][::-1]])
            R[idx] = blk - np.apply_along_axis(
                lambda v: np.convolve(v, ker, 'valid'), 0, pad)
        else:
            R[idx] = blk - blk.mean(0)
    nv = SRC[VARIANT][2]
    np.save(p, R.reshape(NR, NB, nv))
    return np.load(p, mmap_mode='r')


def gaze_labels(sub, kind='filtered'):
    lab = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/{kind}/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        lab[m[NOFF:, 1:].ravel().astype(bool)] = r
    return lab


def _bin(sub, name):
    p = f'{STUDY}/sub-{sub:03d}/recon/bins/{name}/eMask_th0.75_winLen10.mat'
    m = np.asarray(sio.loadmat(p)['eMaskN']).squeeze().reshape(NSHOT, NSEG)
    return m[NOFF:, 1:].ravel().astype(bool)


def motion_labels(sub, kind=None):
    """0 = eye still, 1 = moving (saccade / blink / tracking loss).

    Unlike gaze position, motion has a direct mechanism: it corrupts the phase of
    the readout acquired during it, whatever direction the eye points. So this
    does not require resolving position -- which is what the k-space cannot do.
    Readouts in neither class are transitional and are left unlabelled.
    """
    lab = np.full(NR, -1, np.int8)
    mot = np.zeros(NR, bool)
    for n in ('saccade-only', 'blink-event', 'tracking-loss'):
        mot |= _bin(sub, n)
    lab[_bin(sub, 'fixation-ok') & ~mot] = 0
    lab[mot] = 1
    return lab


class Windows(torch.utils.data.Dataset):
    def __init__(self, sub, lo, hi, W, stride, kind, sidx, task='gaze'):
        lab = gaze_labels if task == 'gaze' else motion_labels
        self.X, self.y, self.W, self.sidx = detrended(sub), lab(sub, kind), W, sidx
        self.starts = np.arange(lo, max(lo, hi - W + 1), stride)

    def __len__(self):
        return len(self.starts)

    def __getitem__(self, i):
        a = int(self.starts[i]); b = a + self.W
        x = np.asarray(self.X[a:b])
        x = np.concatenate([x.real, x.imag], 2)
        return (torch.from_numpy(x.astype(np.float32)).permute(0, 2, 1),
                torch.from_numpy(DIRS[a:b].astype(np.float32)),
                torch.from_numpy(self.y[a:b].astype(np.int64)), self.sidx)


class Net(nn.Module):
    def __init__(self, nsub, W, d=64, drop=0.3):
        super().__init__()
        self.adapt = nn.Parameter(torch.stack([torch.eye(2 * NV) for _ in range(nsub)]))
        self.dirmlp = nn.Sequential(nn.Linear(3, 32), nn.GELU(), nn.Linear(32, 32), nn.GELU())
        ch = [2 * NV, 16, 32, 64]
        self.convs = nn.ModuleList([nn.Conv1d(ch[i], ch[i+1], 7, 2, 3) for i in range(3)])
        self.norms = nn.ModuleList([nn.GroupNorm(4, ch[i+1]) for i in range(3)])
        self.films = nn.ModuleList([nn.Linear(32, 2*ch[i+1]) for i in range(3)])
        self.proj = nn.Linear(64, d)
        # dilations chosen so 1 + 2*sum(dil) >= W
        dil, s = [], 0
        while 1 + 2 * s < W:
            dil.append(2 ** len(dil)); s += dil[-1]
        self.blocks = nn.ModuleList([nn.Sequential(
            nn.Conv1d(d, d, 3, padding=k, dilation=k), nn.GroupNorm(4, d), nn.GELU(),
            nn.Dropout(drop), nn.Conv1d(d, d, 1), nn.GroupNorm(4, d), nn.GELU())
            for k in dil])
        self.rf = 1 + 2 * s
        self.head = nn.Conv1d(d, NCLS, 1)

    def forward(self, x, dirs, si):
        B, W, C, N = x.shape
        x = torch.einsum('bwcn,bdc->bwdn', x, self.adapt[si])
        e = self.dirmlp(dirs.reshape(B * W, 3))
        h = x.reshape(B * W, C, N)
        for conv, norm, film in zip(self.convs, self.norms, self.films):
            h = norm(conv(h))
            g, b = film(e).chunk(2, -1)
            h = F.gelu(h * (1 + g[..., None]) + b[..., None])
        h = self.proj(h.mean(-1)).reshape(B, W, -1).permute(0, 2, 1)
        for blk in self.blocks:
            h = h + blk(h)
        return self.head(h).permute(0, 2, 1)


def bal_acc(p, t):
    hit = np.zeros(NCLS); tot = np.zeros(NCLS)
    for k in range(NCLS):
        m = t == k; tot[k] = m.sum(); hit[k] = (p[m] == k).sum()
    ok = tot > 0
    return float(np.mean(hit[ok] / tot[ok])) if ok.any() else float('nan')


def evaluate(net, dl, dev):
    net.eval(); P = []; T = []
    with torch.no_grad():
        for x, d, y, si in dl:
            p = net(x.to(dev), d.to(dev), si.to(dev)).argmax(-1).cpu().numpy().ravel()
            t = y.numpy().ravel(); m = t >= 0
            P.append(p[m]); T.append(t[m])
    return bal_acc(np.concatenate(P), np.concatenate(T))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--subs', type=int, nargs='*', default=list(KEEP))
    ap.add_argument('--window', type=int, default=256)
    ap.add_argument('--epochs', type=int, default=40)
    ap.add_argument('--batch', type=int, default=16)
    ap.add_argument('--lr', type=float, default=3e-4)
    ap.add_argument('--drop', type=float, default=0.3)
    ap.add_argument('--kind', default='filtered')
    ap.add_argument('--variant', default='ROI-PCA_8', choices=list(SRC))
    ap.add_argument('--task', default='gaze', choices=['gaze', 'motion'])
    a = ap.parse_args()
    global VARIANT, NV, NCLS
    VARIANT = a.variant; NV = SRC[VARIANT][2]
    NCLS = 4 if a.task == 'gaze' else 2
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    si = {s: i for i, s in enumerate(a.subs)}
    W, st = a.window, a.window // 2
    cut1, cut2 = int(0.6 * NR), int(0.8 * NR)
    mk = lambda lo, hi, stride: torch.utils.data.ConcatDataset(
        [Windows(s, lo, hi, W, stride, a.kind, si[s], a.task) for s in a.subs])
    tr, va, te = mk(0, cut1, st), mk(cut1, cut2, W), mk(cut2, NR, W)
    dl = lambda ds, sh: torch.utils.data.DataLoader(ds, a.batch, shuffle=sh, num_workers=4)

    net = Net(len(a.subs), W, drop=a.drop).to(dev)
    npar = sum(p.numel() for p in net.parameters())
    print(f'window {W} readouts ({W*0.008:.2f} s), TCN receptive field {net.rf} '
          f'({net.rf*0.008:.2f} s), {npar/1e3:.0f}k params')
    print(f'windows: train {len(tr)}, val {len(va)}, test {len(te)}, device {dev}')
    opt = torch.optim.AdamW(net.parameters(), lr=a.lr, weight_decay=1e-3)
    best = (-1, -1, None)
    for ep in range(a.epochs):
        net.train(); t0 = time.time(); s = n = 0
        for x, d, y, sj in dl(tr, True):
            loss = F.cross_entropy(net(x.to(dev), d.to(dev), sj.to(dev)).reshape(-1, NCLS),
                                   y.to(dev).reshape(-1), ignore_index=-1)
            opt.zero_grad(); loss.backward(); opt.step(); s += loss.item(); n += 1
        v = evaluate(net, dl(va, False), dev)
        if v > best[0]:
            best = (v, ep, {k: t.detach().clone() for k, t in net.state_dict().items()})
        print(f'ep {ep:3d} loss {s/n:.4f}  val {v:.4f}   [{time.time()-t0:.0f}s]', flush=True)
    net.load_state_dict(best[2])
    print(f'\nbest val {best[0]:.4f} at epoch {best[1]}')
    print(f'TEST balanced accuracy {evaluate(net, dl(te, False), dev):.4f}  (chance {1/NCLS:.2f})')
    print(f'linear probe reference: 0.28-0.30')


if __name__ == '__main__':
    main()
