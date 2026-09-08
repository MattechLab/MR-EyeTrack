"""Head-to-head: does the choice of coil compression change what can be decoded?

Compares ROI-PCA (eye-ROI energy), SVD (whole-FoV energy) at nv=8 and 16 on the
one task that has ever risen above chance here -- 4-way gaze from windows of
readouts. Everything downstream of the channel map is identical: same per-cluster
local detrend, same 32 PCs, same direction-conditioned templates, same
time-blocked CV, same circular-shift null, same balanced accuracy.

ROVir is not re-run: it exists for one subject only, and the earlier head-to-head
already showed it matching ROI-PCA at 2.5x the channels.
"""
import glob
import sys
import time
import numpy as np
import scipy.io as sio
from pathlib import Path
from scipy.spatial import cKDTree

KB = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
SVD = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/svd')
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NR, NC, NFOLD, KPC, MAVG = 79894, 600, 5, 32, 11
NSEG, NSHOT, NOFF = 44, 1872, 14
W = 256
SHIFTS = (11000, 23000, 37000, 51000)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
DIRS = np.load(KB / 'dirs.npy')
CL = cKDTree(DIRS[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(DIRS)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]

VARIANTS = {
    'ROI-PCA nv=8': lambda s: KB / f'sub-{s:03d}_band160-320_ROI-PCA_8.npy',
    'SVD nv=8':     lambda s: SVD / f'sub-{s:03d}_band160-320_SVD_8.npy',
    'SVD nv=16':    lambda s: SVD / f'sub-{s:03d}_band160-320_SVD_16.npy',
}


def gaze(sub):
    out = []
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/filtered/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        out.append(m[NOFF:, 1:].ravel().astype(bool))
    return np.stack(out)


def prep(path):
    X = np.load(path).reshape(NR, -1)
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
    s = R[np.random.default_rng(1).choice(NR, 12000, replace=False)]
    _, _, Vh = np.linalg.svd(s, full_matrices=False)
    return R @ Vh[:KPC].conj().T


def scores(Z, g):
    fold = np.arange(NR) * NFOLD // NR
    S = np.zeros((NR, 4), np.float32); seen = np.zeros(NR, bool)
    for f in range(NFOLD):
        te = fold == f
        for idx in CIDX:
            tr_i, te_i = idx[~te[idx]], idx[te[idx]]
            if te_i.size == 0 or tr_i.size < 40:
                continue
            mu = Z[tr_i].mean(0); T = [None] * 4; keep = []
            for k in range(4):
                sel = tr_i[g[k][tr_i]]
                if sel.size >= 5:
                    T[k] = Z[sel].mean(0) - mu; keep.append(k)
            if len(keep) < 2:
                continue
            Zt = Z[te_i] - Z[te_i].mean(0)
            sc = np.stack([(Zt @ np.conj(T[k])).real - .5 * np.vdot(T[k], T[k]).real
                           for k in keep], 1)
            sc = (sc - sc.mean()) / (sc.std() + 1e-9)
            S[np.ix_(te_i, np.array(keep))] = sc; seen[te_i] = True
    return S, seen


def bal(S, seen, g):
    lab = np.full(NR, -1, np.int8)
    for k in range(4):
        lab[g[k]] = k
    oh = np.stack([(lab == k) for k in range(4)], 1).astype(np.int32)
    cs = np.concatenate([np.zeros((1, 4)), np.cumsum(np.where(seen[:, None], S, 0), 0)])
    cl = np.concatenate([np.zeros((1, 4), np.int64), np.cumsum(oh, 0)])
    cn = np.concatenate([[0], np.cumsum(seen)])
    hit = np.zeros(4); tot = np.zeros(4)
    for a in range(0, NR - W + 1, W // 2):
        b = a + W
        if cn[b] - cn[a] < 0.5 * W:
            continue
        c = cl[b] - cl[a]; n = c.sum()
        if n < 0.3 * W or c.max() < 0.7 * n:
            continue
        yt = int(np.argmax(c)); tot[yt] += 1
        hit[yt] += (yt == int(np.argmax(cs[b] - cs[a])))
    ok = tot > 0
    return float(np.mean(hit[ok] / tot[ok])) if ok.sum() >= 3 else np.nan


def main():
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    print(f'4-way gaze, {W}-readout windows ({W*0.008:.1f} s), balanced accuracy '
          f'(chance 0.250), {len(SHIFTS)}-shift null\n')
    print(f'  {"sub":8s} ' + '  '.join(f'{v:>17s}' for v in VARIANTS))
    res = {v: {'a': [], 'z': []} for v in VARIANTS}
    for sub in subs:
        g = gaze(sub); row = []
        for v, pathf in VARIANTS.items():
            p = pathf(sub)
            if not p.exists():
                row.append('     missing     '); continue
            Z = prep(p)
            S, sn = scores(Z, g); a = bal(S, sn, g)
            nul = []
            for sh in SHIFTS:
                gs = np.roll(g, sh, axis=1)
                Ss, sns = scores(Z, gs); nul.append(bal(Ss, sns, gs))
            nul = np.array(nul, float)
            z = (a - np.nanmean(nul)) / (np.nanstd(nul, ddof=1) + 1e-9)
            res[v]['a'].append(a); res[v]['z'].append(z)
            row.append(f'{a:.3f} (z{z:+5.1f})')
        print(f'  sub-{sub:03d} ' + '  '.join(f'{r:>17s}' for r in row), flush=True)
    print()
    print(f'  {"variant":18s} {"median acc":>11s} {"median z":>9s} {"subjects z>0":>13s}')
    for v in VARIANTS:
        a = np.array(res[v]['a']); z = np.array(res[v]['z'])
        if a.size:
            print(f'  {v:18s} {np.nanmedian(a):11.4f} {np.nanmedian(z):9.2f} '
                  f'{int((z>0).sum()):8d}/{len(z)}')


if __name__ == '__main__':
    main()
