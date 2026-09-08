"""Where in |k| does the decodable evidence live? Head or eye?

Bulk head motion is a whole-head effect and dominates LOW |k| (coarse spatial
scale). Eye structure -- the globe at 24 mm, the lens at ~9 mm -- lives at
|k| 0.08-0.33. If the decoding evidence sits entirely at low |k|, the decoder is
reading head pose, not eye position, which is what the isotropy result in
anisotropy.py already suggests.

Slicing commutes with the cached detrend (it is applied per sample and channel
independently), so the cache is reused directly.
"""
import glob
import sys
import time
import numpy as np
import scipy.io as sio
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NR, NC, NFOLD, KPC = 79894, 600, 5, 32
NSEG, NSHOT, NOFF = 44, 1872, 14
LO, HI = 160, 320
W = 256
SHIFTS = (11000, 23000, 37000, 51000)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
DIRS = np.load(OUT / 'dirs.npy')
CL = cKDTree(DIRS[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(DIRS)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]
KAX = np.abs(np.arange(LO, HI) - 240) / 240.0

BANDS = {'low |k| <0.06  (>16 mm, bulk head)':      KAX < 0.06,
         'mid |k| .06-.15 (7-16 mm, globe)':        (KAX >= 0.06) & (KAX < 0.15),
         'high |k| .15-.33 (3-7 mm, lens/ONH)':     KAX >= 0.15}


def gaze(sub):
    out = []
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/filtered/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        out.append(m[NOFF:, 1:].ravel().astype(bool))
    return np.stack(out)


def reduce_band(R, sel):
    X = R[:, sel, :].reshape(NR, -1)
    s = X[np.random.default_rng(1).choice(NR, 12000, replace=False)]
    _, _, Vh = np.linalg.svd(s, full_matrices=False)
    return X @ Vh[:KPC].conj().T


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
    print(f'4-way gaze, {W} readout windows ({W*0.008:.1f} s), '
          f'{len(subs)} subjects, {len(SHIFTS)}-shift null\n')
    res = {b: {'t': [], 'z': []} for b in BANDS}
    for sub in subs:
        t0 = time.time()
        R = np.load(OUT / 'detrended' / f'sub-{sub:03d}.npy', mmap_mode='r')
        R = np.asarray(R)
        g = gaze(sub)
        line = []
        for bname, sel in BANDS.items():
            Z = reduce_band(R, sel)
            S, sn = scores(Z, g); t = bal(S, sn, g)
            nul = []
            for sh in SHIFTS:
                gs = np.roll(g, sh, axis=1)
                Ss, sns = scores(Z, gs); nul.append(bal(Ss, sns, gs))
            nul = np.array(nul, float)
            z = (t - np.nanmean(nul)) / (np.nanstd(nul, ddof=1) + 1e-9)
            res[bname]['t'].append(t); res[bname]['z'].append(z)
            line.append(f'{t:.3f}(z{z:+.1f})')
        print(f'  sub-{sub:03d}  ' + '  '.join(line) + f'   [{time.time()-t0:.0f}s]', flush=True)
    print(f'\n{"band":38s} {"median acc":>11s} {"median z":>9s} {"subjects z>0":>13s}')
    for b in BANDS:
        t = np.array(res[b]['t']); z = np.array(res[b]['z'])
        print(f'{b:38s} {np.nanmedian(t):11.4f} {np.nanmedian(z):9.2f} '
              f'{int((z>0).sum()):8d}/{len(z)}')
    print('\nchance 0.25. If low |k| carries it and high |k| does not, the decoder')
    print('is reading bulk head pose rather than the eye.')


if __name__ == '__main__':
    main()
