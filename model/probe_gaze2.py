"""Positive control, v2. Two fixes over probe_gaze.py, which failed its null.

v1 decoded *time*, not gaze: a circular-shift null scored 0.33 against 0.25
chance, because slow drift makes a test block resemble whichever training
blocks sit near it, and any temporally blocky label inherits that.

  fix 1  per-cluster local detrend. A cluster's ~133 readouts are spread evenly
         over the scan (~5 s apart, ~one per gaze block), so subtracting a
         moving average over 11 cluster-neighbours removes anything slower than
         ~55 s while leaving the gaze-block rate intact.
  fix 2  reduce to K principal components before building templates. v1 fit a
         2560-real-dimensional template from ~27 readouts per class per fold,
         which is essentially all noise.

The circular-shift null is the acceptance test: real signal must beat it.
"""
import sys
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
NC, NFOLD, MAVG, KPC = 600, 5, 11, 32
WINDOWS = (32, 128, 256, 512, 1024)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)

dirs = np.load(OUT / 'dirs.npy')
NR = dirs.shape[0]
CL = cKDTree(dirs[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(dirs)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]


def detrend_and_reduce(X):
    """Per-cluster local detrend, then project onto K global principal components."""
    R = np.empty_like(X)
    for idx in CIDX:                                   # already in time order
        blk = X[idx]
        n = len(idx)
        w = min(MAVG, n if n % 2 else n - 1)
        if w < 3:
            R[idx] = blk - blk.mean(0); continue
        pad = np.concatenate([blk[:w//2][::-1], blk, blk[-(w//2):][::-1]])
        ker = np.ones(w) / w
        loc = np.apply_along_axis(lambda v: np.convolve(v, ker, 'valid'), 0, pad)
        R[idx] = blk - loc
    sub = R[np.random.default_rng(1).choice(NR, 12000, replace=False)]
    _, _, Vh = np.linalg.svd(sub, full_matrices=False)
    P = Vh[:KPC].conj().T
    return R @ P                                        # [NR, KPC] complex


def readout_scores(Z, gaze, shift=0):
    g = np.roll(gaze, shift, axis=1) if shift else gaze
    fold = np.arange(NR) * NFOLD // NR
    S = np.full((NR, 4), np.nan, np.float32)
    for f in range(NFOLD):
        te = fold == f
        for c in range(NC):
            idx = CIDX[c]
            tr_i, te_i = idx[~te[idx]], idx[te[idx]]
            if te_i.size == 0 or tr_i.size < 40:
                continue
            mu = Z[tr_i].mean(0)
            T, keep = [None]*4, []
            for k in range(4):
                sel = tr_i[g[k][tr_i]]
                if sel.size >= 5:
                    T[k] = Z[sel].mean(0) - mu; keep.append(k)
            if len(keep) < 2:
                continue
            Zt = Z[te_i] - Z[te_i].mean(0)
            for k in keep:
                S[te_i, k] = (Zt @ np.conj(T[k])).real - 0.5 * np.vdot(T[k], T[k]).real
            kk = np.array(keep)
            blk = S[te_i][:, kk]
            S[np.ix_(te_i, kk)] = (blk - blk.mean()) / (blk.std() + 1e-9)
    return S


def window_eval(S, gaze, W, stride=32, purity=0.6, min_lab=0.3):
    ok = ~np.isnan(S).any(1)
    cs = np.concatenate([np.zeros((1, 4)), np.cumsum(np.nan_to_num(S), 0)])
    cn = np.concatenate([[0], np.cumsum(ok)])
    lab = np.full(NR, -1, np.int8)
    for k in range(4):
        lab[gaze[k]] = k
    cl = np.concatenate([np.zeros((1, 4), np.int64),
                         np.cumsum(np.stack([(lab == k) for k in range(4)], 1).astype(np.int32), 0)])
    yh, yt = [], []
    for a in range(0, NR - W + 1, stride):
        b = a + W
        if cn[b] - cn[a] < 0.5 * W:
            continue
        counts = cl[b] - cl[a]; n = counts.sum()
        if n < min_lab * W or counts.max() < purity * n:
            continue
        yh.append(int(np.argmax(cs[b] - cs[a]))); yt.append(int(np.argmax(counts)))
    return np.array(yt), np.array(yh)


def run(sub, shift=0, _cache={}):
    if sub not in _cache:
        X = np.load(OUT / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy').reshape(NR, -1)
        _cache.clear()
        _cache[sub] = detrend_and_reduce(X / (np.abs(X).mean() + 1e-12))
    Z = _cache[sub]
    gaze = np.load(OUT / f'sub-{sub:03d}_labels.npz')['gaze']
    S = readout_scores(Z, gaze, shift)
    return {W: float((lambda t, h: (t == h).mean() if t.size else np.nan)(*window_eval(S, gaze, W)))
            for W in WINDOWS}


if __name__ == '__main__':
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    print('4-way gaze accuracy, within subject, time-blocked CV (chance 0.25)')
    print('window   ' + ' '.join(f'{W*0.008:5.1f}s' for W in WINDOWS) + '\n')
    real, null = {W: [] for W in WINDOWS}, {W: [] for W in WINDOWS}
    for s in subs:
        r = run(s); n = run(s, shift=20000)
        print(f'sub-{s:03d}  true ' + ' '.join(f'{r[W]:6.3f}' for W in WINDOWS) +
              '   | null ' + ' '.join(f'{n[W]:6.3f}' for W in WINDOWS), flush=True)
        for W in WINDOWS:
            real[W].append(r[W]); null[W].append(n[W])
    print('\nmedian   true ' + ' '.join(f'{np.nanmedian(real[W]):6.3f}' for W in WINDOWS) +
          '   | null ' + ' '.join(f'{np.nanmedian(null[W]):6.3f}' for W in WINDOWS))
