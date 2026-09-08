"""Direct measurement: how separable are the gaze classes, per readout?

No classifier, no cross-validation, no class priors -- all three were what made
probe_gaze v1/v2 uninterpretable. Just a two-sample separation with its
estimation bias subtracted in closed form.

Within a direction cluster, for classes A and B in a D-real-dimensional space,

    E[ ||mu_A - mu_B||^2 ]  =  ||true separation||^2 + D*sigma^2*(1/nA + 1/nB)

so subtracting the second term gives an unbiased estimate of the real signal,
which can legitimately come out negative when the signal is zero. Pooling that
numerator over clusters gives per-readout d'^2. Shifted labels are the control:
the estimator must return ~0 there.
"""
import sys
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
NC, MAVG, KPC = 600, 11, 32
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
dirs = np.load(OUT / 'dirs.npy'); NR = dirs.shape[0]
CL = cKDTree(dirs[np.random.default_rng(0).choice(NR, NC, replace=False)]).query(dirs)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]


def prep(sub, detrend=True, kpc=KPC):
    X = np.load(OUT / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy').reshape(NR, -1)
    X = X / (np.abs(X).mean() + 1e-12)
    R = np.empty_like(X)
    for idx in CIDX:
        blk = X[idx]
        if detrend and len(idx) >= MAVG:
            w = MAVG
            pad = np.concatenate([blk[:w//2][::-1], blk, blk[-(w//2):][::-1]])
            loc = np.apply_along_axis(lambda v: np.convolve(v, np.ones(w)/w, 'valid'), 0, pad)
            R[idx] = blk - loc
        else:
            R[idx] = blk - blk.mean(0)
    if kpc:
        sub_ = R[np.random.default_rng(1).choice(NR, 12000, replace=False)]
        _, _, Vh = np.linalg.svd(sub_, full_matrices=False)
        R = R @ Vh[:kpc].conj().T
    return np.column_stack([R.real, R.imag]).astype(np.float64)   # D real dims


def dprime(Z, a_mask, b_mask, shift=0):
    if shift:
        a_mask, b_mask = np.roll(a_mask, shift), np.roll(b_mask, shift)
    D = Z.shape[1]
    num = den = 0.0
    used = 0
    for idx in CIDX:
        A = idx[a_mask[idx]]; B = idx[b_mask[idx]]
        if A.size < 8 or B.size < 8:
            continue
        za, zb = Z[A], Z[B]
        s2 = (za.var(0, ddof=1).sum() * (za.shape[0]-1) +
              zb.var(0, ddof=1).sum() * (zb.shape[0]-1)) / ((za.shape[0]+zb.shape[0]-2) * D)
        d2 = ((za.mean(0) - zb.mean(0))**2).sum() - D * s2 * (1/za.shape[0] + 1/zb.shape[0])
        num += d2; den += s2; used += 1
    # pooled over clusters: sum of unbiased separations / sum of variances
    return (num / den) if used else np.nan, used


if __name__ == '__main__':
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    print("per-readout d'^2 for gaze contrasts (unbiased; ~0 means no separation)")
    print(f'{"sub":8s} {"up/down":>10s} {"left/right":>11s} | '
          f'{"null u/d":>9s} {"null l/r":>9s}   clusters')
    for s in subs:
        Z = prep(s)
        g = np.load(OUT / f'sub-{s:03d}_labels.npz')['gaze']
        ud, n1 = dprime(Z, g[0], g[1])
        lr, _ = dprime(Z, g[2], g[3])
        udn, _ = dprime(Z, g[0], g[1], shift=20000)
        lrn, _ = dprime(Z, g[2], g[3], shift=20000)
        print(f'sub-{s:03d} {ud:10.4f} {lr:11.4f} | {udn:9.4f} {lrn:9.4f}   {n1}', flush=True)
