"""Does the gaze evidence point the way the anatomy says it should?

A radial spoke measures the Fourier transform of the 1D projection of the object
onto that spoke's direction. A projection onto direction d is sensitive only to
the component of a displacement that lies *along* d -- displacement
perpendicular to d leaves the projection unchanged.

The lens and the optic nerve head are what move. So if the decoded evidence is
really gaze, discriminability must be anisotropic in a specific way:

    left/right evidence  concentrated in spokes aligned with the L-R axis
    up/down evidence     concentrated in spokes aligned with the S-I axis

and those two axes must be *different and orthogonal*. If the evidence is
isotropic, or if both contrasts prefer the same axis, then whatever is being
decoded is not eye position.

The test does not assume which array axis is which -- it measures the dependence
on all three and checks whether the two contrasts pick different ones.

Clusters are shared across subjects (the trajectory is bit-identical), so the
same cluster can be averaged over all 11 subjects for a sqrt(11) gain.
"""
import glob
import sys
import time
import numpy as np
import scipy.io as sio
from pathlib import Path
from dprime import prep, CIDX, OUT, NR, NC

NFOLD = 5
NSEG, NSHOT, NOFF = 44, 1872, 14
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
SHIFTS = (11000, 23000, 37000, 51000)
DIRS = np.load(OUT / 'dirs.npy')


def gaze(sub, kind='filtered'):
    out = []
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/{kind}/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        out.append(m[NOFF:, 1:].ravel().astype(bool))
    return np.stack(out)


def auc(s, a, b):
    x = np.concatenate([s[a], s[b]])
    o = np.argsort(x, kind='mergesort'); sx = x[o]
    r = np.empty(len(x)); r[o] = np.arange(len(x), dtype=float)
    i = 0
    while i < len(sx):
        j = i
        while j + 1 < len(sx) and sx[j + 1] == sx[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j) / 2
        i = j + 1
    n1, n0 = int(a.sum()), int(b.sum())
    return (r[:n1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)


def cluster_auc(Z, ga, gb):
    """Per-cluster discriminability of a two-class contrast, time-blocked CV."""
    fold = np.arange(NR) * NFOLD // NR
    s = np.full(NR, np.nan, np.float32)
    for f in range(NFOLD):
        te = fold == f
        for idx in CIDX:
            tr_i, te_i = idx[~te[idx]], idx[te[idx]]
            if te_i.size == 0 or tr_i.size < 40:
                continue
            A, B = tr_i[ga[tr_i]], tr_i[gb[tr_i]]
            if A.size < 8 or B.size < 8:
                continue
            T = Z[A].mean(0) - Z[B].mean(0)
            Zt = Z[te_i] - Z[te_i].mean(0)
            s[te_i] = (Zt @ np.conj(T)).real
    out = np.full(NC, np.nan)
    for c, idx in enumerate(CIDX):
        ok = idx[~np.isnan(s[idx])]
        a, b = ok[ga[ok]], ok[gb[ok]]
        if a.size >= 8 and b.size >= 8:
            m = np.zeros(NR, bool); m[a] = True
            n = np.zeros(NR, bool); n[b] = True
            out[c] = auc(s, m, n) - 0.5
    return out


def main():
    subs = [int(x) for x in sys.argv[1:]] or list(KEEP)
    cen = np.array([DIRS[idx].mean(0) for idx in CIDX])
    cen /= np.linalg.norm(cen, axis=1, keepdims=True)
    contrasts = {'left/right': (2, 3), 'up/down': (0, 1)}
    acc = {k: {'true': [], 'null': []} for k in contrasts}
    for sub in subs:
        t0 = time.time()
        Z = prep(sub)
        g = gaze(sub)
        for nm, (i, j) in contrasts.items():
            acc[nm]['true'].append(cluster_auc(Z, g[i], g[j]))
            for sh in SHIFTS:
                gs = np.roll(g, sh, axis=1)
                acc[nm]['null'].append(cluster_auc(Z, gs[i], gs[j]))
        print(f'  sub-{sub:03d} done [{time.time()-t0:.0f}s]', flush=True)

    print(f'\nper-cluster discriminability (AUC-0.5), averaged over {len(subs)} subjects')
    print('binned by |component of the spoke direction| along each array axis\n')
    for nm in contrasts:
        T = np.nanmean(np.stack(acc[nm]['true']), 0)
        Nl = np.stack(acc[nm]['null'])
        nmu, nsd = np.nanmean(Nl, 0), np.nanstd(Nl, 0)
        print(f'  {nm}:')
        for ax, axname in enumerate(['axis0', 'axis1', 'axis2']):
            comp = np.abs(cen[:, ax])
            edges = np.quantile(comp[~np.isnan(T)], [0, .25, .5, .75, 1.0])
            cells = []
            for k in range(4):
                m = (comp >= edges[k]) & (comp <= edges[k + 1]) & ~np.isnan(T)
                z = np.nanmean(T[m] - nmu[m]) / (np.nanmean(nsd[m]) / np.sqrt(m.sum()))
                cells.append(f'{np.nanmean(T[m])*100:+5.2f} (z{z:+5.1f})')
            print(f'    |d_{axname}| quartiles: ' + '  '.join(cells))
        print()
    print('  reading: entries are 100*(AUC-0.5) with a z against the shifted-label null.')
    print('  a real gaze effect rises across the quartiles of ONE axis, and')
    print('  left/right and up/down must prefer DIFFERENT axes.')


if __name__ == '__main__':
    main()
