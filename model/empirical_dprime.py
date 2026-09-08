"""Per-readout discriminability, measured directly. No forward model.

The ceiling calculation keeps returning d' of 8-25 per readout, which would mean
one readout classifies gaze perfectly -- flatly contradicted by every decode we
have run. Rather than trust a simulation, measure the same quantity from the
acquired data.

A radial spoke's signal depends on its direction, so two readouts are only
comparable if they point the same way. The trajectory revisits directions
(~12 spokes within 1 degree of any given one), and those revisits land far apart
in time and therefore often in different gaze states. So:

    take pairs of readouts with near-identical spoke direction
    split them into same-gaze-bin and different-gaze-bin pairs
    compare  ||y_a - y_b||^2  between the two groups

Same-bin pairs carry noise, drift and everything that is not gaze. Cross-bin
pairs carry all of that plus the gaze difference. The excess is the gaze signal,
and its ratio to the baseline gives the per-readout d' -- with no forward model,
no coil maps and no NUFFT.
"""
import glob

import numpy as np
import scipy.io as sio
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NR, NSEG, NSHOT, NOFF = 79894, 44, 1872, 14
# The trajectory's median nearest-neighbour spacing is 0.71 deg, so a tolerance
# below that finds nothing. At 1.0 deg there are ~225k pairs. The residual
# angular mismatch (~1.4 samples of 240 at the band's outer |k|) affects same-bin
# and cross-bin pairs identically, so the ratio between them stays meaningful.
MAXDEG = 1.0


def gaze(sub):
    lab = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/clean/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        lab[m[NOFF:, 1:].ravel().astype(bool)] = r
    return lab


def main(sub=15):
    Y = np.load(OUT / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy').reshape(NR, -1)
    D = np.load(OUT / 'dirs.npy')
    lab = gaze(sub)

    tree = cKDTree(D)
    pairs = tree.query_pairs(r=2 * np.sin(np.radians(MAXDEG) / 2), output_type='ndarray')
    a, b = pairs[:, 0], pairs[:, 1]
    ok = (lab[a] >= 0) & (lab[b] >= 0)
    a, b = a[ok], b[ok]
    same = lab[a] == lab[b]
    print(f'sub-{sub:03d}: {len(a)} labelled readout pairs within {MAXDEG} deg')
    print(f'  {int(same.sum())} same-bin, {int((~same).sum())} cross-bin; '
          f'median time gap {np.median(np.abs(a - b) * 0.008):.0f} s')

    d2 = (np.abs(Y[a] - Y[b]) ** 2).sum(1)
    ms, mc = np.median(d2[same]), np.median(d2[~same])
    print(f'\n  median ||y_a - y_b||^2   same-bin {ms:.5g}   cross-bin {mc:.5g}')
    print(f'  ratio cross/same {mc / ms:.4f}')
    excess = mc - ms
    print(f'  excess attributable to gaze: {excess:+.4g} '
          f'({100 * excess / ms:+.2f}% of the same-bin baseline)')
    if excess > 0:
        print(f"  -> per-readout d' = {np.sqrt(excess / ms):.3f}")
    else:
        print('  -> no excess: cross-bin pairs are no more different than same-bin pairs')

    print(f'\n  {"contrast":14s} {"n pairs":>8s} {"ratio to same-bin":>18s}')
    for nm, (i, j) in {'up/down': (0, 1), 'left/right': (2, 3),
                       'up/left': (0, 2), 'down/right': (1, 3)}.items():
        m = ((lab[a] == i) & (lab[b] == j)) | ((lab[a] == j) & (lab[b] == i))
        if m.sum() > 50:
            print(f'  {nm:14s} {int(m.sum()):8d} {np.median(d2[m]) / ms:18.4f}')


if __name__ == '__main__':
    main()
