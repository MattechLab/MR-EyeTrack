"""Is there ANY per-readout gaze signal, before compression?

empirical_dprime.py ran this on ROI-PCA nv=8 and found no excess. Compression is
a linear map, so it can only lose information -- but that means the compressed
result cannot by itself prove the information was never there. Running the same
test on the uncompressed 52 channels bounds what *any* channel recombination
could ever deliver.

Method (identical to empirical_dprime.py, so the two numbers are comparable):
take readout pairs pointing the same way (the trajectory revisits directions,
~80 s apart, hence usually in different gaze states), and compare
||y_a - y_b||^2 for same-bin against cross-bin pairs. Same-bin pairs carry noise
and drift; cross-bin pairs carry that plus any gaze difference.

The 52-channel band is 5 GiB per subject, so pair distances are accumulated in
chunks rather than materialising the gathered arrays.
"""
import glob
import sys

import numpy as np
import scipy.io as sio
from pathlib import Path
from scipy.spatial import cKDTree

KB = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
RAW = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband_raw52')
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NR, NSEG, NSHOT, NOFF = 79894, 44, 1872, 14
MAXDEG = 1.0
CHUNK = 4000


def gaze(sub):
    lab = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/clean/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        lab[m[NOFF:, 1:].ravel().astype(bool)] = r
    return lab


def pair_d2(Y, a, b):
    """||Y[a]-Y[b]||^2 in chunks -- the gathered arrays would be tens of GiB."""
    out = np.empty(len(a), np.float64)
    for i in range(0, len(a), CHUNK):
        j = slice(i, min(i + CHUNK, len(a)))
        d = Y[a[j]] - Y[b[j]]
        out[j] = np.einsum('ij,ij->i', d, np.conj(d)).real
    return out


def run(sub, which):
    p = (RAW / f'sub-{sub:03d}_band160-320_raw52.npy' if which == 'raw52'
         else KB / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy')
    Y = np.load(p, mmap_mode='r').reshape(NR, -1)
    D = np.load(KB / 'dirs.npy')
    lab = gaze(sub)
    pairs = cKDTree(D).query_pairs(r=2 * np.sin(np.radians(MAXDEG) / 2),
                                   output_type='ndarray')
    a, b = pairs[:, 0], pairs[:, 1]
    ok = (lab[a] >= 0) & (lab[b] >= 0)
    a, b = a[ok], b[ok]
    same = lab[a] == lab[b]
    d2 = pair_d2(Y, a, b)
    ms, mc = np.median(d2[same]), np.median(d2[~same])
    rng = np.random.default_rng(0)
    r = []
    for _ in range(300):
        i = rng.integers(0, same.sum(), same.sum())
        j = rng.integers(0, (~same).sum(), (~same).sum())
        r.append(np.median(d2[~same][j]) / np.median(d2[same][i]))
    lo, hi = np.percentile(r, [2.5, 97.5])
    return dict(nsame=int(same.sum()), ncross=int((~same).sum()),
                ratio=mc / ms, lo=lo, hi=hi, chans=Y.shape[1])


if __name__ == '__main__':
    subs = [int(x) for x in sys.argv[1:]] or [15, 13, 10]
    print('per-readout gaze excess: cross-bin vs same-bin readout pairs')
    print('ratio 1.000 = no gaze difference; >1 = cross-bin pairs differ more\n')
    print(f'  {"sub":8s} {"representation":16s} {"features":>9s} {"ratio":>8s} '
          f'{"95% CI":>18s}  verdict')
    for s in subs:
        for which, name in (('roipca8', 'ROI-PCA nv=8'), ('raw52', 'uncompressed 52')):
            try:
                r = run(s, which)
            except FileNotFoundError:
                print(f'  sub-{s:03d} {name:16s}  missing'); continue
            v = 'excess present' if r['lo'] > 1.0 else 'no excess'
            print(f'  sub-{s:03d} {name:16s} {r["chans"]:9d} {r["ratio"]:8.4f} '
                  f'  [{r["lo"]:.4f}, {r["hi"]:.4f}]  {v}', flush=True)
    print('\n  If the uncompressed row also shows no excess, no choice of ROI mask')
    print('  or channel recombination can help: the information is not in the readouts.')
