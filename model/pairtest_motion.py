"""Is there a per-readout MOTION signal? The same bound, on the other task.

The gaze pair test found no per-readout difference between gaze bins, and that
result does NOT carry over to motion. Gaze position needs the eye's *location*
resolved -- a sub-millimetre structural difference. A saccade or blink instead
corrupts the phase of the readout being acquired *while it happens*, whatever
direction the eye points. Different mechanism, so it needs its own measurement.

Method mirrors pairtest_raw.py so the numbers are comparable: take readout pairs
with near-identical spoke direction, and compare

    still x motion   pairs   (one readout corrupted)
    still x still    pairs   (baseline: noise, drift, angular mismatch)

If motion leaves a per-readout signature, the first group differs more.

Uses the `_sync` labels: winLen=3 (24 ms window, saccade retention 70% vs 38%)
with the per-subject timing correction applied.
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
TAG = 'winLen3_sync'


def _mask(sub, name):
    p = f'{STUDY}/sub-{sub:03d}/recon/bins/{name}/eMask_th0.75_{TAG}.mat'
    m = np.asarray(sio.loadmat(p)['eMaskN']).squeeze().reshape(NSHOT, NSEG)
    return m[NOFF:, 1:].ravel().astype(bool)


def motion_labels(sub):
    """0 = still, 1 = moving, -1 = transitional (excluded)."""
    mot = np.zeros(NR, bool)
    for n in ('saccade-only', 'blink-event', 'tracking-loss'):
        try:
            mot |= _mask(sub, n)
        except Exception:
            pass
    fix = _mask(sub, 'fixation-ok')
    lab = np.full(NR, -1, np.int8)
    lab[fix & ~mot] = 0
    lab[mot] = 1
    return lab


def pair_d2(Y, a, b):
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
    lab = motion_labels(sub)
    D = np.load(KB / 'dirs.npy')
    pairs = cKDTree(D).query_pairs(r=2 * np.sin(np.radians(MAXDEG) / 2),
                                   output_type='ndarray')
    a, b = pairs[:, 0], pairs[:, 1]
    la, lb = lab[a], lab[b]
    ss = (la == 0) & (lb == 0)                       # still x still
    sm = ((la == 0) & (lb == 1)) | ((la == 1) & (lb == 0))   # still x motion
    if sm.sum() < 200:
        return None
    d2ss = pair_d2(Y, a[ss], b[ss])
    d2sm = pair_d2(Y, a[sm], b[sm])
    gap_ss = np.median(np.abs(a[ss] - b[ss])) * 0.008
    gap_sm = np.median(np.abs(a[sm] - b[sm])) * 0.008
    rng = np.random.default_rng(0)
    r = []
    for _ in range(300):
        i = rng.integers(0, len(d2ss), len(d2ss))
        j = rng.integers(0, len(d2sm), len(d2sm))
        r.append(np.median(d2sm[j]) / np.median(d2ss[i]))
    lo, hi = np.percentile(r, [2.5, 97.5])
    return dict(nss=int(ss.sum()), nsm=int(sm.sum()),
                ratio=np.median(d2sm) / np.median(d2ss), lo=lo, hi=hi,
                gap_ss=gap_ss, gap_sm=gap_sm)


if __name__ == '__main__':
    subs = [int(x) for x in sys.argv[1:]] or [15, 13, 10, 1]
    print(f'per-readout MOTION excess, labels = {TAG}')
    print('still x motion pairs vs still x still pairs, matched spoke direction\n')
    print(f'  {"sub":8s} {"repr":14s} {"n s-s":>8s} {"n s-m":>7s} {"ratio":>8s} '
          f'{"95% CI":>18s} {"gap s-s/s-m":>13s}  verdict')
    for s in subs:
        for which, name in (('roipca8', 'ROI-PCA nv=8'), ('raw52', 'uncompressed')):
            try:
                r = run(s, which)
            except FileNotFoundError:
                continue
            if r is None:
                print(f'  sub-{s:03d} {name:14s}  too few motion pairs'); continue
            v = 'EXCESS' if r['lo'] > 1.0 else 'no excess'
            print(f'  sub-{s:03d} {name:14s} {r["nss"]:8d} {r["nsm"]:7d} '
                  f'{r["ratio"]:8.4f}   [{r["lo"]:.4f}, {r["hi"]:.4f}] '
                  f'{r["gap_ss"]:5.0f}/{r["gap_sm"]:5.0f}s  {v}', flush=True)
    print('\n  ratio > 1 with the CI clear of 1 means a readout acquired during')
    print('  motion differs from a still one -- i.e. training is worth the run.')
