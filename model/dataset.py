"""Readout-level dataset for eye-state classification from compressed k-space.

The index arithmetic between the eye-tracking masks and the exported k-space is
the one thing here that fails silently, so it lives in a single place and is
asserted rather than assumed.

Acquisition layout (pulseq main sequence, all 15 subjects):

    nSeg = 44, nShot = 1872        -> 82 368 acquired readouts
    segment 1 of every shot        -> SI navigator, dropped by the reader
    first nShotOff = 14 shots      -> warm-up, dropped by the reader
    43 x 1858                      -> 79 894 readouts, the exported `y`

The ET masks in recon/bins/ are 82 368 long and in *acquisition* order, so the
same two drops map them onto the export. MATLAB reshapes column-major, which is
why the numpy reshape is (nShot, nSeg) and not the other way round.

Readout geometry: LIBRE spokes are full-diameter projections. |k| runs
1 -> 0 -> 1 across the 480 samples, so sample 240 is k = 0 and the sample index
is simultaneously a k-space coordinate and an echo time. That second role is
why antipodal spokes are *not* interchangeable, even though they cover the same
k-space line.
"""
from pathlib import Path
import numpy as np
import h5py
import scipy.io as sio

NSEG, NSHOT, NSHOT_OFF = 44, 1872, 14
NSAMP, NREADOUT = 480, 79894
DC_SAMPLE = 240
BASE = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')

CLASSES = ('fixation-ok', 'saccade-only', 'blink-event', 'tracking-loss')


def _drop_to_export(mask_82368):
    """Acquisition-order mask -> export-order mask (79 894)."""
    m = np.asarray(mask_82368).squeeze()
    assert m.size == NSEG * NSHOT, f'expected {NSEG*NSHOT} readouts, got {m.size}'
    m = m.reshape(NSHOT, NSEG)          # MATLAB [nSeg x nShot], column-major
    out = m[NSHOT_OFF:, 1:].ravel()     # drop warm-up shots, then the SI segment
    assert out.size == NREADOUT
    return out


def load_labels(sub, name):
    """Per-readout binary mask aligned to the exported k-space."""
    if name in ('fixation-ok', 'saccade-only', 'blink-event', 'tracking-loss', 'no-mo'):
        p = BASE / f'sub-{sub:03d}/recon/bins/{name}/eMask_th0.75_winLen10.mat'
    else:                                # gaze bins: 'clean:0' .. 'clean:3'
        mask_type, region = name.split(':')
        p = BASE / f'sub-{sub:03d}/recon/bins/{mask_type}/eMask_th0.75_region{region}.mat'
    return _drop_to_export(sio.loadmat(p)['eMaskN']).astype(bool)


def label_frame(sub):
    """All ET classes at once, plus the readouts no class claims.

    The classes are not a partition: `gap` holds readouts the 0.75 sliding
    window would not commit to, which sit mostly at event boundaries.
    """
    d = {c: load_labels(sub, c) for c in CLASSES}
    d['gap'] = ~np.logical_or.reduce(list(d.values()))
    return d


def load_band(sub, lo=176, hi=304, variant='ROI-PCA', nv=8, fname=None):
    """k-space samples [lo, hi) of every readout, plus the spoke directions.

    Returns
        y     [nReadout, hi-lo, nv] complex64
        dirs  [nReadout, 3]         float64, unit vectors
    """
    fname = fname or f'kspace_{variant}_{nv}_woBin.mat'
    p = BASE / f'sub-{sub:03d}/recon/{variant}/{fname}'
    with h5py.File(p, 'r') as f:
        ds = f['y']
        nch = ds.shape[0]
        y = np.empty((NREADOUT, hi - lo, nch), np.complex64)
        for r0 in range(0, NREADOUT, 2000):     # gzip chunks, so read in blocks
            r1 = min(r0 + 2000, NREADOUT)
            b = ds[:, r0 * NSAMP:r1 * NSAMP]
            b = (b['real'] + 1j * b['imag']).astype(np.complex64)
            y[r0:r1] = np.transpose(b.reshape(nch, r1 - r0, NSAMP)[:, :, lo:hi], (1, 2, 0))
        d = f['t'][:, 0, :].astype(np.float64)  # sample 0 sits at |k| = 1
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return y, d


def k_axis(lo=176, hi=304):
    """Normalised |k| for each retained sample. 1.0 == Nyquist for 1 mm voxels."""
    return np.abs(np.arange(lo, hi) - DC_SAMPLE) / float(DC_SAMPLE)


if __name__ == '__main__':
    for sub in (15, 5):
        L = label_frame(sub)
        print(f'sub-{sub:03d}: ' + '  '.join(f'{k} {v.mean():.3f}' for k, v in L.items()))
