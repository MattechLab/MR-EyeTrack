"""Per-readout features from the UNCOMPRESSED 52-channel raw data.

The compressed (ROI-PCA nv=8) version of these features gave chance on motion
detection both across and within subjects, while the CBMS 2025 paper reached
98% on uncompressed 20-channel data. The working hypothesis is that ROI-PCA
destroys the coil phase structure the method depends on. This extracts the same
23 statistics x {re, im, mag, phase} from all 52 physical coils:

    23 x 4 x 52 = 4784 features per readout   (1.5 GiB per subject)

Readout ordering matches monalisa exactly: the twix mdb list is in acquisition
order, reshaped (nShot=1872, nSeg=44) with segment fastest, then segment 0 (the
SI navigator) and the first 14 warm-up shots are dropped -> 43 x 1858 = 79 894,
the same readouts as every other dataset in model/.

Order statistics are computed from a single sort per readout rather than
independent percentile calls, which is what makes 52 channels tractable.
"""
import contextlib
import io
import sys
import time
from pathlib import Path

import numpy as np
import scipy.integrate as _si

_si.cumtrapz = getattr(_si, 'cumtrapz', _si.cumulative_trapezoid)   # twixtools needs the old name
sys.path.insert(0, '/home/debi/jaime/repos/MR-EyeTrack/old_study/code/twixtools')
import twixtools  # noqa: E402

STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/features_raw')
NSEG, NSHOT, NOFF = 44, 1872, 14
NS, NCH = 480, 52
NR = (NSHOT - NOFF) * (NSEG - 1)          # 79894
NSTAT = 23
NAMES = ['min', 'max', 'mean', 'median', 'std', 'iqr', 'skew', 'kurt', 'ptp',
         'mean_abs_dev', 'mean_abs_diff', 'mean_diff', 'med_abs_dev',
         'med_abs_diff', 'med_diff', 'rms', 'autocorr', 'travelled',
         'entropy', 'zero_cross', 'pos_turn', 'neg_turn', 'ecdf_slope']
MEASURES = ['re', 'im', 'mag', 'phase']


def stats(x):
    """23 statistics along the last axis. x [..., NS] float32 -> [..., 23]"""
    n = x.shape[-1]
    xs = np.sort(x, axis=-1)                      # one sort serves every order statistic
    d = np.diff(x, axis=-1)
    m2 = (n - 1) // 2
    ds = np.partition(np.abs(d), m2, axis=-1)     # median only: partition beats a full sort
    mu = x.mean(-1, keepdims=True)
    sd = x.std(-1)
    med = xs[..., n // 2]
    q1, q3 = xs[..., n // 4], xs[..., (3 * n) // 4]
    p10, p90 = xs[..., n // 10], xs[..., (9 * n) // 10]
    z = (x - mu) / (sd[..., None] + 1e-12)
    out = [xs[..., 0], xs[..., -1], mu[..., 0], med, sd, q3 - q1,
           (z ** 3).mean(-1), (z ** 4).mean(-1) - 3.0, xs[..., -1] - xs[..., 0],
           np.abs(x - mu).mean(-1), np.abs(d).mean(-1), d.mean(-1),
           np.partition(np.abs(x - med[..., None]), n // 2, axis=-1)[..., n // 2],
           ds[..., m2], np.partition(d, m2, axis=-1)[..., m2],
           np.sqrt((x ** 2).mean(-1)),
           ((x[..., :-1] - mu) * (x[..., 1:] - mu)).mean(-1) / (sd ** 2 + 1e-12),
           np.abs(d).sum(-1)]
    lo, hi = xs[..., :1], xs[..., -1:]
    idx = np.clip(((x - lo) / (hi - lo + 1e-12) * 32).astype(np.int8), 0, 31)
    flat = idx.reshape(-1, n).astype(np.int32)
    nrow = flat.shape[0]
    offs = flat + 32 * np.arange(nrow, dtype=np.int32)[:, None]
    cnt = np.bincount(offs.ravel(), minlength=32 * nrow).reshape(nrow, 32)
    p = cnt / n
    out.append((-(p * np.log2(p + 1e-12)).sum(-1)).reshape(x.shape[:-1]))
    sg = np.sign(x)
    out.append((sg[..., :-1] != sg[..., 1:]).sum(-1).astype(np.float32))
    dd = np.sign(d)
    out.append(((dd[..., :-1] > 0) & (dd[..., 1:] < 0)).sum(-1).astype(np.float32))
    out.append(((dd[..., :-1] < 0) & (dd[..., 1:] > 0)).sum(-1).astype(np.float32))
    out.append(0.8 / (p90 - p10 + 1e-12))
    return np.stack(out, -1).astype(np.float32)


def feature_names():
    return [f'{m}_{s}_ch{c}' for c in range(NCH) for m in MEASURES for s in NAMES]


def keep_index():
    """Acquisition-order indices of the readouts that survive to the export."""
    k = np.arange(NSHOT * NSEG).reshape(NSHOT, NSEG)
    return k[NOFF:, 1:].ravel()


def run(sub, block=1000):
    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / f'sub-{sub:03d}_raw_features.npy'
    if dst.exists():
        print(f'sub-{sub:03d} already done', flush=True); return
    src = next((STUDY / f'sub-{sub:03d}/rawdata').glob('*_T1wLIBRE.dat'))
    t0 = time.time()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        tw = twixtools.read_twix(str(src))
    img = [m for m in tw[-1]['mdb'] if m.is_image_scan()]
    assert len(img) == NSHOT * NSEG, f'{len(img)} image scans, expected {NSHOT*NSEG}'
    keep = keep_index()
    F = np.empty((NR, NCH * 4 * NSTAT), np.float32)
    for b0 in range(0, NR, block):
        b1 = min(b0 + block, NR)
        z = np.stack([img[i].data for i in keep[b0:b1]])          # [n, NCH, NS]
        cols = []
        for meas in (z.real, z.imag, np.abs(z), np.angle(z)):
            cols.append(stats(meas.astype(np.float32)))           # [n, NCH, 23]
        # interleave to channel-major order: ch0(re,im,mag,ph), ch1(...), ...
        S = np.concatenate([c[:, :, None, :] for c in cols], 2)   # [n, NCH, 4, 23]
        F[b0:b1] = S.reshape(b1 - b0, -1)
        if b0 % 20000 == 0:
            print(f'  sub-{sub:03d} {b0}/{NR}  {time.time()-t0:.0f}s', flush=True)
    np.save(dst, F)
    print(f'sub-{sub:03d}  {F.shape}  {F.nbytes/2**30:.2f} GiB  {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
    for s in ([int(a) for a in sys.argv[1:]] or list(KEEP)):
        run(s)
