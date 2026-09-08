"""Per-readout summary features, following Najdenovska et al., CBMS 2025.

Each k-space point is a complex number. Represent it four ways -- real,
imaginary, magnitude, phase -- and describe each readout by 23 statistics of
each representation. With nv = 8 virtual channels that is

    23 stats x 4 measures x 8 channels = 736 features per readout

against the paper's 23 x 4 x 20 = 1840 for its 20 physical channels.

The paper computed these on the full 384-sample readout; here the readout is
480 samples. Features are computed on the whole readout, not the |k| band used
elsewhere in model/, because these statistics are meant to summarise the entire
frequency sweep.
"""
import sys
import time
import numpy as np
import h5py
from pathlib import Path

STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/features')
NR, NS, NV = 79894, 480, 8

NAMES = ['min', 'max', 'mean', 'median', 'std', 'iqr', 'skew', 'kurt', 'ptp',
         'mean_abs_dev', 'mean_abs_diff', 'mean_diff', 'med_abs_dev',
         'med_abs_diff', 'med_diff', 'rms', 'autocorr', 'travelled',
         'entropy', 'zero_cross', 'pos_turn', 'neg_turn', 'ecdf_slope']
MEASURES = ['re', 'im', 'mag', 'phase']


def stats(x):
    """23 statistics along the last axis. x [..., NS] -> [..., 23]"""
    d = np.diff(x, axis=-1)
    mu = x.mean(-1, keepdims=True)
    med = np.median(x, -1, keepdims=True)
    sd = x.std(-1)
    q1, q3 = np.percentile(x, [25, 75], axis=-1)
    p10, p90 = np.percentile(x, [10, 90], axis=-1)
    z = (x - mu) / (sd[..., None] + 1e-12)
    sgn = np.sign(x)
    out = [
        x.min(-1), x.max(-1), mu[..., 0], med[..., 0], sd, q3 - q1,
        (z ** 3).mean(-1), (z ** 4).mean(-1) - 3.0, np.ptp(x, -1),
        np.abs(x - mu).mean(-1), np.abs(d).mean(-1), d.mean(-1),
        np.median(np.abs(x - med), -1), np.median(np.abs(d), -1), np.median(d, -1),
        np.sqrt((x ** 2).mean(-1)),
        ((x[..., :-1] - mu) * (x[..., 1:] - mu)).mean(-1) / (sd ** 2 + 1e-12),
        np.abs(d).sum(-1),
    ]
    # Shannon entropy of a 32-bin histogram, per readout
    lo, hi = x.min(-1, keepdims=True), x.max(-1, keepdims=True)
    idx = np.clip(((x - lo) / (hi - lo + 1e-12) * 32).astype(np.int16), 0, 31)
    flat = idx.reshape(-1, x.shape[-1])
    cnt = np.zeros((flat.shape[0], 32), np.float32)
    np.add.at(cnt, (np.arange(flat.shape[0])[:, None], flat), 1.0)
    p = cnt / x.shape[-1]
    ent = -(p * np.log2(p + 1e-12)).sum(-1).reshape(x.shape[:-1])
    out.append(ent)
    out.append((sgn[..., :-1] != sgn[..., 1:]).sum(-1).astype(np.float32))
    dd = np.sign(d)
    out.append(((dd[..., :-1] > 0) & (dd[..., 1:] < 0)).sum(-1).astype(np.float32))
    out.append(((dd[..., :-1] < 0) & (dd[..., 1:] > 0)).sum(-1).astype(np.float32))
    out.append(0.8 / (p90 - p10 + 1e-12))            # ECDF slope over the 10-90 range
    return np.stack(out, -1)


def feature_names():
    return [f'{m}_{s}_ch{c}' for c in range(NV) for m in MEASURES for s in NAMES]


def run(sub, block=2000):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / f'sub-{sub:03d}_readout_features.npy'
    if p.exists():
        print(f'sub-{sub:03d} already done'); return
    src = STUDY / f'sub-{sub:03d}/recon/ROI-PCA/kspace_ROI-PCA_8_woBin.mat'
    F = np.empty((NR, NV * 4 * 23), np.float32)
    t0 = time.time()
    with h5py.File(src, 'r') as f:
        y = f['y']
        for r0 in range(0, NR, block):
            r1 = min(r0 + block, NR)
            b = y[:, r0 * NS:r1 * NS]
            b = (b['real'] + 1j * b['imag']).reshape(NV, r1 - r0, NS)
            cols = []
            for c in range(NV):
                z = b[c]
                for m in (z.real, z.imag, np.abs(z), np.angle(z)):
                    cols.append(stats(m.astype(np.float64)).astype(np.float32))
            F[r0:r1] = np.concatenate(cols, -1)
    np.save(p, F)
    print(f'sub-{sub:03d}  {F.shape}  {F.nbytes/2**20:.0f} MiB  {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
    for s in ([int(a) for a in sys.argv[1:]] or list(KEEP)):
        run(s)
