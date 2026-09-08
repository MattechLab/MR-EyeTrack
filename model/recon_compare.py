"""Gridded reconstruction from ROI-PCA vs SVD compressed k-space, side by side.

Both are 52 -> 8 channel maps, so the compression is *identical* (6.5x); they
differ only in which 8-dimensional subspace they keep. This reconstructs each
the same way -- adjoint NUFFT per virtual channel with density compensation,
then root-sum-of-squares across channels -- so any visible difference is the
subspace, not the recon.

RSS is used rather than a coil-sensitivity combination because SVD virtual
channels have no precomputed sensitivity maps; RSS needs none and is identical
in form for both variants.
"""
import sys, time, contextlib, io
import numpy as np, h5py, finufft
from pathlib import Path
import scipy.integrate as _si
_si.cumtrapz = getattr(_si, 'cumtrapz', _si.cumulative_trapezoid)
sys.path.insert(0, '/home/debi/jaime/repos/MR-EyeTrack/old_study/code/twixtools')
import twixtools

B = Path('/home/debi/jaime/repos/MR-EyeTrack')
OUT = B / 'data/derived/recon_compare'
NS, NCH, N = 480, 52, 240
NSEG, NSHOT, NOFF = 44, 1872, 14
NR = (NSHOT - NOFF) * (NSEG - 1)
SUB = 15


def keep_index():
    k = np.arange(NSHOT * NSEG).reshape(NSHOT, NSEG)
    return k[NOFF:, 1:].ravel()


def grid_rss(y, t, ve, label):
    """y [nPt, nv] -> RSS image on the 240^3 grid."""
    kx, ky, kz = [np.ascontiguousarray(np.pi * t[:, i]) for i in range(3)]
    acc = np.zeros((N, N, N), np.float64)
    for c in range(y.shape[1]):
        t0 = time.time()
        img = finufft.nufft3d1(kx, ky, kz,
                               np.ascontiguousarray((y[:, c] * ve).astype(np.complex128)),
                               (N, N, N), isign=1, eps=1e-4)
        acc += np.abs(img) ** 2
        print(f'  {label} ch{c} {time.time()-t0:.0f}s', flush=True)
    return np.sqrt(acc).astype(np.float32)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = B / f'data/study/sub-{SUB:03d}/recon/ROI-PCA/kspace_ROI-PCA_8_woBin.mat'
    with h5py.File(p, 'r') as f:
        t = f['t'][:].reshape(-1, 3).astype(np.float64)
        ve = f['ve'][:].ravel().astype(np.float64)
        yp = f['y'][:]
        yp = (yp['real'] + 1j * yp['imag']).T.astype(np.complex64)     # [nPt, 8]
    print(f'ROI-PCA k-space {yp.shape}, traj {t.shape}')
    np.save(OUT / 'roipca8_rss.npy', grid_rss(yp, t, ve, 'ROI-PCA'))
    del yp

    # SVD: project the raw 52 channels with the basis computed in svd_compress
    sys.path.insert(0, str(B / 'model'))
    from svd_compress import open_raw, basis
    img = open_raw(SUB); keep = keep_index()
    V, w = basis(img, keep, every=20)
    print(f'SVD basis: nv=8 retains {100*np.cumsum(w)[7]/w.sum():.2f}% of total energy')
    ys = np.empty((NR * NS, 8), np.complex64)
    for j, i in enumerate(keep):
        ys[j * NS:(j + 1) * NS] = (V[:, :8].conj().T @ img[i].data).T
    del img
    np.save(OUT / 'svd8_rss.npy', grid_rss(ys, t, ve, 'SVD'))
    print('done ->', OUT)


if __name__ == '__main__':
    main()
