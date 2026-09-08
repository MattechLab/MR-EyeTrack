"""Generic SVD coil compression, as a control against the ROI-targeted variants.

ROI-PCA and ROVir both pick a channel subspace using an *eye ROI* criterion.
Plain SVD compression picks the subspace that retains the most signal energy
over the whole FoV, with no ROI involved. It is the standard method, and it is
the right control: if eye information survives SVD but not ROI-PCA, the problem
was the ROI targeting; if it survives neither, compression per se is not the
issue.

    C = sum_pts y y^H      [52 x 52]      -> eigenvectors, descending eigenvalue
    y_svd = V[:, :nv]^H y                  [nv x nPt]

Outputs go to data/derived/svd/ and touch nothing in recon/ROI-PCA or
recon/ROVir.
"""
import contextlib
import io
import sys
import time
from pathlib import Path

import numpy as np
import scipy.integrate as _si

_si.cumtrapz = getattr(_si, 'cumtrapz', _si.cumulative_trapezoid)
sys.path.insert(0, '/home/debi/jaime/repos/MR-EyeTrack/old_study/code/twixtools')
import twixtools  # noqa: E402

STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/svd')
NSEG, NSHOT, NOFF = 44, 1872, 14
NS, NCH = 480, 52
NR = (NSHOT - NOFF) * (NSEG - 1)
LO, HI = 160, 320                      # same |k| band as every other dataset in model/
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)


def keep_index():
    k = np.arange(NSHOT * NSEG).reshape(NSHOT, NSEG)
    return k[NOFF:, 1:].ravel()


def open_raw(sub):
    src = next((STUDY / f'sub-{sub:03d}/rawdata').glob('*_T1wLIBRE.dat'))
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        tw = twixtools.read_twix(str(src))
    img = [m for m in tw[-1]['mdb'] if m.is_image_scan()]
    assert len(img) == NSHOT * NSEG, f'{len(img)} image scans'
    return img


def basis(img, keep, every=20):
    """Channel covariance from a readout subsample -> eigenvectors + energy curve."""
    C = np.zeros((NCH, NCH), np.complex128)
    for i in keep[::every]:
        y = img[i].data                     # [NCH, NS]
        C += y @ y.conj().T
    C = (C + C.conj().T) / 2
    w, V = np.linalg.eigh(C)
    o = np.argsort(w)[::-1]
    w, V = w[o].real, V[:, o]
    return V, w


def run(sub, nvs=(8, 16), every=20):
    OUT.mkdir(parents=True, exist_ok=True)
    if all((OUT / f'sub-{sub:03d}_band{LO}-{HI}_SVD_{nv}.npy').exists() for nv in nvs):
        print(f'sub-{sub:03d} already done', flush=True); return
    t0 = time.time()
    img = open_raw(sub)
    keep = keep_index()
    V, w = basis(img, keep, every)
    cum = np.cumsum(w) / w.sum()
    np.savez(OUT / f'sub-{sub:03d}_svd_basis.npz', V=V, eigval=w, cum_energy=cum)
    print(f'sub-{sub:03d} basis in {time.time()-t0:.0f}s   energy retained: ' +
          '  '.join(f'nv={nv}: {100*cum[nv-1]:.2f}%' for nv in nvs), flush=True)

    bands = {nv: np.empty((NR, HI - LO, nv), np.complex64) for nv in nvs}
    for j, i in enumerate(keep):
        y = img[i].data[:, LO:HI]           # [NCH, band]
        for nv in nvs:
            bands[nv][j] = (V[:, :nv].conj().T @ y).T
        if j % 20000 == 0:
            print(f'  sub-{sub:03d} {j}/{NR}  {time.time()-t0:.0f}s', flush=True)
    for nv in nvs:
        np.save(OUT / f'sub-{sub:03d}_band{LO}-{HI}_SVD_{nv}.npy', bands[nv])
    print(f'sub-{sub:03d} done  {time.time()-t0:.0f}s  '
          f'({", ".join(f"nv={nv}: {bands[nv].nbytes/2**30:.2f} GiB" for nv in nvs)})',
          flush=True)


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or list(KEEP)):
        run(s)
