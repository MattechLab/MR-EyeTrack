"""52-channel k-space band, uncompressed, for the FiLM/TCN model (step 4).

Same band and same readout selection as every other dataset in model/, but with
no channel compression at all: [79894, 160, 52] complex64, 5.1 GiB per subject.

This is the last representation left to try. ROI-PCA (eye-ROI subspace), SVD
(whole-FoV subspace) and 16 channels instead of 8 all decode identically at
~0.287, so the channel map is not the lever; what changes here is that the model
sees all 52 physical coils, including the phase diversity that any compression
mixes away.
"""
import sys
import time
import numpy as np
from pathlib import Path
from svd_compress import open_raw, keep_index, NS, NCH, NR, KEEP

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband_raw52')
LO, HI = 160, 320


def run(sub):
    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / f'sub-{sub:03d}_band{LO}-{HI}_raw52.npy'
    if dst.exists():
        print(f'sub-{sub:03d} already done', flush=True); return
    t0 = time.time()
    img = open_raw(sub)
    keep = keep_index()
    band = np.empty((NR, HI - LO, NCH), np.complex64)
    for j, i in enumerate(keep):
        band[j] = img[i].data[:, LO:HI].T
    np.save(dst, band)
    print(f'sub-{sub:03d}  {band.shape}  {band.nbytes/2**30:.2f} GiB  '
          f'{time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    for s in ([int(a) for a in sys.argv[1:]] or list(KEEP)):
        run(s)
