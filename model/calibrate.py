"""Forward-simulation calibration: what is the per-readout d' between gaze states?

Every probe so far measured a confounded quantity. This measures the physics
directly and has nothing to confound:

  1. take the four binned reconstructions as the best available estimates of the
     four gaze-state images,
  2. push each through the *real* forward model -- multiply by each virtual coil
     map, then NUFFT onto the actual spoke trajectory -- giving noiseless
     k-space per gaze state,
  3. measure the noise per k-space sample from the acquired data itself,
  4. d'_r = ||k_a(r) - k_b(r)|| / sigma over the 480 samples x 8 channels of one
     readout.

k-space is the right domain for this: the noise there is white and per-sample,
so every sample is an independent measurement. In the image domain the noise of
an undersampled radial recon is spatially correlated and the effective degrees
of freedom are unknown, which is what made the earlier image-domain attempt
uninterpretable.

Self-validating: step 2 must reproduce the *measured* k-space. If the predicted
and acquired k-space do not correlate, the forward model (orientation, scaling,
coil maps) is wrong and no number below means anything.
"""
import argparse
import itertools
import numpy as np
import h5py
import scipy.io as sio
from pathlib import Path
from scipy.ndimage import zoom
import finufft

STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
NR, NS, N = 79894, 480, 240


def load_recon(sub, r, kind='steva'):
    if kind == 'steva':
        p = STUDY / f'sub-{sub:03d}/recon/clean/x/x_steva_regionidx_{r}_nIter_20_delta_1.000.mat'
        d = sio.loadmat(p)
        k = [x for x in d if not x.startswith('__')][0]
        return np.asarray(d[k]).squeeze().astype(np.complex128)
    p = STUDY / f'sub-{sub:03d}/recon/clean/x0/x0_regionidx{r}.mat'
    with h5py.File(p, 'r') as f:
        a = f[f['x0'][0, 0]][:]
    return (a['real'] + 1j * a['imag']).astype(np.complex128)


def load_coils(sub, nv=8):
    with h5py.File(STUDY / f'sub-{sub:03d}/recon/ROI-PCA/C_rovir_{nv}.mat', 'r') as f:
        a = f['C_rovir'][:]
    C = (a['real'] + 1j * a['imag']).transpose(3, 2, 1, 0)     # -> [48,48,48,nv]
    return np.stack([zoom(C[..., c], N / 48, order=1) for c in range(C.shape[-1])], -1)


def load_kspace(sub, ro):
    """Measured y for the selected readouts. -> [nro, NS, nv]"""
    p = STUDY / f'sub-{sub:03d}/recon/ROI-PCA/kspace_ROI-PCA_8_woBin.mat'
    with h5py.File(p, 'r') as f:
        y = f['y']
        nv = y.shape[0]
        out = np.empty((len(ro), NS, nv), np.complex64)
        for i, r in enumerate(ro):
            b = y[:, r * NS:(r + 1) * NS]
            out[i] = (b['real'] + 1j * b['imag']).T
        t = f['t'][ro].astype(np.float64)                      # [nro, NS, 3]
    return out, t


def project(vol, kpts):
    """Type-2 NUFFT: image on the 240^3 grid evaluated at the spoke samples."""
    return finufft.nufft3d2(kpts[0], kpts[1], kpts[2],
                            np.ascontiguousarray(vol), isign=-1, eps=1e-5)


def orient(v, axes, flips):
    v = np.transpose(v, axes)
    for f in flips:
        v = np.flip(v, f)
    return np.ascontiguousarray(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sub', type=int, default=15)
    ap.add_argument('--every', type=int, default=20, help='readout subsampling')
    ap.add_argument('--kind', default='steva', choices=['steva', 'x0'])
    a = ap.parse_args()

    ro = np.arange(0, NR, a.every)
    print(f'sub-{a.sub:03d}  {len(ro)} readouts ({a.every}x subsample)  recon={a.kind}')
    y, t = load_kspace(a.sub, ro)
    kx, ky, kz = [np.ascontiguousarray(np.pi * t[..., i].ravel()) for i in range(3)]
    V = [load_recon(a.sub, r, a.kind) for r in range(4)]
    C = load_coils(a.sub)
    mean_img = np.mean(V, 0)

    # ---- step 1: find the orientation that reproduces the measured k-space ----
    probe = slice(0, 400 * NS)
    kp = (kx[probe], ky[probe], kz[probe])
    meas = y[:400, :, 0].ravel()
    best = None
    for axes in itertools.permutations((0, 1, 2)):
        for flips in [(), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]:
            pred = project(orient(mean_img * C[..., 0], axes, flips), kp)
            mid = np.abs(pred) > 0
            c = abs(np.corrcoef(np.abs(pred[mid]), np.abs(meas[mid]))[0, 1])
            if best is None or c > best[0]:
                best = (c, axes, flips)
    corr, axes, flips = best
    print(f'  best forward-model orientation: transpose {axes} flip {flips}  '
          f'magnitude corr {corr:.4f}')
    if corr < 0.7:
        print('  *** forward model does not reproduce the data -- numbers below are void')

    # ---- step 2: noise per k-space sample, from the acquired data ------------
    kr = np.abs(np.arange(NS) - NS // 2) / (NS // 2)
    outer = kr > 0.92                       # outermost |k|: object signal is minimal
    sigma = np.sqrt(np.mean(np.abs(y[:, outer, :]) ** 2) / 2)   # per real component
    print(f'  noise sigma per k-space sample (per component): {sigma:.4g}')

    # ---- step 3: project each gaze state through the full forward model ------
    K = np.empty((4, len(ro), NS, C.shape[-1]), np.complex128)
    for g in range(4):
        for c in range(C.shape[-1]):
            K[g, :, :, c] = project(orient(V[g] * C[..., c], axes, flips),
                                    (kx, ky, kz)).reshape(len(ro), NS)
        print(f'    projected gaze {g}', flush=True)

    # single complex scale putting the simulation on the measured scale
    alpha = np.linalg.norm(y.ravel()) / np.linalg.norm(K.mean(0).ravel())
    K *= alpha
    print(f'  magnitude scale alpha {alpha:.4g}')

    # ---- step 4: per-readout d' between gaze states --------------------------
    hdr = "median d' per readout"
    print(f'\n  {"contrast":12s} {hdr:>22s} {"W for AUC 0.9":>15s} {"= seconds":>10s}')
    pairs = {'up/down': (0, 1), 'left/right': (2, 3), 'up/left': (0, 2),
             'up/right': (0, 3), 'down/left': (1, 2), 'down/right': (1, 3)}
    got = {}
    for nm, (g1, g2) in pairs.items():
        D = K[g1] - K[g2]                                  # [nro, NS, nv]
        d = np.sqrt((np.abs(D) ** 2).sum((1, 2))) / (sigma * np.sqrt(2))
        md = float(np.median(d))
        W = (2 * 1.2816 / md) ** 2 if md > 0 else np.inf
        got[nm] = md
        print(f'  {nm:12s} {md:22.4f} {W:15.0f} {W*0.008:10.1f}')
    v = np.array(list(got.values()))
    print(f'\n  spread across the six pairs: {v.min():.3f} to {v.max():.3f} '
          f'(ratio {v.max()/max(v.min(),1e-9):.2f})')
    print('  geometry predicts left/right (2.9 mm of lens travel) clearly exceeds')
    print('  up/down (1.6 mm). If all six pairs are equal, this is recon noise,')
    print('  not gaze -- the same test that flagged the image-domain attempt.')


if __name__ == '__main__':
    main()
