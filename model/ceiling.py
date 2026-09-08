"""The ceiling: per-readout detectability of the eye rotation we actually measured.

Synthesise the two gaze states instead of taking them from noisy binned recons
-- that is what made the earlier attempt (calibrate.py) return an impossible
d' of 125. Here:

  1. take the woBin reconstruction, the highest-SNR image available;
  2. rotate the eye region about the globe centre by +/- theta/2, with a radial
     taper that reproduces the measured nerve profile (full rotation inside the
     globe, decaying to zero by the orbital apex -- the nerve is tethered);
  3. forward-project both states through the real forward model (virtual coil
     maps, then NUFFT onto the actual spokes);
  4. d'_r = ||K1_r - K2_r|| / sigma, with sigma measured from the acquired data.

The difference cancels any per-readout phase term, which is why the phase
mismatch in the forward model (magnitude corr 0.85, complex ~0) does not matter.

Two angles are run: the ET-implied excursion (best case) and the displacement
actually measured from the nerve profile (~30% of it, the realistic case).
"""
import numpy as np, h5py, nibabel as nib, finufft
from pathlib import Path
from scipy.ndimage import map_coordinates, zoom

B = Path('/home/debi/jaime/repos/MR-EyeTrack')
LAB = B / 'analysis/test_orientation/tissue_check/results/aeye_sub015/labels_in_libre.nii.gz'
KSP = B / 'data/study/sub-015/recon/ROI-PCA/kspace_ROI-PCA_8_woBin.mat'
COIL = B / 'data/study/sub-015/recon/ROI-PCA/C_rovir_8.mat'
WOB = B / 'data/study/sub-015/recon/woBin/x0.mat'
N, NS = 240, 480
EVERY = 4                       # readout subsampling for the projection


def nii_to_mat(a):
    """NIfTI RAS index order -> the .mat array order the trajectory lives in."""
    return np.transpose(np.flip(np.flip(a, 0), 1), (1, 0, 2))


def load_wobin():
    with h5py.File(WOB, 'r') as f:
        d = f['x0']
        a = f[d[0, 0]][:] if d.shape == (1, 1) and d.dtype == object else d[:]
    return (a['real'] + 1j * a['imag']).astype(np.complex64)


def eye_weight(lab, centre, r_zero=30.0):
    """Where the rotation applies, and how strongly.

    Only the eye rotates. A radial taper over everything within 30 mm would also
    rotate orbital fat, the recti, bone and part of the brain -- which is what
    made the first run return d' = 18-57 per readout. So the weight comes from
    the segmentation: 1 inside globe and lens, tapering along the tethered nerve,
    0 everywhere else. Rotating the globe's content about its own centre leaves
    its spherical boundary fixed, so no discontinuity is introduced there.
    """
    from scipy.ndimage import gaussian_filter
    w = np.zeros(lab.shape, np.float32)
    w[(lab == 1) | (lab == 2)] = 1.0
    nerve = np.array(np.nonzero(lab == 3))
    if nerve.size:
        d = np.linalg.norm(nerve - np.asarray(centre, np.float32)[:, None], axis=0)
        w[tuple(nerve)] = np.clip((r_zero - d) / (r_zero - 12.0), 0, 1)
    return gaussian_filter(w, 0.8)


def rotate_eye(vol, centre, axis, deg, weight):
    """Rotate voxel content about `centre` by `deg`, scaled by `weight`."""
    g = np.stack(np.meshgrid(*[np.arange(N)] * 3, indexing='ij')).astype(np.float32)
    rel = g - np.array(centre, np.float32)[:, None, None, None]
    th = np.radians(deg) * weight
    u = np.asarray(axis, np.float32); u = u / np.linalg.norm(u)
    # small-angle inverse rotation: p_src = p - theta * (u x p)
    cross = np.stack([u[1] * rel[2] - u[2] * rel[1],
                      u[2] * rel[0] - u[0] * rel[2],
                      u[0] * rel[1] - u[1] * rel[0]])
    src = g - th * cross
    out = np.empty_like(vol)
    for part, f in (('real', np.real), ('imag', np.imag)):
        out = out + (map_coordinates(f(vol).astype(np.float32), src, order=1,
                                     mode='nearest') * (1 if part == 'real' else 1j))
    return out.astype(np.complex64)


def main():
    lab = nii_to_mat(np.asanyarray(nib.load(LAB).dataobj).astype(np.uint8))
    vol = load_wobin()
    with h5py.File(COIL, 'r') as f:
        a = f['C_rovir'][:]
    C = (a['real'] + 1j * a['imag']).transpose(3, 2, 1, 0)
    C = np.stack([zoom(C[..., c], N / C.shape[0], order=1) for c in range(C.shape[-1])], -1)

    with h5py.File(KSP, 'r') as f:
        t = f['t'][::EVERY].astype(np.float64)                 # [nro, NS, 3]
        y = f['y'][:, :]
    nro = t.shape[0]
    sel = np.concatenate([np.arange(r * NS, (r + 1) * NS) for r in range(0, y.shape[1] // NS, EVERY)])
    yy = y[:, sel]; yy = (yy['real'] + 1j * yy['imag'])
    kr = np.abs(np.arange(NS) - NS // 2) / (NS // 2)
    ymat = yy.reshape(yy.shape[0], nro, NS)
    sigma = np.sqrt(np.mean(np.abs(ymat[:, :, kr > 0.92]) ** 2))   # complex noise scale
    print(f'{nro} readouts, noise sigma per k-space sample {sigma:.4g}')

    # globes and their optical axes from the labels
    gl = np.array(np.nonzero(lab == 2)); mid = int(np.median(gl[0]))
    kx, ky, kz = [np.ascontiguousarray(np.pi * t[..., i].ravel()) for i in range(3)]

    hdr = "median d' per readout"
    print(f'\n{"eye":6s} {"contrast":12s} {"angle":>7s} {hdr:>22s} '
          f'{"W for AUC 0.9":>14s} {"= seconds":>10s}')
    for side, f_ in (('right', lambda a: a > mid), ('left', lambda a: a < mid)):
        keep = f_(gl[0]); gc = gl[:, keep].mean(1)
        # rotation applies to this eye's structures only
        eyelab = np.where(f_(np.arange(N))[:, None, None], lab, 0)
        wgt = eye_weight(eyelab, gc)
        print(f'  {side} eye: rotating {int((wgt > 0.05).sum())} voxels '
              f'(globe+lens+tapered nerve)')
        for cname, axis, deg_et in (('left vs right', (0, 0, 1), 18.00),
                                    ('up vs down', (0, 1, 0), 7.90)):
            for tag, deg in (('ET-implied', deg_et), ('measured ~30%', deg_et * 0.30)):
                K = []
                for s in (+deg / 2, -deg / 2):
                    v = rotate_eye(vol, gc, axis, s, wgt)
                    acc = np.empty((nro * NS, C.shape[-1]), np.complex128)
                    for c in range(C.shape[-1]):
                        acc[:, c] = finufft.nufft3d2(
                            kx, ky, kz,
                            np.ascontiguousarray((v * C[..., c]).astype(np.complex128)),
                            isign=-1, eps=1e-5)
                    K.append(acc)
                D = (K[0] - K[1]).reshape(nro, NS, -1)
                # put the simulation on the measured magnitude scale
                Km = ((K[0] + K[1]) / 2).reshape(nro, NS, -1)
                alpha = np.linalg.norm(ymat.transpose(1, 2, 0).ravel()) / np.linalg.norm(Km.ravel())
                d = np.sqrt((np.abs(D * alpha) ** 2).sum((1, 2))) / sigma
                md = float(np.median(d))
                W = (2 * 1.2816 / md) ** 2 if md > 0 else np.inf
                print(f'{side:6s} {cname:12s} {deg:6.2f}d {md:22.4f} {W:14.0f} {W*0.008:10.1f}'
                      f'   [{tag}]', flush=True)


if __name__ == '__main__':
    main()
