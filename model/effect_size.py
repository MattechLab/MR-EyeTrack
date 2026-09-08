"""How big is the gaze effect, physically, and what does that cost per readout?

Every k-space probe I built was confounded (drift, class priors, or template
noise). This measures the same quantity in the image domain instead, where the
existing binned reconstructions already separate the gaze states.

    d'_total = ||m_up - m_down||_ROI / sigma_of_that_difference

is the discriminability actually delivered by the n_a and n_b readouts that
went into the two bins. Per readout it is d'_total / sqrt(n_eff), with
n_eff = n_a*n_b/(n_a+n_b), because evidence adds as sqrt(N).

The control that makes this interpretable: the same statistic in a static
region of the head. The two bins have different sampling patterns, so their
difference contains streaks and noise everywhere. Only the excess of the ROI
over the static control is gaze.
"""
import numpy as np
import h5py
import scipy.io as sio
from pathlib import Path

BASE = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
NSEG, NSHOT, NOFF = 44, 1872, 14


def vol(sub, region, kind='x0'):
    p = (BASE / f'sub-{sub:03d}/recon/clean/x0/x0_regionidx{region}.mat' if kind == 'x0'
         else BASE / f'sub-{sub:03d}/recon/clean/x/x_steva_regionidx_{region}_nIter_20_delta_1.000.mat')
    with h5py.File(p, 'r') as f:
        key = 'x0' if kind == 'x0' else list(k for k in f if not k.startswith('#'))[0]
        d = f[key]
        a = f[d[0, 0]][:] if d.shape == (1, 1) and d.dtype == object else d[:]
    return (a['real'] + 1j * a['imag']).astype(np.complex64) if a.dtype.names else a


def masks(sub):
    with h5py.File(BASE / f'sub-{sub:03d}/recon/ROI-PCA/masks.mat', 'r') as f:
        # h5py hands back MATLAB arrays with the axes reversed; without this
        # transpose the "ROI" lands in the middle of the head, not the orbits
        m = {k: f[k][:].T.astype(bool) for k in ('roiMask', 'intMask', 'headMask')}
    up = lambda a: np.repeat(np.repeat(np.repeat(a, 5, 0), 5, 1), 5, 2)   # 48 -> 240
    return {k: up(v) for k, v in m.items()}


def nreadouts(sub, region):
    p = BASE / f'sub-{sub:03d}/recon/bins/clean/eMask_th0.75_region{region}.mat'
    m = np.asarray(sio.loadmat(p)['eMaskN']).squeeze().reshape(NSHOT, NSEG)
    return int(m[NOFF:, 1:].sum())


def run(sub, kind='x0'):
    M = masks(sub)
    roi, head = M['roiMask'], M['headMask']
    bg = ~head
    # static control: head tissue away from the orbits, matched in voxel count
    ctrl = M['intMask'] & ~roi
    idx = np.flatnonzero(ctrl.ravel())
    keep = np.random.default_rng(0).choice(idx, min(roi.sum(), idx.size), replace=False)
    ctrl = np.zeros(ctrl.size, bool); ctrl[keep] = True; ctrl = ctrl.reshape(roi.shape)

    V = [vol(sub, r, kind) for r in range(4)]
    # intensity-match on static tissue, so a global scale difference is not "signal"
    ref = np.median(np.abs(V[0][M['intMask'] & ~roi]))
    V = [v * (ref / np.median(np.abs(v[M['intMask'] & ~roi]))) for v in V]

    sig_bg = np.sqrt(np.mean(np.abs(V[0][bg]) ** 2))
    print(f'sub-{sub:03d} ({kind}): ROI {int(roi.sum())} vox, '
          f'control {int(ctrl.sum())} vox, bg sigma {sig_bg:.4g}')
    print(f'  {"contrast":12s} {"||d||_ROI":>10s} {"||d||_ctrl":>11s} {"excess":>8s} '
          f'{"d_total":>9s} {"n_eff":>7s} {"d per RO":>9s} {"W for AUC .9":>13s}')
    for name, (a, b) in {'up/down': (0, 1), 'left/right': (2, 3),
                         'up/left': (0, 2), 'down/right': (1, 3)}.items():
        D = V[a] - V[b]
        e_roi = np.sqrt(np.mean(np.abs(D[roi]) ** 2))
        e_ctl = np.sqrt(np.mean(np.abs(D[ctrl]) ** 2))
        exc = e_roi / e_ctl
        # gaze-specific energy = ROI excess over the artefact/noise floor
        sig2 = max(e_roi ** 2 - e_ctl ** 2, 0.0)
        d_tot = np.sqrt(sig2 * roi.sum()) / (e_ctl + 1e-12)
        na, nb = nreadouts(sub, a), nreadouts(sub, b)
        neff = na * nb / (na + nb)
        d1 = d_tot / np.sqrt(neff)
        W = (2 * 1.2816 / d1) ** 2 if d1 > 0 else np.inf     # AUC 0.9 -> d' = 1.812
        print(f'  {name:12s} {e_roi:10.4g} {e_ctl:11.4g} {exc:8.2f} '
              f'{d_tot:9.1f} {neff:7.0f} {d1:9.4f} {W:13.0f}')


if __name__ == '__main__':
    import sys
    for s in ([int(a) for a in sys.argv[1:]] or [15]):
        run(s, 'x0'); print()
