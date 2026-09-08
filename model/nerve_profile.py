"""Optic nerve displacement along its length, per gaze bin.

Two fixes over nerve_shift.py, both aimed at the nerve specifically:

  * split the nerve by CONNECTED COMPONENT, not by the midline. The nerves
    converge medially toward the chiasm and cross the midline, which is why a
    midline split gave 106 voxels on one side and 315 on the other and made the
    two sides incomparable. Components are 454 / 419 -- properly symmetric.
  * measure displacement in BANDS along the nerve rather than as one number.
    The nerve is tethered at the apex and moves with the globe at its head, so a
    real rotation produces a monotone gradient: large at the ONH, ~zero at the
    apex. Bulk motion or drift shifts every band equally. The profile shape is
    therefore a signature that no single number can give.
"""
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy import ndimage as ndi
from nerve_shift import weighted_shift

B = Path('/home/debi/jaime/repos/MR-EyeTrack')
LAB = B / 'analysis/test_orientation/tissue_check/results/aeye_sub015/labels_in_libre.nii.gz'
BINS = [B / f'data/study/sub-015/recon/clean/x/x_steva_regionidx_{r}_nIter_20_delta_1.000.nii.gz'
        for r in range(4)]
NBAND = 4


def main():
    L = np.asanyarray(nib.load(LAB).dataobj).astype(np.uint8)
    V = [np.asanyarray(nib.load(p).dataobj).astype(np.float32) for p in BINS]
    M = np.mean(V, 0)
    grad = np.sqrt(sum(ndi.sobel(M, axis=a).astype(np.float32) ** 2 for a in range(3)))

    # globes and nerves, both split by connected component
    gl, ng = ndi.label(L == 2)
    gsz = ndi.sum(L == 2, gl, range(1, ng + 1))
    gkeep = np.argsort(gsz)[::-1][:2] + 1
    nv, nn = ndi.label(L == 3)
    nsz = ndi.sum(L == 3, nv, range(1, nn + 1))
    nkeep = np.argsort(nsz)[::-1][:2] + 1

    print('sub-015 — optic nerve displacement along its length')
    print('(NIfTI RAS: axis0 R->L, axis2 I->S; 1 voxel = 1 mm)\n')
    for ni in nkeep:
        nerve = nv == ni
        ncen = np.array(ndi.center_of_mass(nerve))
        # pair this nerve with its nearer globe
        gcs = [np.array(ndi.center_of_mass(gl == g)) for g in gkeep]
        gi = int(np.argmin([np.linalg.norm(ncen - c) for c in gcs]))
        gc = gcs[gi]
        side = 'right' if gc[0] > 120 else 'left'
        idx = np.array(np.nonzero(nerve))
        d = np.linalg.norm(idx - gc[:, None], axis=0)
        edges = np.quantile(d, np.linspace(0, 1, NBAND + 1))
        print(f'{side} nerve ({int(nerve.sum())} vox), distance from globe centre '
              f'{d.min():.0f}-{d.max():.0f} mm')
        print(f'  {"band":22s} {"vox":>5s} {"L-R shift":>10s} {"r":>6s} '
              f'{"U-D shift":>10s} {"r":>6s}')
        for b in range(NBAND):
            m = (d >= edges[b]) & (d <= edges[b + 1])
            sub = np.zeros_like(nerve)
            sub[tuple(idx[:, m])] = True
            w = ndi.gaussian_filter(ndi.binary_dilation(sub, iterations=2).astype(np.float32), 1.0)
            w = w * (grad / (grad[w > .05].mean() + 1e-9))
            bb = np.array(np.nonzero(w > .05))
            lo = np.maximum(bb.min(1) - 4, 0); hi = bb.max(1) + 5
            sl = tuple(slice(l, h) for l, h in zip(lo, hi))
            out = []
            for (a, bb_) , ax in (((2, 3), 0), ((0, 1), 2)):
                sh, r = weighted_shift(V[a][sl], V[bb_][sl], w[sl])
                out += [sh[ax], r]
            print(f'  {edges[b]:5.0f}-{edges[b+1]:5.0f} mm from globe {int(m.sum()):5d} '
                  f'{out[0]:+10.2f} {out[1]:6.3f} {out[2]:+10.2f} {out[3]:6.3f}')
        print()
    print('A rotation gives a monotone gradient, largest near the globe.')
    print('Uniform shift across bands = bulk motion, not eye rotation.')


if __name__ == '__main__':
    main()
