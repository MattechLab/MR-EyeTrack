"""Lens and optic-nerve displacement across gaze bins, inside A-eye label ROIs.

Supersedes globe_shift.py, which correlated 15^3 boxes around a clicked point --
mostly static sclera and fat, which drags any real shift toward zero. Here the
correlation is weighted by the segmentation, so only the moving structure counts.

Everything runs in the NIfTI RAS grid the labels were transformed into:
axis0 = R->L, axis1 = P->A, axis2 = I->S. So left/right gaze displaces along
axis0 and up/down along axis2.

The nerve is split along its length. It is tethered at the orbital apex and
moves with the globe at the head, so a true rotation produces a **gradient** --
largest at the ONH, ~zero at the apex. A nuisance (head motion, drift) displaces
it uniformly. That gradient is a second, independent signature of real rotation.
"""
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy import ndimage as ndi

B = Path('/home/debi/jaime/repos/MR-EyeTrack')
LAB = B / 'analysis/test_orientation/tissue_check/results/aeye_sub015/labels_in_libre.nii.gz'
BINS = [B / f'data/study/sub-015/recon/clean/x/x_steva_regionidx_{r}_nIter_20_delta_1.000.nii.gz'
        for r in range(4)]
NAMES = ['up', 'down', 'left', 'right']
EDGE = True   # weight the correlation by image gradient


def weighted_shift(A, B_, w, rad=5):
    """Sub-voxel displacement of B_ vs A, weighted by w. Returns (d0,d1,d2) voxels."""
    wa = w * (A - (A * w).sum() / w.sum())
    n = 2 * rad + 1
    cc = np.zeros((n, n, n))
    for i, d0 in enumerate(range(-rad, rad + 1)):
        for j, d1 in enumerate(range(-rad, rad + 1)):
            for k, d2 in enumerate(range(-rad, rad + 1)):
                Bs = np.roll(np.roll(np.roll(B_, d0, 0), d1, 1), d2, 2)
                bs = w * (Bs - (Bs * w).sum() / w.sum())
                cc[i, j, k] = (wa * bs).sum() / np.sqrt((wa**2).sum() * (bs**2).sum() + 1e-20)
    p = np.unravel_index(np.argmax(cc), cc.shape)
    out = []
    for ax in range(3):
        i = p[ax]
        if 0 < i < n - 1:
            s = list(p); s[ax] = i - 1; ym = cc[tuple(s)]
            s[ax] = i + 1; yp = cc[tuple(s)]
            den = ym - 2 * cc[p] + yp
            out.append(i - rad + (0.5 * (ym - yp) / den if abs(den) > 1e-12 else 0.0))
        else:
            out.append(float(i - rad))
    return np.array(out), float(cc[p])


def main():
    lab = nib.load(LAB)
    L = np.asanyarray(lab.dataobj).astype(np.uint8)
    V = [np.asanyarray(nib.load(p).dataobj).astype(np.float32) for p in BINS]
    gl = np.array(np.nonzero(L == 2))
    mid = int(np.median(gl[0]))

    # ET-measured excursion -> predicted displacement at radius r
    PRED = {'L-R': np.radians(18.00), 'S-I': np.radians(7.90)}

    print('sub-015 — displacement between gaze bins inside A-eye ROIs')
    print('(NIfTI RAS grid: axis0 R->L, axis2 I->S; 1 voxel = 1 mm)\n')
    rois = {}
    for side, sel in (('right', lambda a: a > mid), ('left', lambda a: a < mid)):
        gsel = sel(gl[0])
        gc = gl[:, gsel].mean(1)
        for name, val in (('lens', 1), ('nerve', 3)):
            m = (L == val)
            idx = np.array(np.nonzero(m))
            keep = sel(idx[0])
            mm = np.zeros_like(m)
            mm[tuple(idx[:, keep])] = True
            if name == 'nerve':
                # split along the nerve's long axis (distance from globe centre)
                d = np.linalg.norm(np.array(np.nonzero(mm)) - gc[:, None], axis=0)
                med = np.median(d)
                for part, cond in (('nerve-head', d <= med), ('nerve-apex', d > med)):
                    sub = np.zeros_like(m)
                    sub[tuple(np.array(np.nonzero(mm))[:, cond])] = True
                    rois[(side, part)] = (sub, gc)
            else:
                rois[(side, name)] = (mm, gc)
    # Anterior cap of the globe: the corneal bulge.
    # A sphere rotating about its own centre is invisible -- the sclera maps onto
    # itself. Only features that break spherical symmetry displace: the cornea,
    # the lens, the ONH. So take the globe's outer shell and keep the cap within
    # CAP_DEG of the optical axis (globe centre -> lens centre).
    CAP_DEG = 45.0
    for side, sel in (('right', lambda a: a > mid), ('left', lambda a: a < mid)):
        gidx = np.array(np.nonzero(L == 2)); gk = sel(gidx[0])
        gc = gidx[:, gk].mean(1)
        lidx = np.array(np.nonzero(L == 1)); lk = sel(lidx[0])
        lc = lidx[:, lk].mean(1)
        u = (lc - gc) / np.linalg.norm(lc - gc)
        g = np.zeros_like(L, bool); g[tuple(gidx[:, gk])] = True
        shell = ndi.binary_dilation(g, iterations=2) & ~ndi.binary_erosion(g, iterations=1)
        sidx = np.array(np.nonzero(shell))
        vec = sidx - gc[:, None]
        cosang = (u[:, None] * vec).sum(0) / (np.linalg.norm(vec, axis=0) + 1e-9)
        cap = np.zeros_like(L, bool)
        cap[tuple(sidx[:, cosang > np.cos(np.radians(CAP_DEG))])] = True
        rois[(side, 'cornea cap')] = (cap, gc)

    # a static control far from the orbit
    ctrl = np.zeros_like(L, bool); ctrl[100:140, 60:100, 100:140] = True
    rois[('--', 'brain control')] = (ctrl, None)

    print(f'  {"eye":6s} {"structure":14s} {"vox":>6s} {"contrast":12s} {"axis":5s} '
          f'{"shift":>7s} {"r":>6s} {"predicted":>10s}')
    res = {}
    for (side, name), (mask, gc) in rois.items():
        w = ndi.binary_dilation(mask, iterations=2).astype(np.float32)
        w = ndi.gaussian_filter(w, 1.0)
        if EDGE:
            # displacement is carried by edges: rotating the globe leaves the
            # homogeneous vitreous unchanged, so weighting by |grad| removes the
            # static interior that otherwise drags the estimate toward zero
            M = np.mean(V, 0)
            g = np.sqrt(sum(ndi.sobel(M, axis=a).astype(np.float32) ** 2 for a in range(3)))
            w = w * (g / (g[w > 0.05].mean() + 1e-9))
        bb = np.array(np.nonzero(w > 0.05))
        lo = bb.min(1) - 4; hi = bb.max(1) + 5
        sl = tuple(slice(max(l, 0), h) for l, h in zip(lo, hi))
        ww = w[sl]
        for cname, (a, b), ax in (('left vs right', (2, 3), 0), ('up vs down', (0, 1), 2)):
            d, r = weighted_shift(V[a][sl], V[b][sl], ww)
            pred = ''
            if gc is not None:
                rad = np.linalg.norm(np.array(np.nonzero(mask)).mean(1) - gc)
                key = 'L-R' if ax == 0 else 'S-I'
                pred = f'{PRED[key] * rad:.2f}'
            res[(side, name, cname)] = d[ax]
            print(f'  {side:6s} {name:14s} {int(mask.sum()):6d} {cname:12s} '
                  f'{"R-L" if ax == 0 else "S-I":5s} {d[ax]:+7.2f} {r:6.3f} {pred:>10s}')
    print('\n  cornea-cap-minus-nerve-head differential (the rotation signature):')
    for side in ('right', 'left'):
        for cname in ('left vs right', 'up vs down'):
            dl = res.get((side, 'cornea cap', cname)); dn = res.get((side, 'nerve-head', cname))
            if dl is not None and dn is not None:
                print(f'    {side:6s} {cname:14s} {dl - dn:+6.2f} mm')
    print('\n  nerve gradient (head minus apex) — a rotation bends the nerve, '
          'a nuisance translates it:')
    for side in ('right', 'left'):
        for cname in ('left vs right', 'up vs down'):
            dh = res.get((side, 'nerve-head', cname)); da = res.get((side, 'nerve-apex', cname))
            if dh is not None and da is not None:
                print(f'    {side:6s} {cname:14s} {dh - da:+6.2f} mm')


if __name__ == '__main__':
    main()
