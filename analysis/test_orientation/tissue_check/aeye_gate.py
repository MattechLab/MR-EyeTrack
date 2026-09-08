"""Drop A-eye label components that are too far from the globes to be orbital.

`crop_quadrant` in run_aeye.py halves the volume in x and y but keeps the FULL
z range, so each quadrant is a column running from the skull vertex down through
the jaw. nnUNet then has the mouth inside its search volume and emits fat labels
there. On sub-015 that is 534 voxels (2.1%), all fat, 73-95 mm from the nearest
globe -- while every correct component sits within 21 mm.

Gating on distance from the detected globe is preferable to adding a z crop:
it needs no assumption about field of view or head position, and the 21 mm vs
73 mm margin means the threshold needs no tuning.

    python aeye_gate.py --labels sub-015_aeye_both.nii.gz [--max-mm 50]
"""
import argparse
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi


def gate(L, zooms, max_mm=50.0, min_vox=5):
    """Keep only components whose centroid lies within max_mm of a globe centre."""
    lab, n = ndi.label(L == 2)                       # globe is the anchor label
    if n == 0:
        raise SystemExit('no globe label found - cannot gate')
    sizes = ndi.sum(L == 2, lab, range(1, n + 1))
    anchors = np.array([ndi.center_of_mass(L == 2, lab, int(i) + 1)
                        for i in np.argsort(sizes)[::-1][:2]])
    out = np.zeros_like(L)
    dropped = {}
    for v in np.unique(L):
        if v == 0:
            continue
        m = L == v
        cl, k = ndi.label(m)
        for c in range(1, k + 1):
            comp = cl == c
            sz = int(comp.sum())
            p = np.array(ndi.center_of_mass(m, cl, c))
            d = np.min(np.linalg.norm((anchors - p) * zooms, axis=1))
            if d <= max_mm and sz >= min_vox:
                out[comp] = v
            else:
                dropped[v] = dropped.get(v, 0) + sz
    return out, dropped, anchors


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', required=True)
    ap.add_argument('--out', default=None)
    ap.add_argument('--max-mm', type=float, default=50.0)
    a = ap.parse_args()
    im = nib.load(a.labels)
    L = np.asanyarray(im.dataobj).astype(np.uint8)
    zooms = np.array(im.header.get_zooms()[:3], float)
    clean, dropped, anchors = gate(L, zooms, a.max_mm)
    out = a.out or a.labels.replace('.nii.gz', '_gated.nii.gz')
    nib.save(nib.Nifti1Image(clean, im.affine, im.header), out)
    print(f'globe anchors (voxel): {np.round(anchors,1).tolist()}')
    print(f'kept {int((clean>0).sum())} voxels, dropped {int((L>0).sum()-(clean>0).sum())}')
    for v, n in sorted(dropped.items()):
        print(f'  label {v}: dropped {n} voxels')
    print(f'-> {out}')
