#!/usr/bin/env python
"""Does A-eye know left from right?  A decisive test.

Run whole-head (no quadrant cropping) A-eye picks exactly one orbit.  If that
choice is driven by genuine anatomical chirality, it can settle the left-right
question the overlap metrics cannot.  If it is driven by position in the voxel
array -- which is what a CNN learns when every training image shares one storage
orientation -- it cannot, because the bug being hunted is "array mirrored,
affine unchanged", and any function of the array alone returns the same answer
for both.

The test feeds A-eye the same head twice: once as stored, once with the data
array mirrored along the world left-right axis and the affine left untouched.
Then it asks where the segmented orbit landed.

  same array slot both times      -> positional.  The affine reports the same
                                     world side either way, so the test has NO
                                     power to detect a flip.
  mirrored array slot             -> chirality-aware.  A-eye tracks the anatomy
                                     rather than the slot, and its answer IS a
                                     left-right check.

Usage
-----
  python aeye_chirality_test.py --source MPRAGE.nii.gz \
      --orig OUT/orig.nii.gz --flip OUT/flip.nii.gz
"""

import argparse

import nibabel as nib
import numpy as np


def lr_axis(affine):
    return int(np.argmax(np.abs(affine[0, :3])))


def summarise(seg_path, ax, n_lr):
    """Centroid of everything A-eye segmented, in array index and world x."""
    img = nib.load(seg_path)
    d = np.asanyarray(img.dataobj)
    m = d > 0
    if not m.any():
        return None
    idx = np.array(np.nonzero(m), dtype=float)
    mean_idx = idx.mean(axis=1)
    world = img.affine[:3, :3] @ mean_idx + img.affine[:3, 3]
    return {
        "voxels": int(m.sum()),
        "labels": sorted(int(v) for v in np.unique(d) if v != 0),
        "array_idx_lr": float(mean_idx[ax]),
        "world_x": float(world[0]),
        "mirrored_idx_lr": float(n_lr - 1 - mean_idx[ax]),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", required=True, help="the image A-eye was run on")
    p.add_argument("--orig", required=True, help="A-eye output for the as-stored image")
    p.add_argument("--flip", required=True, help="A-eye output for the mirrored image")
    args = p.parse_args()

    src = nib.load(args.source)
    ax = lr_axis(src.affine)
    n_lr = src.shape[ax]
    sign = np.sign(src.affine[0, ax])

    print("=== Does A-eye track anatomy or array position? ===")
    print(f"  source axcodes {nib.aff2axcodes(src.affine)}, world L-R = array axis {ax} "
          f"(size {n_lr}, +x at {'high' if sign > 0 else 'low'} index)\n")

    a = summarise(args.orig, ax, n_lr)
    b = summarise(args.flip, ax, n_lr)
    for tag, r in (("as stored", a), ("mirrored ", b)):
        if r is None:
            print(f"  {tag}: nothing segmented")
            continue
        side = "RIGHT" if r["world_x"] > 0 else "LEFT"
        print(f"  {tag}: {r['voxels']:6d} vox, labels {r['labels']}")
        print(f"             array index along L-R = {r['array_idx_lr']:6.1f}"
              f"   world x = {r['world_x']:+7.2f} mm  -> reported {side}")

    if a is None or b is None:
        print("\n  inconclusive: one run segmented nothing")
        return

    d_same = abs(b["array_idx_lr"] - a["array_idx_lr"])
    d_mirror = abs(b["array_idx_lr"] - a["mirrored_idx_lr"])
    print(f"\n  distance to the SAME array slot     : {d_same:6.1f} voxels")
    print(f"  distance to the MIRRORED array slot : {d_mirror:6.1f} voxels")

    same_world_side = (a["world_x"] > 0) == (b["world_x"] > 0)
    if d_same < d_mirror:
        print("\n  VERDICT: POSITIONAL. A-eye picked the same slot in the array both\n"
              "  times, so it followed the voxel grid rather than the anatomy.")
        print("  The affine then reports "
              f"{'the same world side' if same_world_side else 'different world sides'}"
              " for both, which means this\n  test CANNOT detect a left-right flip.")
    else:
        print("\n  VERDICT: CHIRALITY-AWARE. A-eye followed the anatomy to the other\n"
              "  side of the array, so which world side it reports genuinely depends\n"
              "  on whether the data is mirrored -- this IS a usable left-right check.")


if __name__ == "__main__":
    main()
