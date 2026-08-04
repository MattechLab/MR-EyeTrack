#!/usr/bin/env python
"""Side-by-side panels for settling left-right by eye.

The automated checks cannot resolve an L-R flip: the brain and the pair of orbits
are close to bilaterally symmetric, so overlap- and intensity-based measures
score the flipped volume as well as the true one (on Yiwei's T2w it scored
higher).  What they cannot use, a human can -- an individual's brain is not
perfectly symmetric, and a distinctive ventricle horn, vessel or sulcal pattern
is a subject-specific asymmetric fact.  Dice averages that away over the whole
mask; the eye finds it.

This renders the reference, the test volume, and a mirrored copy of the test
volume at the *same world coordinates*, so the comparison is like for like.  If
you can find one clearly asymmetric feature in the reference, the row it matches
tells you whether the test volume's left-right is correct:

  matches TEST            -> left-right is correct
  matches TEST MIRRORED   -> the test volume is left-right flipped
  cannot tell             -> inconclusive; do not record it as verified

Slices default to the lateral ventricles, which are individually asymmetric far
more often than the cortex is.

Two assumptions worth stating.  The reference's own left-right is taken as
ground truth -- it comes from DICOM through dcm2niix, whose handling of the
patient coordinate system is well tested.  And "no visible difference" is not
evidence of correctness: if the subject's anatomy really is symmetric at these
levels, the panel proves nothing either way.

Usage
-----
  python lr_inspect.py --ref MPRAGE.nii.gz --test recon.nii.gz --out DIR
                       [--levels -10,0,10,20] [--label NAME]

--levels are offsets in mm along S from the brain-mask centroid.
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np

from orientation_check import prepare, resample_to


def lr_axis(affine):
    """Array axis most aligned with the world left-right direction."""
    return int(np.argmax(np.abs(affine[0, :3])))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--test", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--label", default=None)
    p.add_argument("--rigid", default=None,
                   help="ANTs rigid0GenericAffine.mat from orientation_check.py; "
                        "aligns the planes so fine asymmetries are comparable")
    p.add_argument("--mask", default=None,
                   help="reference brain mask (SynthStrip ref_mask.nii.gz); strongly "
                        "recommended, it is what puts the slices at the right level")
    p.add_argument("--levels", default="-5,5,15,25",
                   help="S offsets in mm from the brain centroid (default -5,5,15,25, "
                        "which lands on the lateral ventricles)")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    work = Path(args.out).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    label = args.label or Path(args.test).name.replace(".nii.gz", "")
    levels = [float(x) for x in args.levels.split(",")]

    ref_img = prepare(args.ref, work / "lr_ref.nii.gz")
    test_img = prepare(args.test, work / "lr_test.nii.gz")

    ref_data = np.asanyarray(ref_img.dataobj).astype(np.float32)

    if args.rigid:
        # Match the anatomical planes before comparing.  The two scans differ by
        # a few degrees of head rotation, so slices at equal world coordinates are
        # not the same plane, which blurs exactly the fine asymmetries this panel
        # is for.  A rigid transform has determinant +1, so it cannot create or
        # undo a mirror -- using it here matches planes without touching the
        # question being asked.
        import subprocess
        from orientation_check import ants
        warped = work / "lr_test_rigid.nii.gz"
        subprocess.run([ants("antsApplyTransforms"), "-d", "3",
                        "-i", str(work / "lr_test.nii.gz"),
                        "-r", str(work / "lr_ref.nii.gz"),
                        "-t", str(args.rigid), "-n", "Linear",
                        "-o", str(warped)], check=True, capture_output=True)
        test_on_ref = np.asanyarray(nib.load(str(warped)).dataobj).astype(np.float32)
        print("  test rigidly aligned to the reference (planes now match)")
    else:
        test_on_ref = resample_to(test_img, ref_img, order=1)
        print("  no --rigid given; slices are at equal world coordinates but the "
              "residual head rotation means they are not the same anatomical plane")

    # canonical RAS so axis 0 = R, axis 1 = A, axis 2 = S and the panels can be
    # labelled with confidence rather than by guessing the storage order
    def canon(arr):
        return np.asanyarray(
            nib.as_closest_canonical(nib.Nifti1Image(arr, ref_img.affine)).dataobj)

    ref_c, test_c = canon(ref_data), canon(test_on_ref)
    # The mirrored row is the aligned test flipped along the canonical R axis, so
    # both test rows show identical planes and differ only by the mirror -- which
    # is precisely the question.
    flip_c = np.flip(test_c, 0)
    aff_c = nib.as_closest_canonical(
        nib.Nifti1Image(ref_data, ref_img.affine)).affine

    # Brain centroid in world coordinates.  An intensity threshold is a poor
    # substitute here: an MPRAGE FoV runs well down into the neck, and the face
    # and shoulders drag the centroid ~40mm inferior, putting the panels through
    # the skull base instead of the ventricles.  Use the SynthStrip mask when it
    # is available.
    if args.mask:
        m_img = nib.load(args.mask)
        m_c = np.asanyarray(nib.as_closest_canonical(m_img).dataobj) > 0
        nz = np.array(np.nonzero(m_c), dtype=float)
    else:
        print("  no --mask given; estimating the brain centre from intensity, "
              "which biases inferior on whole-head FoVs")
        nz = np.array(np.nonzero(ref_c > np.percentile(ref_c, 60)), dtype=float)
    centre_idx = nz.mean(axis=1)
    centre_world = aff_c[:3, :3] @ centre_idx + aff_c[:3, 3]
    print(f"  brain centre S = {centre_world[2]:+.1f} mm")

    rows = [("REFERENCE (ground truth)", ref_c),
            ("TEST (as written)", test_c),
            ("TEST MIRRORED", flip_c)]

    def norm(v):
        return np.clip(v / (np.percentile(v, 99.5) + 1e-9), 0, 1)

    fig, axes = plt.subplots(3, len(levels), figsize=(3.6 * len(levels), 11))
    if len(levels) == 1:
        axes = axes[:, None]

    for col, dz in enumerate(levels):
        world = centre_world.copy()
        world[2] += dz
        k = int(round(np.linalg.solve(aff_c[:3, :3], world - aff_c[:3, 3])[2]))
        k = max(0, min(ref_c.shape[2] - 1, k))
        for row, (title, vol) in enumerate(rows):
            a = axes[row, col]
            a.imshow(np.rot90(norm(vol)[:, :, k]), cmap="gray", vmin=0, vmax=1)
            a.set_xticks([]); a.set_yticks([])
            # canonical RAS + rot90 => anterior up, subject RIGHT on image right
            a.set_ylabel(title, fontsize=9) if col == 0 else None
            if row == 0:
                a.set_title(f"S = {world[2]:+.0f} mm", fontsize=10)
            a.text(0.02, 0.5, "L", color="#ffd23f", fontsize=13, fontweight="bold",
                   transform=a.transAxes, va="center")
            a.text(0.98, 0.5, "R", color="#ffd23f", fontsize=13, fontweight="bold",
                   transform=a.transAxes, va="center", ha="right")

    fig.suptitle(
        f"{label} - left/right inspection (anterior up, subject RIGHT on image right)\n"
        "Find an asymmetric feature in the REFERENCE row, then see which of the two "
        "rows below it matches.\nMatches TEST = correct.  Matches TEST MIRRORED = "
        "flipped.  No visible asymmetry = inconclusive.",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out_png = work / f"lr_inspect_{label}.png"
    fig.savefig(out_png, dpi=115, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure : {out_png}")


if __name__ == "__main__":
    main()
