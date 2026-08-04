#!/usr/bin/env python
"""Decide left-right by image similarity: as written vs mirrored.

Mask overlap cannot answer this.  A brain mask is the outer envelope of the
brain, which is close to mirror-symmetric, so a left-right flip scores within a
fraction of a percent of the truth -- on one dataset here it scored higher.  The
same is true of the orbits, and of anything else built on the shape of a
near-symmetric object.

The images themselves are a different matter.  Ventricles and sulcal patterns
are strongly asymmetric and highly individual, and an intensity metric sees all
of it.  So: rigidly register the test volume to the reference, rigidly register
a mirrored copy of it as well, and compare how well each ends up matching.

A rigid transform has determinant +1 and cannot create or undo a mirror, so
neither registration can rescue the wrong hypothesis -- each simply gets its
best possible alignment, and the one that fits better is the true orientation.

Both registrations are reported so you can confirm neither of them failed: they
should converge to near-identical rotations and translations.  If one did not
converge, its poor similarity says nothing about anatomy and the result must be
discarded.

Usage
-----
  python lr_flip_test.py --ref MPRAGE.nii.gz --test recon.nii.gz --out DIR
                         [--mask ref_mask.nii.gz] [--label NAME]

--mask restricts the metrics to the brain, keeping skull, neck and background
from diluting them.  Without it, SynthStrip is run on the reference.
"""

import argparse
import csv
import json
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from orientation_check import ants, prepare, read_itk_affine, decompose, register, synthstrip


def lr_axis(affine):
    return int(np.argmax(np.abs(affine[0, :3])))


def similarity(ref, moving, mask, metric, bins_or_radius):
    """ANTs image similarity; ANTs returns these negated, so lower is better."""
    out = subprocess.run(
        [ants("MeasureImageSimilarity"), "-d", "3",
         "-m", f"{metric}[{ref},{moving},1,{bins_or_radius}]", "-x", str(mask)],
        capture_output=True, text=True)
    try:
        return float(out.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--test", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mask", default=None, help="reference brain mask")
    p.add_argument("--label", default=None)
    args = p.parse_args()

    work = Path(args.out).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    label = args.label or Path(args.test).name.replace(".nii.gz", "")

    ref_img = prepare(args.ref, work / "ref.nii.gz")
    test_img = prepare(args.test, work / "test.nii.gz")

    ax = lr_axis(test_img.affine)
    d = np.asanyarray(test_img.dataobj)
    nib.save(nib.Nifti1Image(np.ascontiguousarray(np.flip(d, ax)), test_img.affine),
             work / "test_mirrored.nii.gz")

    mask = Path(args.mask) if args.mask else synthstrip(work, "ref")

    print(f"=== left-right test: {label} ===")
    print(f"  world L-R is array axis {ax} of the test volume\n")

    results = {}
    for tag, moving in (("as written", "test"), ("mirrored", "test_mirrored")):
        print(f"  registering {tag}...")
        mat = register(work, work / "ref.nii.gz", work / f"{moving}.nii.gz",
                       "Rigid", f"lr_{moving}")
        A, t, c = read_itk_affine(mat)
        rig = decompose(A, t, c, nib.load(str(mask)))
        warped = work / f"lr_{moving}_warped.nii.gz"
        results[tag] = {
            "mi": similarity(work / "ref.nii.gz", warped, mask, "MI", 32),
            "cc": similarity(work / "ref.nii.gz", warped, mask, "CC", 4),
            "rotation_deg": rig["rotation_deg"],
            "translation_mm": rig["translation_mm"],
        }

    a, b = results["as written"], results["mirrored"]
    print(f"\n  {'':12}{'MI':>12}{'CC':>12}{'rotation':>11}{'translation':>13}")
    for tag, r in results.items():
        print(f"  {tag:12}{r['mi']:12.4f}{r['cc']:12.4f}"
              f"{r['rotation_deg']:10.2f}d{r['translation_mm']:12.2f}mm")
    print("  (ANTs negates both metrics, so lower is a better match)")

    # A large gap in convergence means one registration failed, and a failed
    # registration explains poor similarity without saying anything about anatomy.
    rot_gap = abs(a["rotation_deg"] - b["rotation_deg"])
    converged = rot_gap < 5.0
    print(f"\n  both registrations converged comparably: {converged} "
          f"(rotations differ by {rot_gap:.2f} deg)")

    mi_margin = b["mi"] - a["mi"]     # positive => as-written fits better
    cc_margin = b["cc"] - a["cc"]
    print(f"  MI margin favouring as-written: {mi_margin:+.4f}")
    print(f"  CC margin favouring as-written: {cc_margin:+.4f}")

    if not converged:
        verdict = "INCONCLUSIVE - one registration did not converge comparably"
    elif mi_margin > 0.02 and cc_margin > 0.02:
        verdict = "PASS - left-right is correct as written"
    elif mi_margin < -0.02 and cc_margin < -0.02:
        verdict = "FAIL - the mirrored volume fits better; left-right is flipped"
    else:
        verdict = ("INCONCLUSIVE - margins too small to call; this subject's "
                   "anatomy may be unusually symmetric at this resolution")
    print(f"\n  VERDICT: {verdict}\n")

    results["margins"] = {"mi": mi_margin, "cc": cc_margin,
                          "converged": converged, "verdict": verdict}
    with open(work / f"report_lr_{label}.json", "w") as fh:
        json.dump({"label": label, "ref": args.ref, "test": args.test,
                   "results": results}, fh, indent=2)
    with open(work / f"report_lr_{label}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "mi_as_written", "mi_mirrored", "cc_as_written",
                    "cc_mirrored", "mi_margin", "cc_margin", "rot_as_written",
                    "rot_mirrored", "verdict"])
        w.writerow([label, f"{a['mi']:.4f}", f"{b['mi']:.4f}", f"{a['cc']:.4f}",
                    f"{b['cc']:.4f}", f"{mi_margin:.4f}", f"{cc_margin:.4f}",
                    f"{a['rotation_deg']:.2f}", f"{b['rotation_deg']:.2f}", verdict])
    print(f"  report : {work / f'report_lr_{label}.json'}")
    print(f"  csv    : {work / f'report_lr_{label}.csv'}\n")


if __name__ == "__main__":
    main()
