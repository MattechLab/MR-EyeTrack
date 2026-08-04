#!/usr/bin/env python
"""Run A-eye on both orbits of a reference MPRAGE, via quadrant cropping.

Task313_Eye is a single-orbit model: run on a whole head it segments one eye and
ignores the other.  Both eyes are obtained the way a-eye_web does it -- crop the
left and right anterior quadrants, segment each separately, then uncrop and
merge.  See a-eye_web/package/quadrant_segmentation/quadrant.py.

Difference from that implementation: the input is first reoriented to canonical
RAS.  quadrant.py slices `data[:mid_x, mid_y:, :]` assuming array axis 0 is
left-right and axis 1 is anterior-posterior, which holds for MR-EyeTrack MPRAGEs
(axcodes RAS) but NOT for Yiwei's (axcodes PIL) -- there it would crop an
inferior/posterior quadrant and find no eye at all.  Reorienting first makes the
assumption true for any input.

Output labels follow a-eye_web/package/biomarkers/biomarkers.py:
  1 lens   2 globe   3 optic nerve   4 intraconal fat   5 extraconal fat
  6 lateral rectus   7 medial rectus   8 inferior rectus   9 superior rectus

The merged mask is written in canonical RAS space, not the input grid.  That is
harmless for eye_roi_check.py, which maps masks through their own affine, and
canonical reorientation is a lossless permutation/flip.

Usage
-----
  python run_aeye.py --input MPRAGE.nii.gz --out WORKDIR [--sudo-gpu] [--fold 0]

--sudo-gpu runs the container with `sudo ... --gpus all`.  Plain `--gpus all`
fails on this host (no NVIDIA container toolkit for the unprivileged daemon),
and sudo prompts for a password, so this flag only works when run interactively.
Without it inference is CPU-only: ~1 min per fold per orbit.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

AEYE_IMAGE = os.environ.get("AEYE_IMAGE", "jaimebarran/fw_gear_aeye:0.0.1")

LABELS = {
    1: "lens", 2: "globe", 3: "nerve", 4: "int_fat", 5: "ext_fat",
    6: "lat_mus", 7: "med_mus", 8: "inf_mus", 9: "sup_mus",
}


def crop_quadrant(data, affine, left_side):
    """Crop the anterior-superior quadrant holding one eye, in canonical RAS.

    Unlike quadrant.py the affine is corrected for the crop offset, so the
    intermediate file is geometrically honest.  The merge does not depend on
    this, but a wrong origin on an intermediate is a trap for anyone who opens
    one to debug.
    """
    mid_x, mid_y = data.shape[0] // 2, data.shape[1] // 2
    if left_side:
        sub, offset = data[:mid_x, mid_y:, :], (0, mid_y, 0)
    else:
        sub, offset = data[mid_x:, mid_y:, :], (mid_x, mid_y, 0)
    new_affine = affine.copy()
    new_affine[:3, 3] = affine[:3, :3] @ np.array(offset, dtype=float) + affine[:3, 3]
    return np.ascontiguousarray(sub), new_affine


def uncrop_quadrant(seg, full_shape, left_side):
    mid_x, mid_y = full_shape[0] // 2, full_shape[1] // 2
    full = np.zeros(full_shape, dtype=np.uint8)
    if left_side:
        full[:mid_x, mid_y:, :] = seg
    else:
        full[mid_x:, mid_y:, :] = seg
    return full


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="reference MPRAGE NIfTI")
    p.add_argument("--out", required=True, help="work/output dir (Docker-shareable)")
    p.add_argument("--fold", default=None,
                   help="single fold for speed, e.g. 0; default is the 5-fold ensemble")
    p.add_argument("--sudo-gpu", action="store_true",
                   help="run the container under sudo with --gpus all (prompts)")
    p.add_argument("--tta", action="store_true", help="enable test-time augmentation (8x slower)")
    args = p.parse_args()

    work = Path(args.out).expanduser().resolve()
    (work / "in").mkdir(parents=True, exist_ok=True)
    (work / "out").mkdir(parents=True, exist_ok=True)
    case = Path(args.input).name.split(".")[0]

    img = nib.as_closest_canonical(nib.load(args.input))
    data = np.asanyarray(img.dataobj).astype(np.float32)
    print(f"=== A-eye both orbits: {case} ===")
    print(f"  input reoriented to canonical RAS: {img.shape} {nib.aff2axcodes(img.affine)}")

    for side, left in (("L", True), ("R", False)):
        sub, aff = crop_quadrant(data, img.affine, left)
        nib.save(nib.Nifti1Image(sub, aff), work / "in" / f"{case}{side}_0000.nii.gz")
        print(f"  {side} quadrant {sub.shape}")

    cmd = ["docker", "run", "--rm", "-e", "HOME=/tmp", "-v", f"{work}:/data"]
    if args.sudo_gpu:
        cmd = ["sudo"] + cmd[:2] + ["--gpus", "all"] + cmd[2:]
    else:
        cmd += ["--user", f"{os.getuid()}:{os.getgid()}"]
    cmd += ["--entrypoint", "nnUNet_predict", AEYE_IMAGE,
            "-i", "/data/in", "-o", "/data/out",
            "-tr", "nnUNetTrainerV2", "-m", "3d_fullres",
            "-p", "nnUNetPlansv2.1", "-t", "Task313_Eye"]
    if args.fold is not None:
        cmd += ["-f", args.fold]
    if not args.tta:
        cmd += ["--disable_tta"]

    print(f"  running nnUNet ({'GPU via sudo' if args.sudo_gpu else 'CPU'})...")
    if subprocess.run(cmd).returncode != 0:
        sys.exit("nnUNet_predict failed")

    merged = np.zeros(img.shape[:3], dtype=np.uint8)
    for side, left in (("L", True), ("R", False)):
        seg_path = work / "out" / f"{case}{side}.nii.gz"
        if not seg_path.exists():
            print(f"  WARNING: no output for {side} quadrant")
            continue
        seg = np.asanyarray(nib.load(str(seg_path)).dataobj).astype(np.uint8)
        full = uncrop_quadrant(seg, img.shape[:3], left)
        merged[full > 0] = full[full > 0]

    out_path = work / f"{case}_aeye_both.nii.gz"
    nib.save(nib.Nifti1Image(merged, img.affine), out_path)

    print("\n  merged labels (world centroids):")
    for lab, name in LABELS.items():
        m = merged == lab
        if not m.any():
            print(f"    {lab} {name:8s} absent")
            continue
        for sign, tag in ((1, "R"), (-1, "L")):
            sel = np.zeros_like(m)
            idx = np.array(np.nonzero(m))
            wx = img.affine[0, :3] @ idx + img.affine[0, 3]
            keep = (wx * sign) > 0
            if keep.sum() == 0:
                continue
            sel[tuple(idx[:, keep])] = True
            ci = np.array(np.nonzero(sel), dtype=float).mean(1)
            w = img.affine[:3, :3] @ ci + img.affine[:3, 3]
            print(f"    {lab} {name:8s} {tag}: {int(sel.sum()):6d} vox  "
                  f"RAS [{w[0]:+7.2f} {w[1]:+7.2f} {w[2]:+7.2f}]")

    print(f"\n  mask : {out_path}\n")


if __name__ == "__main__":
    main()
