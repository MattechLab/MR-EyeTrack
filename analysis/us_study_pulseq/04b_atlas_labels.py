#!/usr/bin/env python
"""Multi-atlas A-eye labels for subjects without a usable MPRAGE (sub-010, sub-011).

data/study/sub-010 and sub-011 hold byte-identical copies of sub-001's MPRAGE
(same md5, same acquisition time), so the MPRAGE route of 04_register_masks.py
would paint sub-001's orbits onto them. Instead, the 13 subjects with valid
MPRAGE-derived LIBRE labels serve as atlases: each atlas woBin REGU is
registered to the target woBin REGU (same contrast) -- rigid + affine on the
whole head, then SyN (CC) inside a box around the orbits -- its labels are
carried over, and the 13 label maps are fused by per-voxel majority vote.

`--validate` runs the same procedure leave-one-out on subjects that do have
MPRAGE labels and reports Dice per structure against them, so the fused labels
come with a measured accuracy.

    python 04b_atlas_labels.py                      # sub-010, sub-011
    python 04b_atlas_labels.py --validate sub-004 sub-008 sub-013 sub-015
"""
import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
ANTS = "/usr/local/ants-2.6.3-ubuntu-24.04-X64-gcc/ants-2.6.3/bin"
NO_MPRAGE = ["sub-010", "sub-011"]
NAMES = {1: "lens", 2: "vitreous", 3: "nerve", 4: "int_fat", 5: "ext_fat",
         6: "lat_rect", 7: "med_rect", 8: "inf_rect", 9: "sup_rect"}


def wobin(sub):
    return BIDS / sub / "anat" / f"{sub}_rec-stevaFull_T1w.nii.gz"


def sh(cmd):
    env = dict(os.environ, ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="4")
    r = subprocess.run(["nice", "-n", "10", *cmd], capture_output=True, text=True, env=env)
    if r.returncode:
        raise RuntimeError(f"{cmd[0]}: {r.stderr[-1500:]}")


def orbit_box(labels_img, margin_mm=15):
    """Fixed-space mask: bounding box of the (affinely carried) labels + margin."""
    L = np.asanyarray(labels_img.dataobj) > 0
    idx = np.argwhere(L)
    lo = np.maximum(idx.min(0) - margin_mm, 0)
    hi = np.minimum(idx.max(0) + margin_mm + 1, L.shape)
    box = np.zeros(L.shape, np.uint8)
    box[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = 1
    return box


def register_one(args):
    atlas, target, work = args
    work.mkdir(parents=True, exist_ok=True)
    out = work / f"{atlas}_labels.nii.gz"
    if out.exists():
        return out
    fixed, moving = wobin(target), wobin(atlas)
    labels = ROOT / "masks" / atlas / "aeye_labels.nii.gz"
    pre = str(work / f"{atlas}_")
    sh([f"{ANTS}/antsRegistration", "-d", "3", "--float", "1", "-v", "0",
        "-o", f"[{pre}aff_]", "-r", f"[{fixed},{moving},1]",
        "-t", "Rigid[0.1]", "-m", f"Mattes[{fixed},{moving},1,32,Regular,0.25]",
        "-c", "[1000x500x250,1e-7,10]", "-f", "8x4x2", "-s", "3x2x1vox",
        "-t", "Affine[0.1]", "-m", f"Mattes[{fixed},{moving},1,32,Regular,0.25]",
        "-c", "[1000x500x250,1e-7,10]", "-f", "8x4x2", "-s", "3x2x1vox"])
    sh([f"{ANTS}/antsApplyTransforms", "-d", "3", "-i", str(labels), "-r", str(fixed),
        "-o", f"{pre}aff_labels.nii.gz", "-n", "GenericLabel", "-t", f"{pre}aff_0GenericAffine.mat"])
    ref = nib.load(fixed)
    box = orbit_box(nib.load(f"{pre}aff_labels.nii.gz"))
    nib.save(nib.Nifti1Image(box, ref.affine), f"{pre}box.nii.gz")
    # SyN on the orbit box only: the warp field lives on the crop, and is the
    # identity outside it, where the affine alone applies
    sh([f"{ANTS}/ExtractRegionFromImageByMask", "3", str(fixed), f"{pre}fixed_box.nii.gz", f"{pre}box.nii.gz", "1", "0"])
    sh([f"{ANTS}/antsRegistration", "-d", "3", "--float", "1", "-v", "0",
        "-o", f"[{pre}syn_]", "--initial-moving-transform", f"{pre}aff_0GenericAffine.mat",
        "-t", "SyN[0.1,3,0]", "-m", f"CC[{pre}fixed_box.nii.gz,{moving},1,3]",
        "-c", "[100x70x40,1e-6,10]", "-f", "4x2x1", "-s", "2x1x0vox"])
    sh([f"{ANTS}/antsApplyTransforms", "-d", "3", "-i", str(labels), "-r", str(fixed),
        "-o", str(out), "-n", "GenericLabel",
        "-t", f"{pre}syn_1Warp.nii.gz", "-t", f"{pre}aff_0GenericAffine.mat"])
    return out


def fuse(paths, ref):
    stack = np.stack([np.asanyarray(nib.load(p).dataobj).astype(np.uint8) for p in paths])
    votes = np.stack([(stack == k).sum(0) for k in range(10)])     # 0 = background
    fused = votes.argmax(0).astype(np.uint8)
    return nib.Nifti1Image(fused, nib.load(ref).affine)


def dice(a, b):
    s = a.sum() + b.sum()
    return 2 * (a & b).sum() / s if s else np.nan


def atlases_for(target):
    return [f"sub-{i:03d}" for i in range(1, 16)
            if f"sub-{i:03d}" not in NO_MPRAGE + [target]
            and (ROOT / "masks" / f"sub-{i:03d}" / "aeye_labels.nii.gz").exists()]


def run_target(target, workers, out_path):
    work = ROOT / "atlas_work" / target
    tasks = [(a, target, work) for a in atlases_for(target)]
    with ProcessPoolExecutor(workers) as ex:
        paths = list(ex.map(register_one, tasks))
    img = fuse(paths, wobin(target))
    nib.save(img, out_path)
    return img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.validate:
        rows = []
        for t in a.validate:
            fused = run_target(t, a.workers, ROOT / "atlas_work" / t / "fused_labels.nii.gz")
            F = np.asanyarray(fused.dataobj)
            T = np.asanyarray(nib.load(ROOT / "masks" / t / "aeye_labels.nii.gz").dataobj)
            for k, n in NAMES.items():
                rows.append(dict(subject=t, structure=n, dice=dice(F == k, T == k)))
            rows.append(dict(subject=t, structure="eyeball", dice=dice(np.isin(F, [1, 2]), np.isin(T, [1, 2]))))
            print(t, "validated", flush=True)
        df = pd.DataFrame(rows)
        df.to_csv(ROOT / "metrics/atlas_label_validation.csv", index=False)
        print(df.groupby("structure").dice.agg(["mean", "min", "max"]).round(3))
    else:
        for t in NO_MPRAGE:
            out = ROOT / "masks" / t / "aeye_labels.nii.gz"
            out.parent.mkdir(parents=True, exist_ok=True)
            run_target(t, a.workers, out)
            (out.parent / "aeye_labels_SOURCE.txt").write_text(
                "multi-atlas: 13 LIBRE woBin atlases, SyN in orbit box, majority vote (04b_atlas_labels.py)\n")
            print(t, "fused", flush=True)


if __name__ == "__main__":
    main()
