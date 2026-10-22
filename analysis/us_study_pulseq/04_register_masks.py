#!/usr/bin/env python
"""Carry every region mask onto each subject's LIBRE grid (one grid per subject).

All 34 volumes of a subject share the woBin geometry (01_mat2bids.py), so each
mask is built once, in that grid, and reused for every gaze bin and
undersampling level. That is deliberate: a mask recomputed per image would
change with the image quality being measured (MRIQC does exactly that, which is
why its IQMs are reported separately from these fixed-mask metrics).

Sources
  brain   MRIQC brain_tissue_segmentation of woBin REGU (1 CSF, 2 GM, 3 WM),
          from a dedicated MRIQC pass on the 15 woBin volumes (work_wobin/) so
          the masks do not wait for the 510-volume runs
  eyes    MNI152NLin2009cAsym desc-eye mask through MRIQC's inverse MNI
          transform of woBin REGU -- last year's "Eyes" region
  aeye    A-eye labels on the MPRAGE (03_run_aeye.sh), MPRAGE -> woBin REGU
          by ANTs: whole-head rigid (Mattes MI), then a rigid refinement
          restricted to a 15 mm shell around the orbits. Registered ONCE against
          woBin and applied to all bins, per CLAUDE.md.

Output: results/us_pulseq/masks/sub-XXX/{brain_dseg,eyes_mni,aeye_labels}.nii.gz
"""
import argparse
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
ANTS = "/usr/local/ants-2.6.3-ubuntu-24.04-X64-gcc/ants-2.6.3/bin"
MNI_EYE = Path.home() / ".cache/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-01_desc-eye_mask.nii.gz"


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{cmd[0]} failed:\n{r.stderr[-2000:]}")
    return r.stdout


def mriqc_node(sub, node, pattern, tag="wobin", stem=None):
    """Find a file in an MRIQC work dir for one input image.

    MRIQC keys every per-image node directory by the same hash; the hash is
    recovered from <hash>/conform/<stem>_conformed.nii.gz."""
    stem = stem or f"{sub}_rec-stevaFull_T1w"
    work = ROOT / f"work_{tag}/mriqc_wf/anatMRIQC"
    conf = list(work.glob(f"*/conform/{stem}_conformed.nii.gz"))
    if len(conf) != 1:
        raise FileNotFoundError(f"{stem}: {len(conf)} conformed inputs in {work}")
    h = conf[0].parent.parent.name
    cand = sorted((work / node / h).rglob(pattern))
    if not cand:
        raise FileNotFoundError(f"{stem}: {node}/{h}/**/{pattern} not found")
    return cand[0]


def brain_and_eyes(sub, out):
    ref = BIDS / sub / "anat" / f"{sub}_rec-stevaFull_T1w.nii.gz"
    seg = mriqc_node(sub, "brain_tissue_segmentation", "segment.nii.gz")
    sh([f"{ANTS}/antsApplyTransforms", "-d", "3", "-i", str(seg), "-r", str(ref),
        "-o", str(out / "brain_dseg.nii.gz"), "-n", "GenericLabel"])
    inv = mriqc_node(sub, "SpatialNormalization", "*InverseComposite.h5")
    sh([f"{ANTS}/antsApplyTransforms", "-d", "3", "-i", str(MNI_EYE), "-r", str(ref),
        "-o", str(out / "eyes_mni.nii.gz"), "-n", "GenericLabel", "-t", str(inv)])


def aeye(sub, out):
    fixed = BIDS / sub / "anat" / f"{sub}_rec-stevaFull_T1w.nii.gz"
    moving = REPO / f"data/study/{sub}/dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/{sub}.nii.gz"
    labels = ROOT / f"aeye/{sub}/{sub}_aeye_both_gated.nii.gz"
    m = nib.load(moving)
    if m.ndim == 4:                                   # dcm2niix writes (x,y,z,1)
        m3 = out / "mprage.nii.gz"
        nib.save(nib.Nifti1Image(np.asanyarray(m.dataobj)[..., 0].astype(np.float32), m.affine), m3)
        moving = m3

    # orbit shell in MPRAGE space: A-eye labels dilated by 15 mm
    L = nib.load(labels)
    zoom = np.array(L.header.get_zooms()[:3])
    dist = ndi.distance_transform_edt(np.asanyarray(L.dataobj) == 0, sampling=zoom)
    nib.save(nib.Nifti1Image((dist <= 15).astype(np.uint8), L.affine), out / "orbit_shell_mprage.nii.gz")

    pre = str(out / "mprage2libre_")
    common = ["-d", "3", "--float", "1", "-u", "1", "-z", "1", "--winsorize-image-intensities", "[0.005,0.995]",
              "--interpolation", "Linear", "-v", "0"]
    sh([f"{ANTS}/antsRegistration", *common, "-o", f"[{pre}stage1_,{pre}stage1_Warped.nii.gz]",
        "-t", "Rigid[0.1]", "-m", f"Mattes[{fixed},{moving},1,32,Regular,0.25]",
        "-c", "[1000x500x250x0,1e-7,10]", "-f", "8x4x2x1", "-s", "3x2x1x0vox"])
    sh([f"{ANTS}/antsRegistration", *common, "-o", f"[{pre},{pre}Warped.nii.gz]",
        "--initial-moving-transform", f"{pre}stage1_0GenericAffine.mat",
        "-t", "Rigid[0.05]", "-m", f"Mattes[{fixed},{moving},1,32,Regular,0.5]",
        "-x", f"[NULL,{out / 'orbit_shell_mprage.nii.gz'}]",
        "-c", "[500x250x0,1e-7,10]", "-f", "4x2x1", "-s", "2x1x0vox"])
    sh([f"{ANTS}/antsApplyTransforms", "-d", "3", "-i", str(labels), "-r", str(fixed),
        "-o", str(out / "aeye_labels.nii.gz"), "-n", "GenericLabel", "-t", f"{pre}0GenericAffine.mat"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", nargs="*", default=[f"sub-{i:03d}" for i in range(1, 16)])
    ap.add_argument("--only", choices=["brain", "aeye", "all"], default="all")
    a = ap.parse_args()
    for sub in a.subjects:
        out = ROOT / "masks" / sub
        out.mkdir(parents=True, exist_ok=True)
        if a.only in ("brain", "all"):
            brain_and_eyes(sub, out)
        if a.only in ("aeye", "all"):
            aeye(sub, out)
        print(f"{sub} done", flush=True)


if __name__ == "__main__":
    main()
