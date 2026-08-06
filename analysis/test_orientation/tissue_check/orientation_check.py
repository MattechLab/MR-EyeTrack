#!/usr/bin/env python
"""Validate the orientation of a mat2nii output against its reference MPRAGE.

The test asks one question: with the affines as written in the two headers and
no registration whatsoever, does the brain in the test image land on the brain
in the reference image?  Brain masks come from SynthStrip, which is contrast
agnostic and so works on fat-suppressed LIBRE T1w/T2w as well as on MPRAGE.

Three numbers are produced, and keeping them apart is what stops the test from
being circular:

  header-only Dice   THE TEST.  Test mask resampled into the reference grid
                     through the two affines alone.  Any orientation error
                     destroys it.
  registered Dice    THE CEILING.  Same overlap after rigid registration.  This
                     is the best Dice achievable given segmentation noise and
                     genuine head motion between the two acquisitions, so it is
                     what makes the header-only score interpretable.
  recovered rigid    THE DIAGNOSTIC.  The transform registration had to invent,
  transform          starting from the header alignment.  Near-identity means
                     the header is right; a near-90/180 deg rotation or signed
                     axis permutation means it is not.

Registration is deliberately NOT initialised by centre-of-mass or by moments:
it starts from the identity in physical space, i.e. from what the headers say,
so the recovered transform measures the header error directly.

Blind spot: a brain mask is the outer envelope of the brain, which is close to
mirror-symmetric, so a pure L-R flip scores within a fraction of a percent of
the truth here and is NOT resolved by this script.  Use lr_flip_test.py for
that -- it compares the images rather than their masks, where the ventricles and
sulcal patterns are strongly asymmetric.  This is reported in the output as a
reminder.

Usage
-----
  python orientation_check.py --ref MPRAGE.nii.gz --test recon.nii.gz \
      --out WORKDIR [--label NAME] [--n4] [--controls]

WORKDIR must be under a path Docker is allowed to bind-mount (on debi that
means somewhere under /home/debi -- /tmp is not shared).

--controls runs the negative controls: all 48 signed axis permutations of the
test volume are scored with the header-only metric, which calibrates what a
passing Dice actually means and shows how far the true orientation sits above
the best wrong one.
"""

import argparse
import csv
import itertools
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import io as sio
from scipy import ndimage

ANTS_BIN = os.environ.get(
    "ANTSPATH", "/usr/local/ants-2.6.3-ubuntu-24.04-X64-gcc/ants-2.6.3/bin"
)
SYNTHSTRIP_IMAGE = os.environ.get("SYNTHSTRIP_IMAGE", "freesurfer/synthstrip:latest")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def ants(tool):
    """Absolute path to an ANTs binary, preferring ANTS_BIN over $PATH."""
    cand = Path(ANTS_BIN) / tool
    if cand.exists():
        return str(cand)
    found = shutil.which(tool)
    if not found:
        sys.exit(f"ANTs tool not found: {tool} (looked in {ANTS_BIN} and $PATH)")
    return found


def run(cmd, **kw):
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{res.stdout}\n{res.stderr}")
    return res.stdout


def prepare(path, out_path):
    """Load a volume, reduce it to a real 3D float32 frame, and write it out.

    Reconstructions may arrive complex-valued or 4D; the affine is preserved
    untouched because it is the object under test.
    """
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    if np.iscomplexobj(data):
        data = np.abs(data)
    while data.ndim > 3:
        data = data[..., 0]
    data = np.nan_to_num(data.astype(np.float32))
    nib.save(nib.Nifti1Image(data, img.affine, dtype=np.float32), str(out_path))
    return nib.load(str(out_path))


def synthstrip(work, name):
    """Run SynthStrip on work/<name>.nii.gz, returning the brain mask path.

    Docker mounts the work directory, so it must be a host path Docker is
    configured to share.  --user keeps the outputs owned by the caller.
    """
    mask = work / f"{name}_mask.nii.gz"
    if mask.exists():
        return mask
    cmd = [
        "docker", "run", "--rm",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{work}:/data",
        SYNTHSTRIP_IMAGE,
        "-i", f"/data/{name}.nii.gz",
        "-m", f"/data/{name}_mask.nii.gz",
    ]
    # Kept as cheap insurance: the Docker VM used to have ~7.7 GiB, where a 480^3
    # volume sat close enough to the ceiling that a concurrent job could get this
    # one killed mid-frame. It is now 64 GiB, so this should no longer trigger.
    for attempt in range(3):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and mask.exists():
            return mask
        if attempt < 2:
            print(f"    SynthStrip attempt {attempt + 1} failed, retrying in 30 s "
                  f"(likely Docker memory contention)...")
            time.sleep(30)
    sys.exit(f"SynthStrip failed for {name}:\n{res.stdout}\n{res.stderr}")


def resample_to(src_img, ref_img, order=0):
    """Pull src into ref's voxel grid using the two affines only.

    For every reference voxel the world point is mapped back through the source
    affine, so this is a pure header operation -- no registration, no
    optimisation.  order=0 for masks, order=1 for images.
    """
    ref_shape = ref_img.shape[:3]
    grid = np.indices(ref_shape, dtype=np.float32).reshape(3, -1)
    world = ref_img.affine[:3, :3] @ grid + ref_img.affine[:3, 3:4]
    inv = np.linalg.inv(src_img.affine)
    src_idx = inv[:3, :3] @ world + inv[:3, 3:4]
    data = np.asanyarray(src_img.dataobj).astype(np.float32)
    out = ndimage.map_coordinates(data, src_idx, order=order, mode="constant", cval=0.0)
    return out.reshape(ref_shape)


def dice(a, b):
    a, b = a > 0, b > 0
    denom = a.sum() + b.sum()
    return float(2.0 * (a & b).sum() / denom) if denom else 0.0


def centroid_world(mask_img):
    d = np.asanyarray(mask_img.dataobj) > 0
    if not d.any():
        return np.full(3, np.nan)
    idx = np.array(np.nonzero(d), dtype=np.float64).mean(axis=1)
    return mask_img.affine[:3, :3] @ idx + mask_img.affine[:3, 3]


def surface_distances(a, b, spacing):
    """Symmetric surface distances (mm) between two binary masks on one grid."""
    a, b = a > 0, b > 0
    if not a.any() or not b.any():
        return np.array([np.nan])
    surf = lambda m: m ^ ndimage.binary_erosion(m, ndimage.generate_binary_structure(3, 1))
    sa, sb = surf(a), surf(b)
    da = ndimage.distance_transform_edt(~sa, sampling=spacing)
    db = ndimage.distance_transform_edt(~sb, sampling=spacing)
    return np.concatenate([db[sa], da[sb]])


def n4(in_path, out_path):
    run([ants("N4BiasFieldCorrection"), "-d", "3", "-i", str(in_path),
         "-o", str(out_path), "-s", "4", "-c", "[50x50x50x50,1e-7]"])
    return out_path


# ----------------------------------------------------------------------------
# registration
# ----------------------------------------------------------------------------

def register(work, fixed, moving, transform="Rigid", prefix="reg"):
    """Register moving to fixed with Mattes MI, starting from the identity.

    No --initial-moving-transform: the optimisation begins exactly where the
    headers put the two images, so whatever it recovers is the header error
    (plus real inter-scan motion).  MI is used because the two images are
    different contrasts by construction.
    """
    out = work / prefix
    run([
        ants("antsRegistration"),
        "--dimensionality", "3", "--float", "1", "--verbose", "0",
        "--output", f"[{out},{out}_warped.nii.gz]",
        "--interpolation", "Linear",
        "--winsorize-image-intensities", "[0.005,0.995]",
        "--use-histogram-matching", "0",
        "--transform", f"{transform}[0.1]",
        "--metric", f"MI[{fixed},{moving},1,32,Regular,0.25]",
        "--convergence", "[1000x500x250x100,1e-6,10]",
        "--shrink-factors", "8x4x2x1",
        "--smoothing-sigmas", "3x2x1x0vox",
    ])
    return Path(f"{out}0GenericAffine.mat")


def read_itk_affine(mat_path):
    """Read an ANTs GenericAffine.mat and return (matrix, offset, centre) in RAS.

    ITK works in LPS, so the transform is conjugated by diag(-1,-1,1) to land in
    the RAS convention NIfTI affines use.
    """
    m = sio.loadmat(str(mat_path))
    # key is AffineTransform_float_3_3 or _double_ depending on --float
    key = next(k for k in m if k.startswith("AffineTransform"))
    params = np.asarray(m[key]).ravel()
    centre = np.asarray(m["fixed"]).ravel()
    A_lps = params[:9].reshape(3, 3)
    t_lps = params[9:12]
    D = np.diag([-1.0, -1.0, 1.0])
    return D @ A_lps @ D, D @ t_lps, D @ centre


def decompose(A, t, centre, ref_mask_img):
    """Turn a recovered transform into interpretable numbers.

    Rotation angle comes from the trace; scales come from the singular values
    (only meaningful for the 12-DOF affine run).  The displacement statistics
    are computed over the actual brain voxels rather than at the origin, which
    is what makes them physically meaningful -- a small rotation about a distant
    centre still moves the brain a long way.
    """
    U, S, Vt = np.linalg.svd(A)
    R = U @ Vt
    if np.linalg.det(R) < 0:                     # keep it a proper rotation
        U[:, -1] *= -1
        R = U @ Vt
    angle = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))

    d = np.asanyarray(ref_mask_img.dataobj) > 0
    idx = np.array(np.nonzero(d), dtype=np.float64)
    if idx.shape[1] > 200000:                    # subsample, purely for speed
        idx = idx[:, :: idx.shape[1] // 200000 + 1]
    pts = ref_mask_img.affine[:3, :3] @ idx + ref_mask_img.affine[:3, 3:4]
    moved = A @ (pts - centre[:, None]) + centre[:, None] + t[:, None]
    disp = np.linalg.norm(moved - pts, axis=0)

    return {
        "rotation_deg": float(angle),
        "translation_mm": float(np.linalg.norm(t)),
        "translation_ras_mm": [float(x) for x in t],
        "scales": [float(x) for x in S],
        "brain_displacement_mean_mm": float(disp.mean()),
        "brain_displacement_p95_mm": float(np.percentile(disp, 95)),
        "brain_displacement_max_mm": float(disp.max()),
    }


# ----------------------------------------------------------------------------
# negative controls
# ----------------------------------------------------------------------------

def signed_permutations():
    """The 48 signed axis permutations -- every way a reorientFcn can be wrong."""
    for perm in itertools.permutations(range(3)):
        for flips in itertools.product([False, True], repeat=3):
            yield perm, flips


def run_controls(test_mask_img, ref_img, ref_mask):
    """Score every signed permutation of the test mask with the header-only metric.

    The permutation is applied to the mask rather than re-running SynthStrip on
    48 permuted volumes.  That assumes SynthStrip is permutation-equivariant,
    which is not exactly true of a CNN but is more than close enough to
    calibrate the metric, and it turns 8 minutes of GPU-less inference into a
    couple of seconds.

    The affine is held fixed while the data moves, which is precisely the bug
    being simulated: mat2nii_twix writes direction cosines from the reference
    regardless of how the data array is ordered, so a wrong reorientFcn shows up
    as exactly this -- permuted data under an unchanged affine.
    """
    data = np.asanyarray(test_mask_img.dataobj)
    results = []
    for perm, flips in signed_permutations():
        v = np.transpose(data, perm)
        for ax, f in enumerate(flips):
            if f:
                v = np.flip(v, ax)
        img = nib.Nifti1Image(np.ascontiguousarray(v), test_mask_img.affine)
        d = dice(resample_to(img, ref_img, order=0), ref_mask)
        results.append({
            "perm": list(perm),
            "flips": [bool(f) for f in flips],
            "identity": perm == (0, 1, 2) and not any(flips),
            "dice": d,
        })
    results.sort(key=lambda r: -r["dice"])
    return results


# ----------------------------------------------------------------------------
# QC figure
# ----------------------------------------------------------------------------

def qc_figure(ref_img, ref_mask, test_on_ref, test_mask_on_ref, out_png, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref = np.asanyarray(ref_img.dataobj).astype(np.float32)
    idx = np.array(np.nonzero(ref_mask), dtype=np.float64).mean(axis=1).astype(int)

    def norm(v):
        hi = np.percentile(v, 99.5)
        return np.clip(v / (hi + 1e-9), 0, 1)

    ref_n, test_n = norm(ref), norm(test_on_ref)
    views = [("Sagittal", 0), ("Coronal", 1), ("Axial", 2)]

    fig, axes = plt.subplots(2, 3, figsize=(13, 9))
    for col, (name, ax_i) in enumerate(views):
        sl = [slice(None)] * 3
        sl[ax_i] = idx[ax_i]
        sl = tuple(sl)
        for row, (bg, title) in enumerate([(ref_n, "REF"), (test_n, "TEST")]):
            a = axes[row, col]
            a.imshow(np.rot90(bg[sl]), cmap="gray", vmin=0, vmax=1)
            a.contour(np.rot90(ref_mask[sl].astype(float)), levels=[0.5],
                      colors="#39d353", linewidths=1.1)
            a.contour(np.rot90(test_mask_on_ref[sl].astype(float)), levels=[0.5],
                      colors="#ff5c5c", linewidths=1.1)
            a.set_title(f"{title} - {name}", fontsize=10)
            a.axis("off")
    fig.suptitle(f"{label}\ngreen = reference brain mask, red = test brain mask "
                 f"(header-only resampling, no registration)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="reference NIfTI (MPRAGE)")
    p.add_argument("--test", required=True, help="NIfTI produced by mat2nii")
    p.add_argument("--out", required=True, help="work/output dir (must be Docker-shareable)")
    p.add_argument("--label", default=None, help="name for the report and figure")
    p.add_argument("--n4", action="store_true",
                   help="N4 bias correct the test volume before registration")
    p.add_argument("--controls", action="store_true",
                   help="score all 48 signed axis permutations to calibrate the metric")
    p.add_argument("--affine-check", action="store_true",
                   help="also fit a 12-DOF affine to check voxel-size/FoV scaling")
    args = p.parse_args()

    work = Path(args.out).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    label = args.label or Path(args.test).name.replace(".nii.gz", "")

    if work.is_relative_to(Path("/tmp")):
        print("WARNING: Docker on this host does not share /tmp; "
              "SynthStrip will fail. Use a directory under your home.\n")

    print(f"=== orientation check: {label} ===")
    print(f"  ref  : {args.ref}")
    print(f"  test : {args.test}")
    print(f"  work : {work}\n")

    ref_img = prepare(args.ref, work / "ref.nii.gz")
    test_img = prepare(args.test, work / "test.nii.gz")

    print(f"  ref  shape {ref_img.shape} zooms {np.round(ref_img.header.get_zooms()[:3], 3)} "
          f"axcodes {nib.aff2axcodes(ref_img.affine)}")
    print(f"  test shape {test_img.shape} zooms {np.round(test_img.header.get_zooms()[:3], 3)} "
          f"axcodes {nib.aff2axcodes(test_img.affine)}\n")

    if args.n4:
        print("  N4 bias correcting test volume...")
        n4(work / "test.nii.gz", work / "test_n4.nii.gz")
        moving_for_reg = work / "test_n4.nii.gz"
    else:
        moving_for_reg = work / "test.nii.gz"

    print("  running SynthStrip on reference...")
    ref_mask_img = nib.load(str(synthstrip(work, "ref")))
    print("  running SynthStrip on test...")
    test_mask_img = nib.load(str(synthstrip(work, "test")))

    ref_mask = np.asanyarray(ref_mask_img.dataobj) > 0
    spacing = np.abs(ref_img.header.get_zooms()[:3])
    # each image gets its own voxel volume: the two grids need not match, and
    # for the 0.5 mm reconstructions they do not
    ref_vox_ml = np.abs(np.linalg.det(ref_img.affine[:3, :3])) / 1000.0
    test_vox_ml = np.abs(np.linalg.det(test_img.affine[:3, :3])) / 1000.0

    # ---- THE TEST: header-only overlap -------------------------------------
    test_mask_on_ref = resample_to(test_mask_img, ref_img, order=0) > 0
    header_dice = dice(test_mask_on_ref, ref_mask)

    c_ref = centroid_world(ref_mask_img)
    c_test = centroid_world(test_mask_img)
    centroid_mm = float(np.linalg.norm(c_ref - c_test))

    sd = surface_distances(ref_mask, test_mask_on_ref, spacing)
    hd95 = float(np.nanpercentile(sd, 95))

    # ---- THE CEILING: overlap after rigid registration ---------------------
    print("  rigid registration (MI, header-initialised)...")
    rigid_mat = register(work, work / "ref.nii.gz", moving_for_reg, "Rigid", "rigid")
    warped_mask = work / "test_mask_rigid.nii.gz"
    run([ants("antsApplyTransforms"), "-d", "3",
         "-i", str(work / "test_mask.nii.gz"), "-r", str(work / "ref.nii.gz"),
         "-t", str(rigid_mat), "-n", "NearestNeighbor", "-o", str(warped_mask)])
    reg_dice = dice(np.asanyarray(nib.load(str(warped_mask)).dataobj) > 0, ref_mask)

    # ---- THE DIAGNOSTIC: what registration had to invent -------------------
    A, t, centre = read_itk_affine(rigid_mat)
    rigid = decompose(A, t, centre, ref_mask_img)

    affine_info = None
    if args.affine_check:
        print("  12-DOF affine (scale diagnostic)...")
        aff_mat = register(work, work / "ref.nii.gz", moving_for_reg, "Affine", "affine")
        Aa, ta, ca = read_itk_affine(aff_mat)
        affine_info = decompose(Aa, ta, ca, ref_mask_img)

    # ---- report -------------------------------------------------------------
    gap = reg_dice - header_dice
    print("\n" + "=" * 68)
    print("  BRAIN MASK OVERLAP")
    print(f"    header-only Dice        {header_dice:6.4f}   <- THE TEST")
    print(f"    registered Dice         {reg_dice:6.4f}   <- the ceiling")
    print(f"    gap                     {gap:6.4f}   (header error + head motion)")
    print(f"    centroid distance       {centroid_mm:6.2f} mm  (world coords, no resampling)")
    print(f"    surface HD95            {hd95:6.2f} mm")
    print(f"    mask volumes            ref {ref_mask.sum() * ref_vox_ml:.1f} mL / "
          f"test {(np.asanyarray(test_mask_img.dataobj) > 0).sum() * test_vox_ml:.1f} mL")
    print("\n  RECOVERED RIGID TRANSFORM (identity expected if header is correct)")
    print(f"    rotation                {rigid['rotation_deg']:6.2f} deg")
    print(f"    translation             {rigid['translation_mm']:6.2f} mm  "
          f"RAS [{rigid['translation_ras_mm'][0]:+.2f} "
          f"{rigid['translation_ras_mm'][1]:+.2f} "
          f"{rigid['translation_ras_mm'][2]:+.2f}]")
    print(f"    brain displacement      mean {rigid['brain_displacement_mean_mm']:.2f} mm, "
          f"p95 {rigid['brain_displacement_p95_mm']:.2f} mm, "
          f"max {rigid['brain_displacement_max_mm']:.2f} mm")
    if affine_info:
        print(f"    12-DOF scales           "
              f"[{affine_info['scales'][0]:.4f} {affine_info['scales'][1]:.4f} "
              f"{affine_info['scales'][2]:.4f}]  (1.0 = correct voxel size)")

    verdict = interpret(header_dice, reg_dice, rigid)
    print(f"\n  VERDICT: {verdict}")
    print("  NOTE: mask overlap cannot resolve a left-right flip -- run "
          "lr_flip_test.py for that.")
    print("=" * 68 + "\n")

    results = {
        "label": label,
        "ref": str(args.ref),
        "test": str(args.test),
        "header_only_dice": header_dice,
        "registered_dice": reg_dice,
        "dice_gap": gap,
        "centroid_distance_mm": centroid_mm,
        "hd95_mm": hd95,
        "rigid": rigid,
        "affine": affine_info,
        "verdict": verdict,
    }

    if args.controls:
        print("  scoring 48 signed axis permutations...")
        ctl = run_controls(test_mask_img, ref_img, ref_mask)
        results["controls"] = ctl
        true_rank = next(i for i, r in enumerate(ctl) if r["identity"])
        best_wrong = next(r for r in ctl if not r["identity"])
        print(f"\n  NEGATIVE CONTROLS (header-only Dice for every wrong reorientFcn)")
        print(f"    as-written orientation ranks {true_rank + 1} of 48 "
              f"(Dice {ctl[true_rank]['dice']:.4f})")
        print(f"    best wrong orientation  Dice {best_wrong['dice']:.4f} "
              f"(perm {best_wrong['perm']}, flips {best_wrong['flips']})")
        print(f"    margin over best wrong  {ctl[true_rank]['dice'] - best_wrong['dice']:+.4f}")
        print(f"    median across all 48    {np.median([r['dice'] for r in ctl]):.4f}")
        results["control_summary"] = {
            "true_rank": true_rank + 1,
            "margin_over_best_wrong": ctl[true_rank]["dice"] - best_wrong["dice"],
            "median_dice": float(np.median([r["dice"] for r in ctl])),
        }
        print()

    print("  writing QC figure...")
    test_on_ref = resample_to(test_img, ref_img, order=1)
    qc_figure(ref_img, ref_mask, test_on_ref, test_mask_on_ref,
              work / f"qc_{label}.png", label)

    with open(work / f"report_{label}.json", "w") as fh:
        json.dump(results, fh, indent=2)

    # one row per dataset, so several runs concatenate into a study-level table
    with open(work / f"report_{label}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "header_only_dice", "registered_dice", "dice_gap",
                    "centroid_distance_mm", "hd95_mm", "rotation_deg",
                    "translation_mm", "brain_disp_mean_mm", "brain_disp_p95_mm",
                    "min_affine_scale", "control_true_rank",
                    "control_margin_over_best_wrong", "verdict"])
        cs = results.get("control_summary", {})
        w.writerow([label, f"{header_dice:.4f}", f"{reg_dice:.4f}", f"{gap:.4f}",
                    f"{centroid_mm:.2f}", f"{hd95:.2f}",
                    f"{rigid['rotation_deg']:.2f}", f"{rigid['translation_mm']:.2f}",
                    f"{rigid['brain_displacement_mean_mm']:.2f}",
                    f"{rigid['brain_displacement_p95_mm']:.2f}",
                    f"{min(affine_info['scales']):.4f}" if affine_info else "",
                    cs.get("true_rank", ""),
                    f"{cs['margin_over_best_wrong']:.4f}" if cs else "",
                    verdict])

    if args.controls:
        with open(work / f"controls_{label}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["dataset", "rank", "perm", "flips", "is_as_written", "dice"])
            for i, c in enumerate(results["controls"]):
                w.writerow([label, i + 1, "".join(map(str, c["perm"])),
                            "".join("1" if f else "0" for f in c["flips"]),
                            int(c["identity"]), f"{c['dice']:.4f}"])

    print(f"  figure : {work / f'qc_{label}.png'}")
    print(f"  report : {work / f'report_{label}.json'}")
    print(f"  csv    : {work / f'report_{label}.csv'}"
          + (f" + controls_{label}.csv" if args.controls else "") + "\n")


def interpret(header_dice, reg_dice, rigid):
    """Classify the result: correct header, header bug, or plain misalignment.

    The distinction that matters is between a header bug and genuine inter-scan
    head motion.  A bug lands on a near-exact signed axis permutation, so it
    shows up as a rotation close to a multiple of 90 deg; motion is a small
    arbitrary transform.  The gap to the registered ceiling separates the two
    when the rotation is ambiguous.
    """
    rot = rigid["rotation_deg"]
    disp = rigid["brain_displacement_mean_mm"]
    near_right_angle = min(abs(rot - k) for k in (90, 180, 270)) < 25

    if header_dice < 0.5 or near_right_angle:
        return ("FAIL - orientation looks wrong "
                f"(rotation {rot:.1f} deg near a signed axis permutation)")
    if disp > 15 or header_dice < 0.65:
        return ("INCONCLUSIVE - large residual misalignment; check whether this "
                "is real head motion between scans or an origin/FoV error")
    if reg_dice - header_dice > 0.10:
        return ("MARGINAL - orientation is broadly right but the header-only "
                "overlap sits well below the ceiling")
    return "PASS - header-only overlap is at the registration ceiling"


if __name__ == "__main__":
    main()
