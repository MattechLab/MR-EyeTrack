#!/usr/bin/env python
"""Localised orientation check at the orbit, using A-eye masks from the MPRAGE.

A-eye is trained on standard T1w and is out of distribution on fat-suppressed
LIBRE, so it is never run on the test image.  Instead it is run on the
reference MPRAGE, where it works, and the resulting masks are pushed into the
test volume through the two affines alone -- no registration.  The question is
then simply: did they land on the structure?

That is answered without segmenting the test image at all.  Each mask is given
a small local translation search, and the shift that best locks it onto the
underlying structure is reported in millimetres.  The objective is the contrast
between the mask interior and a shell around it, which is sign-agnostic: the
vitreous is bright on T2w and dark on fat-suppressed T1w, and |mean_in -
mean_shell| does not care which.

A few millimetres of displacement is expected and is not a header error -- the
MPRAGE and the reconstruction are separate acquisitions, so the head moves, and
the globe additionally moves within the orbit because gaze differs between
scans.  Compare the number against the brain-level residual from
orientation_check.py: displacement at the orbit that is comparable to the
brain-level rigid residual is motion, and only a large excess indicates a
geometry problem.

Like the brain-mask test, this cannot detect a left-right flip: it maps the
left globe onto the right globe and every measurement stays small.

Usage
-----
  python eye_roi_check.py --ref MPRAGE.nii.gz --test recon.nii.gz \
      --masks globe_ex_lens.nii.gz lens.nii.gz optic_nerve.nii.gz \
      --out WORKDIR [--search-mm 8]

Masks must be in the reference image's space (that is, A-eye run on --ref).
They may be binary or multi-label; each label, and each disconnected component
within a label, is scored separately so that the two orbits are independent
measurements.
"""

import argparse
import csv
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

from orientation_check import prepare, resample_to

MIN_COMPONENT_VOXELS = 50

# A-eye's Task313_Eye label map, from
# a-eye_web/package/biomarkers/biomarkers.py
AEYE_LABELS = {
    1: "lens", 2: "globe", 3: "nerve", 4: "int_fat", 5: "ext_fat",
    6: "lat_mus", 7: "med_mus", 8: "inf_mus", 9: "sup_mus",
}


def separability(image, idx_in, idx_shell, shift):
    """|mean_in - mean_shell| / pooled sd, for the mask translated by `shift`.

    Sign-agnostic on purpose so the same objective works for a bright vitreous
    (T2w) and a dark one (fat-suppressed T1w).
    """
    a = image[tuple(idx_in[d] + shift[d] for d in range(3))]
    b = image[tuple(idx_shell[d] + shift[d] for d in range(3))]
    pooled = np.sqrt(0.5 * (a.var() + b.var())) + 1e-9
    return float(abs(a.mean() - b.mean()) / pooled)


def local_search(image, mask, search_vox):
    """Best integer translation of `mask` within +/- search_vox, by separability.

    Returns (best_shift, score_at_header_position, score_at_best).
    """
    shell = ndimage.binary_dilation(mask, iterations=3) & ~mask
    idx_in = np.array(np.nonzero(mask))
    idx_shell = np.array(np.nonzero(shell))
    if idx_in.shape[1] < MIN_COMPONENT_VOXELS or idx_shell.shape[1] < MIN_COMPONENT_VOXELS:
        return None, np.nan, np.nan

    # keep every sampled index inside the volume for the whole search range
    lo = np.array([search_vox] * 3)
    hi = np.array(image.shape) - search_vox - 1
    keep = np.all((idx_in >= lo[:, None]) & (idx_in <= hi[:, None]), axis=0)
    idx_in = idx_in[:, keep]
    keep = np.all((idx_shell >= lo[:, None]) & (idx_shell <= hi[:, None]), axis=0)
    idx_shell = idx_shell[:, keep]
    if idx_in.shape[1] < MIN_COMPONENT_VOXELS or idx_shell.shape[1] < MIN_COMPONENT_VOXELS:
        return None, np.nan, np.nan

    rng = range(-search_vox, search_vox + 1)
    base = separability(image, idx_in, idx_shell, (0, 0, 0))
    best, best_shift = -np.inf, (0, 0, 0)
    for dx in rng:
        for dy in rng:
            for dz in rng:
                s = separability(image, idx_in, idx_shell, (dx, dy, dz))
                if s > best:
                    best, best_shift = s, (dx, dy, dz)
    return best_shift, base, best


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="reference NIfTI the masks belong to")
    p.add_argument("--test", required=True, help="NIfTI produced by mat2nii")
    p.add_argument("--masks", required=True, nargs="+", help="A-eye masks in --ref space")
    p.add_argument("--out", required=True, help="work/output dir")
    p.add_argument("--search-mm", type=float, default=8.0,
                   help="half-width of the local translation search (default 8 mm)")
    p.add_argument("--label", default=None)
    p.add_argument("--label-names", default=None,
                   help="comma-separated names for labels 1..N; defaults to the "
                        "A-eye map (lens,globe,nerve,int_fat,ext_fat,lat_mus,"
                        "med_mus,inf_mus,sup_mus)")
    args = p.parse_args()

    if args.label_names:
        label_names = {i + 1: n for i, n in enumerate(args.label_names.split(","))}
    else:
        label_names = dict(AEYE_LABELS)

    work = Path(args.out).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    label = args.label or Path(args.test).name.replace(".nii.gz", "")

    ref_img = prepare(args.ref, work / "eye_ref.nii.gz")
    test_img = prepare(args.test, work / "eye_test.nii.gz")
    test_data = np.asanyarray(test_img.dataobj).astype(np.float32)

    zooms = np.abs(test_img.header.get_zooms()[:3])
    search_vox = int(round(args.search_mm / float(np.min(zooms))))

    print(f"=== eye ROI check: {label} ===")
    print(f"  ref  : {args.ref}")
    print(f"  test : {args.test}")
    print(f"  local search +/-{args.search_mm:.1f} mm ({search_vox} voxels)\n")

    rows = []
    propagated_total = np.zeros(test_img.shape[:3], dtype=bool)

    for mask_path in args.masks:
        name = Path(mask_path).name.replace(".nii.gz", "").replace(".nii", "")
        m_img = nib.load(mask_path)
        m_data = np.asanyarray(m_img.dataobj)

        # The mask carries its own affine and is mapped through it, so a mask in
        # canonical-RAS space (as run_aeye.py produces) works unchanged; only a
        # genuinely different subject/session would be a problem.
        if not np.allclose(m_img.affine, ref_img.affine, atol=1e-3):
            print(f"  note: {name} is on a different grid from --ref "
                  f"({nib.aff2axcodes(m_img.affine)} vs {nib.aff2axcodes(ref_img.affine)}); "
                  f"mapping through its own affine")

        for lab in [v for v in np.unique(m_data) if v != 0]:
            comps, n = ndimage.label(m_data == lab)
            for c in range(1, n + 1):
                comp = comps == c
                if comp.sum() < MIN_COMPONENT_VOXELS:
                    continue

                # centroid in world coords, used only to name the side
                idx = np.array(np.nonzero(comp), dtype=np.float64).mean(axis=1)
                world = m_img.affine[:3, :3] @ idx + m_img.affine[:3, 3]
                side = "R" if world[0] > 0 else "L"
                struct = label_names.get(int(lab), f"{name}[{int(lab)}]") \
                    if len(np.unique(m_data)) > 2 else name
                tag = f"{struct}-{side}"

                # header-only propagation into the test grid
                comp_img = nib.Nifti1Image(comp.astype(np.uint8), m_img.affine)
                prop = resample_to(comp_img, test_img, order=0) > 0
                if prop.sum() < MIN_COMPONENT_VOXELS:
                    print(f"  {tag:28s} propagated outside the test FoV -- skipped")
                    continue
                propagated_total |= prop

                shift, base_score, best_score = local_search(test_data, prop, search_vox)
                if shift is None:
                    print(f"  {tag:28s} too close to the volume edge -- skipped")
                    continue

                # voxel shift -> world displacement
                disp = test_img.affine[:3, :3] @ np.array(shift, dtype=float)
                disp_mm = float(np.linalg.norm(disp))

                rows.append({
                    "structure": tag,
                    "voxels": int(prop.sum()),
                    "displacement_mm": disp_mm,
                    "displacement_ras_mm": [float(x) for x in disp],
                    "separability_at_header": base_score,
                    "separability_at_best": best_score,
                    "centroid_ref_ras": [float(x) for x in world],
                })
                print(f"  {tag:28s} displacement {disp_mm:5.2f} mm  "
                      f"RAS [{disp[0]:+5.1f} {disp[1]:+5.1f} {disp[2]:+5.1f}]  "
                      f"separability {base_score:.2f} -> {best_score:.2f}")

    if rows:
        d = np.array([r["displacement_mm"] for r in rows])
        print(f"\n  {len(rows)} structures | displacement median {np.median(d):.2f} mm, "
              f"max {d.max():.2f} mm")
        print("  Compare against the brain-level rigid residual from orientation_check.py:")
        print("  comparable = inter-scan motion; a large excess = geometry problem.")
        print("  NOTE: a left-right flip is NOT detectable here.\n")

    nib.save(nib.Nifti1Image(propagated_total.astype(np.uint8), test_img.affine),
             work / f"eye_masks_on_test_{label}.nii.gz")
    qc_figure(test_data, propagated_total, work / f"qc_eye_{label}.png", label)

    with open(work / f"report_eye_{label}.json", "w") as fh:
        json.dump({"label": label, "ref": args.ref, "test": args.test,
                   "search_mm": args.search_mm, "structures": rows}, fh, indent=2)

    with open(work / f"report_eye_{label}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "structure", "voxels", "displacement_mm",
                    "disp_R_mm", "disp_A_mm", "disp_S_mm",
                    "separability_at_header", "separability_at_best"])
        for r in rows:
            w.writerow([label, r["structure"], r["voxels"],
                        f"{r['displacement_mm']:.3f}",
                        *[f"{x:.3f}" for x in r["displacement_ras_mm"]],
                        f"{r['separability_at_header']:.4f}",
                        f"{r['separability_at_best']:.4f}"])

    print(f"  figure : {work / f'qc_eye_{label}.png'}")
    print(f"  report : {work / f'report_eye_{label}.json'}")
    print(f"  csv    : {work / f'report_eye_{label}.csv'}\n")


def qc_figure(test_data, prop, out_png, label):
    """Show the propagated masks on the test image, centred on the orbits."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not prop.any():
        return
    idx = np.array(np.nonzero(prop), dtype=np.float64).mean(axis=1).astype(int)
    hi = np.percentile(test_data, 99.5)
    img = np.clip(test_data / (hi + 1e-9), 0, 1)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, (name, ax_i) in zip(axes, [("Sagittal", 0), ("Coronal", 1), ("Axial", 2)]):
        sl = [slice(None)] * 3
        sl[ax_i] = idx[ax_i]
        sl = tuple(sl)
        ax.imshow(np.rot90(img[sl]), cmap="gray", vmin=0, vmax=1)
        ax.contour(np.rot90(prop[sl].astype(float)), levels=[0.5],
                   colors="#ffd23f", linewidths=1.2)
        ax.set_title(f"TEST - {name}", fontsize=10)
        ax.axis("off")
    fig.suptitle(f"{label}\nyellow = A-eye masks from the MPRAGE, "
                 f"propagated through the headers alone", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
