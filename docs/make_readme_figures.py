"""Figures for the top-level README.

    python docs/make_readme_figures.py

Needs the BIDS tree built by analysis/conversion/build_bids_dataset.py (or the
OpenNeuro dataset ds008973 downloaded to the same layout).
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO = Path(__file__).resolve().parents[1]
BIDS = REPO / "data" / "study_bids"
OUT = REPO / "docs" / "img"
SUB = "sub-001"
BINS = ["up", "down", "left", "right"]

# eye-level slices of the 240^3 RAS grid (array indices), globe centre at (88, 188, 122)
GX, GY, GZ = 88, 188, 122


def load(gaze):
    f = BIDS / "derivatives/mreyetrack-recon" / SUB / "anat" / f"{SUB}_acq-{gaze}_rec-regu_T1w.nii.gz"
    return np.asarray(nib.load(f).dataobj)


def gaze_figure():
    vols = {g: load(g) for g in BINS}
    fig, ax = plt.subplots(2, 4, figsize=(10, 3.6), gridspec_kw={"wspace": 0.03, "hspace": 0.04,
                                                                  "height_ratios": [66 / 140, 55 / 75]})
    for j, g in enumerate(BINS):
        v = vols[g]
        axial = v[:, :, GZ].T
        sag = v[GX, :, :].T
        for i, (im, ext, ref) in enumerate([
            (axial[150:216, 50:190], (50, 190, 150, 216), (GX, GY)),
            (sag[95:150, 150:225], (150, 225, 95, 150), (GY, GZ)),
        ]):
            a = ax[i, j]
            a.imshow(im, origin="lower", cmap="gray", extent=ext, vmax=np.percentile(im, 99.7), aspect="equal")
            a.axvline(ref[0], c="#ff5a36", lw=0.6, ls=(0, (3, 3)))
            a.axhline(ref[1], c="#ff5a36", lw=0.6, ls=(0, (3, 3)))
            a.set_xticks([]); a.set_yticks([])
            for s in a.spines.values():
                s.set_visible(False)
        ax[0, j].set_title(f"gaze {g}", fontsize=12)
    ax[0, 0].set_ylabel("axial", fontsize=11)
    ax[1, 0].set_ylabel("sagittal", fontsize=11)
    fig.savefig(OUT / "gaze_bins_eye.png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def pipeline_figure():
    fig, ax = plt.subplots(figsize=(10, 3.1))
    ax.set_xlim(0, 100); ax.set_ylim(0, 31); ax.axis("off")
    blue, orange, green, grey = "#2f6fb5", "#e0702d", "#3c9a63", "#444444"

    def box(x, y, w, h, text, c):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    fc="white", ec=c, lw=1.6))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.6, color=grey)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, lw=1.3, color=grey))

    box(1, 19, 21, 9, "3D radial-spiral k-space\n(pulseq LIBRE, TR 8 ms)", blue)
    box(1, 3, 21, 9, "Eye tracking\n(EyeLink, 1 kHz)", orange)
    box(30, 11, 21, 9, "Sync + sliding-window\ngaze label per readout", green)
    box(59, 11, 18, 9, "k-space split by gaze\nup / down / left / right", green)
    box(85, 11, 14, 9, "4 images\n(one per gaze)", blue)
    arrow(22.6, 23.5, 29.2, 17); arrow(22.6, 7.5, 29.2, 14)
    arrow(51.6, 15.5, 58.2, 15.5); arrow(77.6, 15.5, 84.2, 15.5)
    for x, t in [(40.5, "S2"), (68, "S3"), (92, "S4")]:
        ax.text(x, 22.5 if x != 40.5 else 22.5, t, ha="center", fontsize=9, color=green, fontweight="bold")
    fig.savefig(OUT / "pipeline.png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    pipeline_figure()
    gaze_figure()
