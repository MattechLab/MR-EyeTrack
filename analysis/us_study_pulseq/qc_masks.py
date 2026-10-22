#!/usr/bin/env python
"""QC figure per subject: A-eye labels (after MPRAGE->LIBRE) and brain dseg on
woBin REGU and on the four clean gaze bins, axial through the lens centroid.

    python qc_masks.py [--subjects sub-001 ...]
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
COLORS = {1: "yellow", 2: "cyan", 3: "magenta", 4: "orange", 5: "orange",
          6: "lime", 7: "lime", 8: "lime", 9: "lime"}


def contour(ax, lab2d, colors):
    for v, c in colors.items():
        m = lab2d == v
        if m.any():
            ax.contour(m.T, levels=[0.5], colors=c, linewidths=0.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", nargs="*", default=[f"sub-{i:03d}" for i in range(1, 16)])
    a = ap.parse_args()
    out = ROOT / "qc"
    out.mkdir(exist_ok=True)
    for sub in a.subjects:
        m = ROOT / "masks" / sub
        L = np.asanyarray(nib.load(m / "aeye_labels.nii.gz").dataobj)
        dseg_p = m / "brain_dseg.nii.gz"
        dseg = np.asanyarray(nib.load(dseg_p).dataobj) if dseg_p.exists() else None
        z = int(round(ndi.center_of_mass(L == 1)[2])) if (L == 1).any() else L.shape[2] // 2
        xs = np.nonzero(L.any(axis=(1, 2)))[0]
        ys = np.nonzero(L.any(axis=(0, 2)))[0]
        x0, x1 = max(xs.min() - 15, 0), min(xs.max() + 15, L.shape[0])
        y0, y1 = max(ys.min() - 10, 0), min(ys.max() + 10, L.shape[1])
        vols = [("woBin", BIDS / sub / "anat" / f"{sub}_rec-stevaFull_T1w.nii.gz")] + [
            (p, BIDS / sub / "anat" / f"{sub}_acq-{p}_rec-stevafull_T1w.nii.gz") for p in ["up", "down", "left", "right"]]
        fig, axs = plt.subplots(2, 5, figsize=(17, 7))
        for j, (name, p) in enumerate(vols):
            v = nib.load(p).get_fdata(dtype=np.float32)
            sl = v[x0:x1, y0:y1, z]
            for i in range(2):
                ax = axs[i, j]
                ax.imshow(sl.T, cmap="gray", origin="lower", vmin=0, vmax=np.percentile(sl, 99.5))
                if i == 1:
                    contour(ax, L[x0:x1, y0:y1, z], COLORS)
                ax.set_title(f"{name}" + ("  + A-eye" if i else ""), fontsize=9)
                ax.axis("off")
        fig.suptitle(f"{sub}  axial z={z} (lens centroid)  yellow lens / cyan vitreous / magenta nerve / "
                     f"orange fat / lime recti", fontsize=10)
        fig.tight_layout()
        fig.savefig(out / f"{sub}_aeye_axial.png", dpi=110)
        plt.close(fig)

        if dseg is not None:
            v = nib.load(vols[0][1]).get_fdata(dtype=np.float32)
            fig, axs = plt.subplots(1, 3, figsize=(13, 4.5))
            for ax, (axis, idx) in zip(axs, [(2, 120), (1, 120), (0, 120)]):
                sl = np.take(v, idx, axis=axis)
                ax.imshow(sl.T, cmap="gray", origin="lower", vmin=0, vmax=np.percentile(v, 99.5))
                ax.contour(np.take(dseg, idx, axis=axis).T, levels=[0.5, 1.5, 2.5], colors=["r", "g", "b"],
                           linewidths=0.5)
                eyes = ROOT / "masks" / sub / "eyes_mni.nii.gz"
                if eyes.exists():
                    ax.contour(np.take(np.asanyarray(nib.load(eyes).dataobj), idx, axis=axis).T,
                               levels=[0.5], colors="y", linewidths=0.6)
                ax.axis("off")
            fig.suptitle(f"{sub} woBin REGU: MRIQC dseg (r CSF, g GM, b WM), MNI eye mask (yellow)")
            fig.tight_layout()
            fig.savefig(out / f"{sub}_brain.png", dpi=100)
            plt.close(fig)
        print(sub, "qc ok", flush=True)


if __name__ == "__main__":
    main()
