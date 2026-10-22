#!/usr/bin/env python
"""DCTE on a box around each structure, without zero-filling.

The legacy DCTE (05_region_metrics.py, from analysis/comparison/dcte/dcte.ipynb)
zeroes everything outside the region before the DCT. The step at the mask border
then dominates the high-frequency energy: filling the WM mask of sub-001 with a
constant gives DCTE 21.6 against 21.2 for the real image (GM 54.5 vs 55.5), so the
tissue-level values mostly measure the mask shape, which is the same at every
undersampling level.

Here the same ratio (high/low DCT octant, over the mean intensity) is taken on the
image inside the bounding box of the structure (+3 voxels), with nothing zeroed:
  eyeball  mean over the two eyes, each its own box around lens + vitreous (A-eye)
  brain    box around the MRIQC brain segmentation
Higher = more high-frequency energy = more in focus.

Output: results/us_pulseq/metrics/dcte_crop.csv (one row per volume x region).
"""
import runpy
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.fft import dctn

RM = runpy.run_path(str(Path(__file__).with_name("05_region_metrics.py")), run_name="dcte_crop")
ROOT = RM["ROOT"]
POSITIONS, LEVELS = RM["POSITIONS"], RM["LEVELS"]


def dcte_box(v):
    d = (v - v.min()) / (v.max() - v.min() + 1e-8)
    c = np.abs(dctn(d.astype(np.float64), norm="ortho"))
    h = [s // 2 for s in c.shape]
    return (c[h[0]:, h[1]:, h[2]:].mean() / (c[:h[0], :h[1], :h[2]].mean() + 1e-8)) / (d.mean() + 1e-8)


def boxes(sub, pad=3):
    m = ROOT / "masks" / sub
    L = np.asanyarray(nib.load(m / "aeye_labels.nii.gz").dataobj)
    dseg = np.asanyarray(nib.load(m / "brain_dseg.nii.gz").dataobj)
    lab, n = ndi.label((L == 1) | (L == 2))
    sizes = ndi.sum(lab > 0, lab, range(1, n + 1))
    out = {f"eye{i}": ndi.find_objects((lab == k + 1).astype(int))[0]
           for i, k in enumerate(np.argsort(sizes)[::-1][:2])}
    out["brain"] = ndi.find_objects((dseg > 0).astype(int))[0]
    grow = lambda sl: tuple(slice(max(s.start - pad, 0), min(s.stop + pad, n)) for s, n in zip(sl, L.shape))
    return {k: grow(v) for k, v in out.items()}


def run(task):
    sub, method = task
    B = boxes(sub)
    rows = []
    items = [(p, l) for p in POSITIONS for l in LEVELS] + [("none", "woBin")]
    for pos, lvl in items:
        v = RM["load"](RM["path"](sub, method, None if lvl == "woBin" else pos, lvl))
        base = dict(subject=sub, method=method, position=pos, level=lvl, metric="dcte_crop")
        rows.append({**base, "region": "eyeball", "value": np.mean([dcte_box(v[B["eye0"]]), dcte_box(v[B["eye1"]])])})
        rows.append({**base, "region": "brain", "value": dcte_box(v[B["brain"]])})
    return rows


def main():
    tasks = [(f"sub-{i:03d}", m) for i in range(1, 16) for m in RM["METHODS"]]
    with ProcessPoolExecutor(10) as ex:
        rows = [r for chunk in ex.map(run, tasks) for r in chunk]
    out = ROOT / "metrics/dcte_crop.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
