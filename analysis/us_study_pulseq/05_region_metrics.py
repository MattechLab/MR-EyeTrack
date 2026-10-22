#!/usr/bin/env python
"""Fixed-mask image quality metrics for every volume of the undersampling study.

One row per (subject, method, gaze position, level, region, metric). Masks come
from 04_register_masks.py and are the same for every volume of a subject, so a
metric change across levels is a change in the image, not in the mask.

Regions
  volume                      whole 240^3 FoV (last year's "Whole Volume")
  brain, csf, gm, wm          MRIQC segmentation of woBin REGU
  eyes_mni                    MNI152NLin2009cAsym desc-eye mask, last year's "Eyes".
                              NB: it is not an eye mask -- 1.17 L, the whole
                              anterior-inferior head block below z = 1 mm and in
                              front of y = -4 mm (orbits, nose, sinuses, face and
                              air in front of it), minus the brain.
  lens, vitreous, eyeball,    A-eye labels from the MPRAGE (eyeball = lens +
  nerve, fat, muscle, orbit   vitreous, fat = intra+extraconal, muscle = 4 recti)

Metrics
  ssim        mean of the 3D SSIM map (Gaussian sigma 1.5, Wang 2004) inside the
              region, against the clean (fully binned, no undersampling) volume
              of the same method and gaze bin. Last year's notebook ran skimage
              SSIM on the 1D vector of masked voxels, whose 7-sample window
              follows array order rather than anatomy; this replaces it.
  nrmse       ||x - ref|| / ||ref|| inside the region, same reference.
  ssim_wobin, nrmse_wobin   the same against woBin (all 79 894 readouts): one
              reference shared by every bin and level of a subject.
  dcte        DCT energy ratio exactly as analysis/comparison/dcte/dcte.ipynb
              (min-max normalised volume, zeroed outside the region, high/low
              octant ratio over mean intensity) -- kept for continuity.
  mean, std, cv, p05, p95     intensity statistics (volume / p99.5)
  contrast and CNR between adjacent tissues, edge width at the lens and globe
  boundaries (see tissue_metrics).

Intensities: each volume is divided by its own 99.5th percentile (01_mat2bids
stores p99.5 = 1000), so SSIM uses data_range 1 on [0, 1].
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.fft import dctn
from skimage.metrics import structural_similarity

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
POSITIONS = ["up", "down", "left", "right"]
LEVELS = ["full", "50p", "75p", "95p"]
METHODS = {"REGU": "steva", "GRID": "mathilda"}
DCTE_REGIONS = ["volume", "brain", "csf", "gm", "wm", "eyes_mni", "eyeball", "orbit"]


def path(sub, method, pos, level):
    tag = METHODS[method]
    if level == "woBin":
        return BIDS / sub / "anat" / f"{sub}_rec-{tag}Full_T1w.nii.gz"
    return BIDS / sub / "anat" / f"{sub}_acq-{pos}_rec-{tag}{level}_T1w.nii.gz"


def load(p):
    return nib.load(p).get_fdata(dtype=np.float32)


def regions(sub):
    m = ROOT / "masks" / sub
    dseg = np.asanyarray(nib.load(m / "brain_dseg.nii.gz").dataobj)
    eyes = np.asanyarray(nib.load(m / "eyes_mni.nii.gz").dataobj) > 0
    L = np.asanyarray(nib.load(m / "aeye_labels.nii.gz").dataobj)
    return {
        "volume": np.ones(dseg.shape, bool),
        "brain": dseg > 0, "csf": dseg == 1, "gm": dseg == 2, "wm": dseg == 3,
        "eyes_mni": eyes,
        "lens": L == 1, "vitreous": L == 2, "eyeball": (L == 1) | (L == 2),
        "nerve": L == 3, "fat": (L == 4) | (L == 5), "muscle": (L >= 6) & (L <= 9),
        "orbit": L > 0,
    }


def dcte_legacy(vol, mask, hf_ratio=0.5):
    """analysis/comparison/dcte/dcte.ipynb compute_dcte_3d, on min-max data."""
    d = (vol - vol.min()) / (vol.max() - vol.min() + 1e-8)
    d = np.where(mask, d, 0).astype(np.float64)
    c = np.abs(dctn(d, norm="ortho", workers=2))
    nx, ny, nz = c.shape
    cx, cy, cz = int(nx * (1 - hf_ratio)), int(ny * (1 - hf_ratio)), int(nz * (1 - hf_ratio))
    low, high = c[:cx, :cy, :cz].mean(), c[cx:, cy:, cz:].mean()
    return (high / (low + 1e-8)) / (d.mean() + 1e-8)


def signed_distance(mask):
    """mm from the mask surface, negative inside (1 mm isotropic grid)."""
    return ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)


def edge_width(v, sd, lo, hi, step=0.5):
    """Edge width (mm) from the boundary profile averaged over the whole surface.

    Voxels are binned by signed distance to the structure surface (lo..hi mm,
    `step` bins) and averaged, which removes noise without smoothing the image.
    Width = profile range / steepest profile slope: a Gaussian-blurred step of
    FWHM w gives ~w, and it tolerates the non-monotone globe wall (vitreous ->
    bright choroid/sclera -> dark suppressed fat) where an erf fit would not."""
    sel = (sd >= lo) & (sd <= hi)
    edges = np.arange(lo, hi + step, step)
    k = np.digitize(sd[sel], edges) - 1
    n = np.bincount(k, minlength=len(edges))
    prof = np.bincount(k, weights=v[sel], minlength=len(edges)) / np.maximum(n, 1)
    ok = n[:-1] > 20
    d, prof = (edges[:-1] + step / 2)[ok], prof[:-1][ok]
    slope = np.abs(np.gradient(prof, d)).max()
    return float((prof.max() - prof.min()) / slope) if slope > 0 else np.nan


def tissue_metrics(v, R, cache):
    out = {}
    for name, m in R.items():
        if not m.any():
            continue
        x = v[m]
        mu, sd = float(x.mean()), float(x.std())
        out[(name, "mean")] = mu
        out[(name, "std")] = sd
        out[(name, "cv")] = sd / mu if mu > 0 else np.nan
        out[(name, "p05")], out[(name, "p95")] = map(float, np.percentile(x, [5, 95]))

    def pair(a, b, tag):
        if not (R[a].any() and R[b].any()):
            return
        xa, xb = v[R[a]], v[R[b]]
        d = xa.mean() - xb.mean()
        out[(tag, "contrast")] = float(d / (xa.mean() + xb.mean()))
        out[(tag, "cnr")] = float(abs(d) / np.sqrt(xa.var() + xb.var()))
        out[(tag, "cjv")] = float((xa.std() + xb.std()) / abs(d)) if d else np.nan

    pair("lens", "vitreous", "lens_vitreous")
    pair("nerve", "fat", "nerve_fat")
    pair("muscle", "fat", "muscle_fat")
    pair("wm", "gm", "wm_gm")

    if R["lens"].any():
        out[("lens", "edge_width")] = edge_width(v, cache["sd_lens"], -2.0, 3.0)
    if R["eyeball"].any():
        out[("eyeball", "edge_width")] = edge_width(v, cache["sd_eyeball"], -4.0, 4.0)
    return out


def process(args):
    sub, method = args
    R = regions(sub)
    cache = {"sd_lens": signed_distance(R["lens"]), "sd_eyeball": signed_distance(R["eyeball"])}
    rows = []
    vols = {}
    for pos in POSITIONS:
        for lvl in LEVELS:
            vols[(pos, lvl)] = load(path(sub, method, pos, lvl))
    vols[("none", "woBin")] = load(path(sub, method, None, "woBin"))

    for (pos, lvl), raw in vols.items():
        v = np.clip(raw / 1000.0, 0, 1)
        base = dict(subject=sub, method=method, position=pos, level=lvl)
        for (region, metric), val in tissue_metrics(v, R, cache).items():
            rows.append({**base, "region": region, "metric": metric, "value": val})
        for region in DCTE_REGIONS:
            if R[region].any():
                rows.append({**base, "region": region, "metric": "dcte", "value": dcte_legacy(raw, R[region])})

    # full-reference metrics. "ssim": every level of a bin, and woBin, against that
    # bin's clean volume (last year's design). "ssim_wobin": every binned volume
    # against woBin (79 894 readouts) -- a reference shared by all bins of a
    # subject, exact for static anatomy (brain), motion-confounded at the eye.
    def full_ref(ref, t, base, suffix):
        _, smap = structural_similarity(ref, t, data_range=1.0, gaussian_weights=True, sigma=1.5,
                                        use_sample_covariance=False, full=True)
        for region, m in R.items():
            if not m.any():
                continue
            rows.append({**base, "region": region, "metric": f"ssim{suffix}", "value": float(smap[m].mean())})
            rows.append({**base, "region": region, "metric": f"nrmse{suffix}",
                         "value": float(np.linalg.norm(t[m] - ref[m]) / np.linalg.norm(ref[m]))})

    wobin = np.clip(vols[("none", "woBin")] / 1000.0, 0, 1)
    for pos in POSITIONS:
        ref = np.clip(vols[(pos, "full")] / 1000.0, 0, 1)
        for lvl in ["50p", "75p", "95p", "woBin"]:
            t = wobin if lvl == "woBin" else np.clip(vols[(pos, lvl)] / 1000.0, 0, 1)
            full_ref(ref, t, dict(subject=sub, method=method, position=pos, level=lvl), "")
        for lvl in LEVELS:
            t = np.clip(vols[(pos, lvl)] / 1000.0, 0, 1)
            full_ref(wobin, t, dict(subject=sub, method=method, position=pos, level=lvl), "_wobin")
    print(f"{sub} {method} done", flush=True)
    return rows


def run_cached(task):
    sub, method = task
    out = ROOT / "metrics/region" / f"{sub}_{method}.csv"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(process(task)).to_csv(out, index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", nargs="*", default=[f"sub-{i:03d}" for i in range(1, 16)])
    ap.add_argument("--methods", nargs="*", default=list(METHODS))
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--out", default=str(ROOT / "metrics/region_metrics.csv"))
    a = ap.parse_args()
    tasks = [(s, m) for s in a.subjects for m in a.methods
             if (ROOT / "masks" / s / "aeye_labels.nii.gz").exists()]
    with ProcessPoolExecutor(a.workers) as ex:
        list(ex.map(run_cached, tasks))
    # one cache file per subject x method; the table is rebuilt from all of them
    df = pd.concat([pd.read_csv(f) for f in sorted((ROOT / "metrics/region").glob("*.csv"))],
                   ignore_index=True)
    df.to_csv(a.out, index=False)
    print(f"{len(df)} rows from {df.subject.nunique()} subjects -> {a.out}")


if __name__ == "__main__":
    main()
