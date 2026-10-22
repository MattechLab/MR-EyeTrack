#!/usr/bin/env python
"""MReyeQC eye IQMs (and, optionally, its trained rating model) on every volume.

Runs the IQM engine of /home/debi/jaime/repos/mreyeqc (fetmrqc_sr.SRMetrics)
unchanged: robust preprocessing on, bias correction off (its CLI defaults),
image resampled to 0.8 mm and cropped to the eye labels. The A-eye labels are
the LIBRE-space ones from 04_register_masks.py, split per orbit because A-eye
(Task313) and, as far as we know, the MReyeQC training data are single-orbit.

Label bug (mreyeqc @ 8420cfb): SEGM numbers LENS as 0 and load_and_format_seg
remaps raw label 1 -> 0, so seg_dict["LENS"] = (seg == 0) is the lens PLUS the
whole background. Every lens IQM, and seg_cnr / seg_cjv (lens vs globe), is
computed on background. `--labels fixed` (default) patches SEGM/map_seg to
1-based in this process only; `--labels as-trained` keeps the bug, which is
what the shipped regression model saw in training, and is the only input the
model can be fed consistently.

Needs the `mreyeqc` conda env (python 3.9 + fetal_brain_qc) and the mreyeqc
repo on PYTHONPATH -- nothing is installed into that repo:

    PYTHONPATH=/home/debi/jaime/repos/mreyeqc \
      /home/debi/miniconda3/envs/mreyeqc/bin/python 06_mreyeqc_iqms.py iqms
    ... 06_mreyeqc_iqms.py predict      # after MRIQC: adds the model's rating
"""
import argparse
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
MREYEQC = Path("/home/debi/jaime/repos/mreyeqc")
MODEL_DIR = MREYEQC / "analysis/trained_models_MREye_REGRESSION"
POSITIONS = ["up", "down", "left", "right"]
LEVELS = ["full", "50p", "75p", "95p"]
METHODS = {"REGU": "steva", "GRID": "mathilda"}


def split_orbits(sub):
    """Write aeye_labels_{L,R}.nii.gz: one orbit each, split at the midpoint of
    the two vitreous centroids along world x (RAS: +x is the subject's right)."""
    m = ROOT / "masks" / sub
    out = {s: m / f"aeye_labels_{s}.nii.gz" for s in "LR"}
    if all(p.exists() for p in out.values()):
        return out
    im = nib.load(m / "aeye_labels.nii.gz")
    L = np.asanyarray(im.dataobj).astype(np.uint8)
    lab, n = ndi.label(L == 2)
    sizes = ndi.sum(L == 2, lab, range(1, n + 1))
    cents = [ndi.center_of_mass(L == 2, lab, int(i) + 1) for i in np.argsort(sizes)[::-1][:2]]
    wx = [float(im.affine[0, :3] @ np.array(c) + im.affine[0, 3]) for c in cents]
    mid = np.mean(wx)
    ijk = np.indices(L.shape, sparse=True)
    world_x = im.affine[0, 0] * ijk[0] + im.affine[0, 1] * ijk[1] + im.affine[0, 2] * ijk[2] + im.affine[0, 3]
    for s, keep in (("R", world_x > mid), ("L", world_x <= mid)):
        nib.save(nib.Nifti1Image(np.where(keep, L, 0).astype(np.uint8), im.affine), out[s])
    return out


def images(sub):
    for method, tag in METHODS.items():
        yield method, "none", "woBin", BIDS / sub / "anat" / f"{sub}_rec-{tag}Full_T1w.nii.gz"
        for pos in POSITIONS:
            for lvl in LEVELS:
                yield method, pos, lvl, BIDS / sub / "anat" / f"{sub}_acq-{pos}_rec-{tag}{lvl}_T1w.nii.gz"


FIXED_MAP = [None, 1, 2, 3, 4, 4, 5, 5, 5, 5]
SEGM_AS_SHIPPED = {"LENS": 0, "GLOBE": 1, "OPTIQUE_NERVE": 2, "FAT": 3, "MUSCLE": 4}


def run_one(task):
    sub, method, pos, lvl, img, side, seg, labels = task
    warnings.filterwarnings("ignore")
    from fetmrqc_sr.metrics import SRMetrics
    from fetmrqc_sr.metrics import metrics_sr
    kw = {}
    if labels == "fixed":
        metrics_sr.SEGM = {k: v + 1 for k, v in SEGM_AS_SHIPPED.items()}
        kw["map_seg"] = FIXED_MAP
    else:
        metrics_sr.SEGM = dict(SEGM_AS_SHIPPED)
    devnull = open(os.devnull, "w")
    old = sys.stdout
    sys.stdout = devnull                     # the engine prints per metric
    try:
        res = SRMetrics(verbose=False, robust_preprocessing=True, correct_bias=False, **kw).evaluate_metrics(
            str(img), str(seg), str(seg))
    finally:
        sys.stdout = old
        devnull.close()
    return {"bids_name": img.name.replace(".nii.gz", ""), "subject": sub, "method": method,
            "position": pos, "level": lvl, "side": side, **res}


def iqms(a):
    tasks = []
    for sub in a.subjects:
        segs = split_orbits(sub)
        for method, pos, lvl, img in images(sub):
            for side, seg in segs.items():
                tasks.append((sub, method, pos, lvl, img, side, seg, a.labels))
    out = Path(a.out)
    done = set()
    if out.exists():
        prev = pd.read_csv(out)
        done = set(zip(prev.bids_name, prev.side))
    tasks = [t for t in tasks if (t[4].name.replace(".nii.gz", ""), t[5]) not in done]
    print(f"{len(tasks)} evaluations to run", flush=True)
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(run_one, tasks, chunksize=2)):
            rows.append(r)
            if (i + 1) % 50 == 0 or i + 1 == len(tasks):
                df = pd.DataFrame(rows)
                if out.exists() and done:
                    df = pd.concat([pd.read_csv(out), df], ignore_index=True).drop_duplicates(["bids_name", "side"])
                df.to_csv(out, index=False)
                print(f"  {i + 1}/{len(tasks)}", flush=True)


def predict(a):
    """Apply the trained MReyeQC regressor (rating; higher = better).

    The model also expects MRIQC whole-image IQMs (cjv, cnr, ... wm2max); they
    come from this study's MRIQC group tables. Out of distribution: the model
    was trained on non-fat-suppressed T1w, so treat the output as exploratory."""
    import joblib
    feats = json.loads(next(MODEL_DIR.glob("features_*.json")).read_text())
    model = joblib.load(next(MODEL_DIR.glob("model_*.joblib")))
    scaler = joblib.load(next(MODEL_DIR.glob("scaler_*.joblib")))
    eye = pd.read_csv(a.out)
    mq = pd.concat([pd.read_csv(ROOT / f"mriqc_{t}/group_T1w.tsv", sep="\t") for t in ("regu", "grid", "wobin")
                    if (ROOT / f"mriqc_{t}/group_T1w.tsv").exists()]).drop_duplicates("bids_name")
    df = eye.merge(mq, on="bids_name", how="inner", suffixes=("", "_mriqc"))   # no MRIQC row -> no rating
    missing = [f for f in feats if f not in df.columns]
    if missing:
        print(f"WARNING: {len(missing)} model features absent, set to 0: {missing[:10]}")
        for f in missing:
            df[f] = 0.0
    X = df[feats].astype(float).fillna(0.0)
    X = X.replace({True: 1.0, False: 0.0})
    df["mreyeqc_rating"] = model.predict(scaler.transform(X.values))
    keep = ["bids_name", "subject", "method", "position", "level", "side", "mreyeqc_rating"]
    df[keep].to_csv(Path(a.out).with_name(f"mreyeqc_ratings_{a.labels.replace('-', '')}.csv"), index=False)
    print(df.groupby(["method", "level"]).mreyeqc_rating.describe())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["iqms", "predict"])
    ap.add_argument("--subjects", nargs="*", default=[f"sub-{i:03d}" for i in range(1, 16)])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--labels", choices=["fixed", "as-trained"], default="fixed")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    a.out = a.out or str(ROOT / f"metrics/mreyeqc_iqms_{a.labels.replace('-', '')}.csv")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    {"iqms": iqms, "predict": predict}[a.step](a)


if __name__ == "__main__":
    main()
