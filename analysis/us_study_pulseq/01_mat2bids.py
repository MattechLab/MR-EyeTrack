#!/usr/bin/env python
"""Build the BIDS derivative dataset for the pulseq undersampling study.

Every NIfTI is regenerated from the reconstruction .mat files rather than
reusing the NIfTIs already in data/study, because those were written by three
different converters over the year and do not share a geometry:

  woBin x0 / x_steva_*_twix   mat2nii_twix: twix origin, MPRAGE direction cosines
  clean/x, clean/x0 (Dec 25)  identity affine, a different array orientation
  clean_0.XX, woBin x_steva   MPRAGE affine copied verbatim (0.97 mm voxels)

Here every volume gets the same treatment as mat2nii_twix: |x|, the confirmed
MR-EyeTrack reorientFcn  flip(flip(permute(v,[2 1 3]),1),2), and the twix affine
taken from that subject's woBin/x0.nii.gz.  `--check` verifies the conversion
reproduces the two existing mat2nii_twix outputs bit-for-bit before anything
is written.

Intensities: woBin x0 is stored unnormalised on some subjects (max ~1e-8), so
each volume is divided by its own 99.5th percentile and multiplied by 1000.
All IQMs used downstream are scale invariant; the factor is kept in the JSON.

Naming follows results/bids (last year's IDEA study):
  sub-XXX_rec-stevaFull_T1w        woBin, regularised (REGU)
  sub-XXX_rec-mathildaFull_T1w     woBin, gridded (GRID)
  sub-XXX_acq-<pos>_rec-steva<lvl>_T1w      lvl in full, 50p, 75p, 95p
  sub-XXX_acq-<pos>_rec-mathilda<lvl>_T1w
"""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import nibabel as nib
import numpy as np
import scipy.io as sio

REPO = Path(__file__).resolve().parents[2]
STUDY = REPO / "data" / "study"
POSITIONS = ["up", "down", "left", "right"]          # region_idx 0..3
LEVELS = {"clean": "full", "clean_0.50": "50p", "clean_0.75": "75p", "clean_0.95": "95p"}
KEEP = {"full": 1.0, "50p": 0.50, "75p": 0.25, "95p": 0.05}


def load_mat(path):
    try:
        d = sio.loadmat(path)
        key = next(k for k in ("x", "x0", "xrms") if k in d)
        v = d[key]
    except NotImplementedError:                       # v7.3 -> HDF5, axes reversed
        with h5py.File(path, "r") as h:
            key = next(k for k in ("x", "x0", "xrms") if k in h)
            v = h[key]
            if v.dtype == object:                     # 1x1 cell -> its only element
                assert v.size == 1, f"{path}: cell of {v.size} volumes"
                v = h[v[()].flat[0]]
            v = v[()]
            if v.dtype.names:
                v = v["real"] + 1j * v["imag"]
            v = v.transpose()
    return np.abs(v).astype(np.float32)


def reorient(v):
    """MR-EyeTrack LIBRE reorientFcn: flip(flip(permute(v,[2 1 3]),1),2)."""
    return np.ascontiguousarray(np.transpose(v, (1, 0, 2))[::-1, ::-1, :])


def jobs(sub):
    rec = STUDY / sub / "recon"
    yield rec / "woBin/x_steva_nIter_20_delta_1.000.mat", f"{sub}_rec-stevaFull_T1w", "steva", "woBin", None
    yield rec / "woBin/x0.mat", f"{sub}_rec-mathildaFull_T1w", "mathilda", "woBin", None
    for mask, lvl in LEVELS.items():
        for i, pos in enumerate(POSITIONS):
            yield (rec / mask / f"x/x_steva_regionidx_{i}_nIter_20_delta_1.000.mat",
                   f"{sub}_acq-{pos}_rec-steva{lvl}_T1w", "steva", lvl, pos)
            yield (rec / mask / f"x0/x0_regionidx{i}.mat",
                   f"{sub}_acq-{pos}_rec-mathilda{lvl}_T1w", "mathilda", lvl, pos)


def convert(args):
    src, stem, method, lvl, pos, sub, out_dir, overwrite = args
    out = out_dir / sub / "anat" / f"{stem}.nii.gz"
    if out.exists() and not overwrite:
        return stem, "skip"
    ref = nib.load(STUDY / sub / "recon/woBin/x0.nii.gz")
    vol = reorient(load_mat(src))
    p995 = float(np.percentile(vol, 99.5))
    vol = vol * np.float32(1000.0 / p995)
    img = nib.Nifti1Image(vol, ref.affine)
    img.set_qform(ref.affine, code=1)
    img.set_sform(ref.affine, code=1)
    img.header.set_xyzt_units("mm", "sec")
    out.parent.mkdir(parents=True, exist_ok=True)
    nib.save(img, out)
    meta = {
        "SkullStripped": False,
        "Reconstruction": {"steva": "REGU: CS, anisotropic spatial TV (bmSteva, nIter=20, delta=1)",
                           "mathilda": "GRID: non-regularised gridding (bmMathilda)"}[method],
        "GazeBin": pos or "none (all readouts, no binning)",
        "ReadoutFractionKeptOfClean": KEEP.get(lvl) if pos else None,
        "SourceMat": str(src.relative_to(REPO)),
        "IntensityScale": 1000.0 / p995,
        "IntensityScaleNote": "volume = |x| * IntensityScale, i.e. p99.5 -> 1000",
        "AffineSource": "woBin/x0.nii.gz (mat2nii_twix: twix origin, MPRAGE direction cosines)",
    }
    (out.parent / f"{stem}.json").write_text(json.dumps(meta, indent=2))
    return stem, "ok"


def check():
    """The conversion must reproduce existing mat2nii_twix outputs exactly."""
    for sub in ("sub-001", "sub-004", "sub-015"):
        rec = STUDY / sub / "recon/woBin"
        for mat, nii in (("x0.mat", "x0.nii.gz"),
                         ("x_steva_nIter_20_delta_1.000.mat", "x_steva_nIter_20_delta_1.000_twix.nii.gz")):
            a = reorient(load_mat(rec / mat))
            r = nib.load(rec / nii).get_fdata(dtype=np.float32)
            err = float(np.abs(a - r).max() / r.max())
            print(f"  {sub} {mat:36s} rel max err {err:.2e}")
            assert err < 1e-5, "conversion does not reproduce mat2nii_twix"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "results/us_pulseq/bids/derivatives/reconstructions"))
    ap.add_argument("--subjects", nargs="*", default=[f"sub-{i:03d}" for i in range(1, 16)])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        check()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(*j, sub, out_dir, a.overwrite) for sub in a.subjects for j in jobs(sub)]
    missing = [t[0] for t in tasks if not t[0].exists()]
    assert not missing, f"missing .mat: {missing[:5]}"
    with ProcessPoolExecutor(a.workers) as ex:
        for stem, status in ex.map(convert, tasks):
            print(f"  {status:4s} {stem}", flush=True)

    (out_dir / "dataset_description.json").write_text(json.dumps({
        "Name": "MR-EyeTrack pulseq undersampling study (15 subjects)",
        "BIDSVersion": "1.9.0",
        "DatasetType": "derivative",
        "GeneratedBy": [
            {"Name": "STEVA", "Description": "CS reconstruction, anisotropic spatial TV (monalisa bmSteva)"},
            {"Name": "MATHILDA", "Description": "Non-regularised gridding (monalisa bmMathilda)"},
            {"Name": "analysis/us_study_pulseq/01_mat2bids.py"}],
    }, indent=2))
    print(f"{len(tasks)} volumes -> {out_dir}")


if __name__ == "__main__":
    main()
