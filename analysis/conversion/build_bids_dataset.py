"""Assemble data/study into a BIDS dataset (staging tree for a Zenodo upload).

    python analysis/conversion/build_bids_dataset.py [--mode symlink|hardlink|copy]
                                                     [--subjects 1 2 ...] [--dry-run]
                                                    

Layout (BIDS 1.9: sourcedata/ and derivatives/ live at the dataset root)

    study_bids/
      dataset_description.json  participants.tsv  README.md  CHANGES
      sub-NNN/anat/   sub-NNN_acq-full_rec-regu_T1w.nii.gz   <- LIBRE, no binning (main image)
                      sub-NNN_acq-mprage_T1w.nii.gz            <- IDEA/Siemens reference MPRAGE (sub-010/011: none)
      sourcedata/sub-NNN/{rawdata,et}         raw twix .dat, EDF/psydat/log/csv + BIDS-converted ET
      sourcedata/pulseq/                      the two .seq files (shared by all subjects)
      derivatives/mreyetrack-recon/sub-NNN/anat/
                      sub-NNN_acq-{up,down,left,right}_rec-{regu,grid}_T1w      (clean)
                      sub-NNN_acq-full_rec-grid_T1w                               (no binning)
                      ..._desc-keep{50,25,5}_T1w                                      (random undersampling)
      derivatives/mreyetrack-binning/sub-NNN/{eyemasks,bins}   gaze labels per readout

`acq` = which gaze bin (or `full` = all readouts), `rec` = algorithm (regu = regularised iterative
compressed-sensing/STEVA, grid = gridding/MATHILDA). `desc-keepNN` = NN % of the clean readouts kept.
"""
import argparse
import csv
import json
import os
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "data" / "study"
DST = REPO / "data" / "study_bids"

GAZE = {0: "up", 1: "down", 2: "left", 3: "right"}
# mask_type -> (desc label or None, fraction of the clean readouts that is KEPT)
MASKS = {"clean": (None, 1.0), "clean_0.50": ("keep50", 0.50),
         "clean_0.75": ("keep25", 0.25), "clean_0.95": ("keep5", 0.05)}
RECS = {  # rec label -> (subfolder of recon/<mask>, filename pattern, long name)
    "regu": ("x", "x_steva_regionidx_{r}_nIter_20_delta_1.000.nii.gz", "STEVA"),
    "grid": ("x0", "x0_regionidx{r}.nii.gz", "MATHILDA"),
}
SEQ_MAIN = "yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq"

LICENSE = "CC0"
BIDS_VERSION = "1.9.0"

AUTHORS = ["Jaime Barranco", "Yiwei Jia", "Eleonora Fornari", "Oscar Esteban", "Jean-Baptiste Ledoux",
           "Mauro Leidi", "Helene Vitali", "Cédric Campos Carvalho", "Elena Najdenovska",
           "Laura Elena Raileanu", "Jessica Bastiaansen", "Micah M. Murray", "Benedetta Franceschiello"]
FUNDING = [
    "Fondation Asile des aveugles (grant #232933 to M.M.M.)",
    "A grantor advised by Carigest SA (#232920 to M.M.M.)",
    "Swiss National Science Foundation (grants #169206 to M.M.M., #220433 and #229214 to B.F.)",
    "HES-SO University of Applied Sciences and Arts Western Switzerland (grant 129816/IA-RECHERCHE23-19)",
]
ETHICS = ["MR-Eye Track protocol: positive approval from the local ethics committee via the Swiss BASEC "
          "(Business Administration System for Ethics Committees) on 27 June 2024, protocol number 2024-00165",
          "Original protocol (Franceschiello et al.): approved by the cantonal ethics committee, "
          "protocol number 2018-00240"]

ACQ = {  # shared by every LIBRE image
    "Manufacturer": "Siemens", "ManufacturersModelName": "MAGNETOM Prisma Fit",
    "SoftwareVersions": "syngo MR XA60", "MagneticFieldStrength": 3,
    "MRAcquisitionType": "3D", "PulseSequenceType": "Pulseq 3D radial-spiral LIBRE (T1w)",
    "RepetitionTime": 0.008, "EchoTime": 0.0036, "FlipAngle": 6,
    "SequenceName": SEQ_MAIN, "NumberOfSegmentsPerShot": 44, "NumberOfShots": 1872,
    "WarmUpShotsDiscarded": 14, "ReconMatrixSize": [240, 240, 240],
}


class Linker:
    def __init__(self, mode, dry):
        self.mode, self.dry, self.n, self.missing = mode, dry, 0, []

    def put(self, src, dst):
        src, dst = Path(src), Path(dst)
        if not src.exists():
            self.missing.append(str(src.relative_to(SRC)) if SRC in src.parents else str(src))
            return False
        self.n += 1
        if self.dry:
            return True
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        if self.mode == "symlink":
            dst.symlink_to(src.resolve())
        elif self.mode == "hardlink":
            os.link(src, dst)
        else:
            shutil.copy2(src, dst)
        return True

    def put_3d(self, src, dst):
        """dcm2niix writes some volumes as (x, y, z, 1); BIDS T1w must be 3D, so those are
        rewritten (small files) rather than linked."""
        import nibabel as nib
        if not Path(src).exists():
            return self.put(src, dst)
        img = nib.load(src)
        if img.ndim == 3:
            return self.put(src, dst)
        self.n += 1
        if not self.dry:
            dst = Path(dst)
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            out = nib.Nifti1Image(img.get_fdata().reshape(img.shape[:3]).astype(img.get_data_dtype()),
                                  img.affine, img.header)
            out.header.set_data_dtype(img.get_data_dtype())
            nib.save(out, dst)
        return True

    def tree(self, src, dst):
        for p in sorted(Path(src).rglob("*")):
            if p.is_file():
                self.put(p, Path(dst) / p.relative_to(src))

    def write_json(self, dst, obj):
        if self.dry:
            return
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        Path(dst).write_text(json.dumps(obj, indent=2) + "\n")

    def write_text(self, dst, text):
        if self.dry:
            return
        Path(dst).parent.mkdir(parents=True, exist_ok=True)
        Path(dst).write_text(text)


# (age, sex) read from the twix headers of each subject's own scan (flPatientAge, lPatientSex;
# Siemens codes 1 = female, 2 = male).
DEMOGRAPHICS = {1: (30, "M"), 2: (26, "F"), 3: (44, "M"), 4: (42, "M"), 5: (32, "M"), 6: (40, "M"),
                7: (28, "M"), 8: (40, "M"), 9: (25, "M"), 10: (21, "F"), 11: (28, "M"), 12: (27, "M"),
                13: (27, "M"), 14: (30, "M"), 15: (24, "F")}
BAD_DICOM = {10, 11}  # their dicom/ folder is a byte-identical copy of sub-001's, not their own scan


def build_subject(L, sub):
    s = SRC / sub
    anat = DST / sub / "anat"
    der = DST / "derivatives" / "mreyetrack-recon" / sub / "anat"

    # --- main image: LIBRE without binning, STEVA
    stem = f"{sub}_acq-full_rec-regu_T1w"
    if L.put(s / "recon/woBin/x_steva_nIter_20_delta_1.000.nii.gz", anat / f"{stem}.nii.gz"):
        L.write_json(anat / f"{stem}.json", {**ACQ, "ReconstructionAlgorithm": "STEVA",
                     "GazeBin": "none (all fixation readouts)", "NumberOfIterations": 20, "Delta": 1.0,
                     "SourceDatasets": [{"URL": f"sourcedata/{sub}/rawdata/{sub}_T1wLIBRE.dat"}]})

    stem = f"{sub}_acq-full_rec-grid_T1w"
    if L.put(s / "recon/woBin/x0.nii.gz", anat / f"{stem}.nii.gz"):
        L.write_json(anat / f"{stem}.json", {**ACQ, "ReconstructionAlgorithm": "MATHILDA (gridding)",
                     "GazeBin": "none (all fixation readouts)",
                     "SourceDatasets": [{"URL": f"sourcedata/{sub}/rawdata/{sub}_T1wLIBRE.dat"}]})

    # --- reference MPRAGE converted from the IDEA DICOM
    num = int(sub.split("-")[1])
    mp = next(iter(sorted((s / "dicom").glob("*/*.nii.gz"))), None)
    if mp is not None and num not in BAD_DICOM:
        stem = f"{sub}_acq-mprage_T1w"
        L.put_3d(mp, anat / f"{stem}.nii.gz")
        L.put(mp.with_suffix("").with_suffix(".json"), anat / f"{stem}.json")

    # --- sourcedata
    sd = DST / "sourcedata" / sub
    L.tree(s / "rawdata", sd / "rawdata")
    L.tree(s / "et", sd / "et")

    # --- derivatives: other reconstructions
    for mask, (desc, kept) in MASKS.items():
        for rec, (folder, pat, algo) in RECS.items():
            for r, gaze in GAZE.items():
                stem = f"{sub}_acq-{gaze}_rec-{rec}" + (f"_desc-{desc}" if desc else "") + "_T1w"
                if L.put(s / "recon" / mask / folder / pat.format(r=r), der / f"{stem}.nii.gz"):
                    meta = {**ACQ, "ReconstructionAlgorithm": algo, "GazeBin": gaze,
                            "ETMask": "clean (fixations only, blinks and saccades removed)",
                            "FractionOfCleanReadoutsKept": kept,
                            "Sources": [f"bids::sourcedata/{sub}/rawdata/{sub}_T1wLIBRE.dat",
                                        f"bids::derivatives/mreyetrack-binning/{sub}"]}
                    if rec == "regu":
                        meta.update(NumberOfIterations=20, Delta=1.0)
                    L.write_json(der / f"{stem}.json", meta)
    L.put(s / "recon/C.mat", DST / "derivatives/mreyetrack-recon" / sub / "coilmaps" / f"{sub}_desc-coilsens_coilmap.mat")

    # --- derivatives: binning labels
    b = DST / "derivatives" / "mreyetrack-binning" / sub
    L.tree(s / "eyemasks", b / "eyemasks")
    L.tree(s / "recon/bins", b / "bins")



def top_level(L, subs):
    L.write_json(DST / "dataset_description.json", {
        "Name": "MR-EyeTrack: gaze-resolved 3D T1w MRI of the eye with simultaneous eye tracking",
        "BIDSVersion": BIDS_VERSION, "DatasetType": "raw", "License": LICENSE,
        "Authors": AUTHORS, "Funding": FUNDING, "EthicsApprovals": ETHICS,
        "ReferencesAndLinks": ["https://github.com/MattechLab/MR-EyeTrack"],
    })
    L.write_json(DST / "derivatives/mreyetrack-recon/dataset_description.json", {
        "Name": "MR-EyeTrack gaze-resolved reconstructions", "BIDSVersion": BIDS_VERSION,
        "DatasetType": "derivative", "License": LICENSE,
        "GeneratedBy": [{"Name": "monalisa", "Description": "bmMathilda (gridding) and bmSteva "
                         "(compressed-sensing, 20 iterations, delta 1) on pulseq 3D radial-spiral k-space"},
                        {"Name": "MR-EyeTrack", "CodeURL": "https://github.com/MattechLab/MR-EyeTrack"}],
        "SourceDatasets": [{"URL": "../../"}],
    })
    L.write_json(DST / "derivatives/mreyetrack-binning/dataset_description.json", {
        "Name": "MR-EyeTrack gaze labels per k-space readout", "BIDSVersion": BIDS_VERSION,
        "DatasetType": "derivative", "License": LICENSE,
        "GeneratedBy": [{"Name": "MR-EyeTrack", "Description": "S2 binning: sliding-window gaze labels "
                         "(eMask_th0.75_winLen3_sync = 24 ms window, per-subject ET-MRI timing corrected)"}],
        "SourceDatasets": [{"URL": "../../"}],
    })
    L.put(SRC / "sync_corrections.mat", DST / "derivatives/mreyetrack-binning/sync_corrections.mat")
    L.tree(SRC / "pulseq", DST / "sourcedata" / "pulseq")

    rows = [(f"sub-{n:03d}", *DEMOGRAPHICS[n]) for n in subs]
    if not L.dry:
        with open(DST / "participants.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["participant_id", "age", "sex"])
            w.writerows(rows)
    L.write_json(DST / "participants.json", {
        "age": {"Description": "Age at scan", "Units": "years"},
        "sex": {"Description": "Sex", "Levels": {"M": "male", "F": "female"}}})
    L.write_text(DST / "CHANGES", "1.0.0 2026-10-22\n  - Initial release.\n")
    L.write_text(DST / "README.md", README)


README = """# MR-EyeTrack

3D radial-spiral T1-weighted (pulseq LIBRE) MRI of the eye with simultaneous SR Research EyeLink
eye tracking. Each k-space readout (TR 8 ms) is labelled by gaze direction (up / down / left / right)
and each bin is reconstructed independently, giving four motion-resolved images per subject.

## Layout

| path | content |
| --- | --- |
| `sub-*/anat/*_acq-full_rec-regu_T1w` | LIBRE image from all fixation readouts (no gaze binning), regularised (STEVA) |
| `sub-*/anat/*_acq-full_rec-grid_T1w` | same readouts, gridding (MATHILDA) |
| `sub-*/anat/*_acq-mprage_T1w` | reference MPRAGE from the scanner (Siemens IDEA), converted from DICOM |
| `sourcedata/sub-*/rawdata` | Siemens twix `.dat`: `_T1wLIBRE` (main), `_HC` / `_BC` (head / body coil prescans), `_AdjCoilSens` |
| `sourcedata/sub-*/et` | eye tracking: original EyeLink / PsychoPy files (EDF, csv, log, psydat) and the BIDS conversion (`.tsv.gz` + `.json`) |
| `sourcedata/pulseq` | the `.seq` files defining the k-space trajectory (prescan and main) |
| `derivatives/mreyetrack-recon` | gaze-binned and other reconstructions, coil sensitivity maps |
| `derivatives/mreyetrack-binning` | ET masks and per-readout gaze labels, ET-MRI sync corrections |

## File naming

- `acq-up|down|left|right` gaze bin, `acq-full` all readouts.
- `rec-regu` regularised iterative compressed-sensing reconstruction (STEVA, 20 iterations, delta 1);
  `rec-grid` gridding (MATHILDA, the initial estimate of the iterative recon).
- `desc-keep50|keep25|keep5` random undersampling study: percentage of the clean readouts kept.
- Images are 240 x 240 x 240, about 1 mm isotropic, float32, RAS.

## Readout indexing

Of 82 368 acquired readouts (44 segments x 1872 shots), segment 1 of each shot (SI navigator) and the
first 14 shots (warm-up) are dropped, leaving 79 894 (43 x 1858). Labels use the same ordering.
In numpy: `mask.reshape(1872, 44)[14:, 1:].ravel()`; MATLAB is column-major.

## Authors

Jaime Barranco (1,2,3), Yiwei Jia (1,2,4), Eleonora Fornari (3,5), Oscar Esteban (1), Jean-Baptiste Ledoux (3,5),
Mauro Leidi (1,2), Helene Vitali (6,7), Cédric Campos Carvalho (8), Elena Najdenovska (8),
Laura Elena Raileanu (8), Jessica Bastiaansen (9,10), Micah M. Murray (2,3), Benedetta Franceschiello (1,2).

1. HES-SO University of Applied Sciences and Arts Western Switzerland.
2. The Sense Innovation and Research Centre, Lausanne and Sion, Switzerland.
3. Department of Radiology, Lausanne University Hospital and University of Lausanne, Lausanne, Switzerland.
4. Graduate School for Cellular and Biomedical Sciences, University of Bern, Bern, Switzerland.
5. CIBM Center for Biomedical Imaging, Lausanne, Switzerland.
6. Unit for visually impaired people, Italian Institute of Technology, Genova, Italy.
7. Department of Computer Science, Bioengineering, Robotics and Systems, University of Genoa, Genova, Italy.
8. School of Engineering and Management Vaud (HEIG-VD), HES-SO University of Applied Sciences and Arts Western
   Switzerland, Yverdon-les-Bains, Switzerland.
9. Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital,
   University of Bern, Switzerland.
10. Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland.

Lead contact: Benedetta Franceschiello (benedetta.franceschiello@hevs.ch).
Technical contact: Jaime Barranco (jaime.barrancohernandez@hevs.ch).
Code: https://github.com/MattechLab/MR-EyeTrack

## Ethics, funding and competing interests

The MR-Eye Track protocol was approved by the local ethics committee through the Swiss BASEC
(protocol number 2024-00165, 27 June 2024). The original protocol was approved by the cantonal ethics
committee (protocol number 2018-00240). Funding is listed in `dataset_description.json`.
The authors had a previous patent application related to the original study; it has expired, and the
authors declare no competing interests.

## Not included

The split k-space (mitosius, ~50 GB per subject) and the ROI-PCA compressed k-space are not shared; they
can be rebuilt from `sourcedata` + `derivatives/mreyetrack-binning` with the MR-EyeTrack pipeline.
The reference MPRAGE is not available for sub-010 and sub-011.

## Known issues

- Eye tracking of sub-007 and sub-009 started after the scan began (offsets -451 / -456 ms), unlike the
  other subjects (about +61 ms). All offsets, plus the per-subject EyeLink clock-rate error, are in
  `derivatives/mreyetrack-binning/sync_corrections.mat` and are already applied to the `_sync` labels.
- The `acq-full_rec-grid` images of sub-001 to sub-005 are stored without intensity normalization (values
  around 1e-8, against about 1-10 for every other image). They are intact, but a viewer with a fixed display
  window shows them black; rescale the display window to view them.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["symlink", "hardlink", "copy"], default="symlink")
    ap.add_argument("--subjects", type=int, nargs="*", default=list(range(1, 16)))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    L = Linker(a.mode, a.dry_run)
    for n in a.subjects:
        build_subject(L, f"sub-{n:03d}")
    top_level(L, a.subjects)
    print(f"{'would link' if a.dry_run else 'linked'} {L.n} files into {DST} ({a.mode})")
    if L.missing:
        print(f"{len(L.missing)} expected sources missing:")
        for m in L.missing:
            print("  ", m)


if __name__ == "__main__":
    main()
