# Readout undersampling study — pulseq cohort (15 subjects)

Second iteration of the undersampling sub-study (ISMRM 2026 abstract, 3 IDEA
subjects, `analysis/comparison/`). Question: how many readouts can a gaze bin
lose before the motion-resolved image degrades, for the brain and for each
orbital structure.

Code lives here; data products live in `results/us_pulseq/` (gitignored).

## Pipeline

| step | script | output (`results/us_pulseq/`) |
| --- | --- | --- |
| 1 | `01_mat2bids.py` | `bids/derivatives/reconstructions/` — 510 volumes, 34 per subject |
| 2 | `02_run_mriqc.sh wobin\|regu\|grid` | `mriqc_*/`, `work_*/` (kept: masks come from them) |
| 3 | `03_run_aeye.sh` | `aeye/sub-XXX/*_aeye_both_gated.nii.gz` (MPRAGE space) |
| 4 | `04_register_masks.py` | `masks/sub-XXX/{brain_dseg,eyes_mni,aeye_labels}.nii.gz` (LIBRE space) |
| 5 | `05_region_metrics.py` | `metrics/region_metrics.csv` |
| 5b | `05b_dcte_crop.py` | `metrics/dcte_crop.csv` |
| 6 | `06_mreyeqc_iqms.py iqms\|predict` | `metrics/mreyeqc_iqms_*.csv`, `mreyeqc_ratings_*.csv` |
| 7 | `07_analysis.py` | `analysis/*.csv` |
| 8 | `08_figures.py` | `figs/` |
| QC | `qc_masks.py` | `qc/sub-XXX_aeye_axial.png`, `qc/sub-XXX_brain.png` |

Environments: `usqc` (steps 1, 4, 5, 7, 8 — nibabel, scikit-image, statsmodels,
SimpleITK), `mreyetrack` (A-eye wrapper), `mreyeqc` (python 3.9 +
fetal_brain_qc, step 6, run with `PYTHONPATH=/home/debi/jaime/repos/mreyeqc`).
ANTs 2.6.3 is native; MRIQC and A-eye run in Docker.

## Volumes and naming

Per subject: woBin (all 79 894 readouts) and 4 gaze bins × 4 levels, each as
REGU (`steva`, CS with spatial TV, δ=1, 20 it.) and GRID (`mathilda`).
Names follow last year's `results/bids`:
`sub-XXX_rec-stevaFull_T1w` (woBin), `sub-XXX_acq-<up|down|left|right>_rec-steva<full|50p|75p|95p>_T1w`,
same with `mathilda`. `full` = clean bin; `50p/75p/95p` = 50/75/95 % of the
clean readouts removed at random (`S2_1_eyeMask_undersampling_study.m`).

## Decisions worth knowing

- **All NIfTIs are regenerated from the `.mat` recons.** The NIfTIs already in
  `data/study/sub-XXX/recon/` come from three converters and disagree: the Dec
  2025 clean bins have an identity affine and a different array orientation
  (corr 0.83 with the twix conversion of the same volume), `clean_0.XX` and the
  non-`_twix` woBin carry the MPRAGE affine with 0.97 mm voxels. `01_mat2bids.py`
  applies the confirmed reorientFcn and the twix affine of `woBin/x0.nii.gz` to
  everything and checks it reproduces the existing `mat2nii_twix` outputs
  (relative error < 1e-7).
- **Intensity scale.** `woBin/x0` is unnormalised on some subjects (max ~1e-8),
  so every volume is scaled to p99.5 = 1000. All metrics are scale invariant.
- **One mask set per subject**, built on woBin REGU and reused for all 34
  volumes, so a metric change across levels is a change in the image. MRIQC
  re-segments every volume itself; its IQMs are reported as MRIQC defines them.
- **A-eye labels**: MPRAGE → woBin REGU, rigid (Mattes MI) on the whole head,
  then rigid refinement inside a 15 mm shell around the orbits. One transform
  per subject, applied to every bin (CLAUDE.md: a per-bin registration would
  absorb the eye displacement). The lens label therefore does not follow the
  lens in eccentric gaze bins; lens metrics compare levels *within* a bin.
- **SSIM** is the mean of the 3D SSIM map (Gaussian σ = 1.5) inside each
  region. Last year's notebook ran SSIM on the flattened vector of masked
  voxels, which slides its window along array order; values are not
  comparable to last year's.
- **The MNI "eye" mask is not an eye mask.**
  `tpl-MNI152NLin2009cAsym_res-01_desc-eye_mask.nii.gz` (same md5 in the
  templateflow cache, `~/jaime/data/templateflow` and the A-eye folder) is one
  1.17 L component: the anterior-inferior head block in front of y = −4 mm and
  below z = 1 mm, minus the brain — orbits plus nose, sinuses, face and the air
  in front of it (`results/us_pulseq/qc/mni_eye_mask.png`). Last year's "Eyes"
  region was this block. It is kept here as `eyes_mni` for continuity; the
  eye-specific regions are the A-eye ones.
- **sub-010 and sub-011 have no MPRAGE of their own.** Their
  `dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/*.nii.gz` and `.dcm` are
  byte-identical to sub-001's (same md5, acquisition time 16:22:25), and the
  filer has none either. Consequences: (1) their A-eye labels come from
  multi-atlas fusion of the other 13 subjects' LIBRE labels (`04b_atlas_labels.py`,
  validated leave-one-out against MPRAGE-derived labels); (2) their
  `woBin/x0.nii.gz` affines — and so every NIfTI here — carry sub-001's
  direction cosines. Voxel-wise metrics are unaffected; anything that relies on
  the absolute orientation of those two subjects is not.
- **DCTE: use `dcte_crop`, not the legacy `dcte`.** The legacy DCTE (exactly as in
  `analysis/comparison/dcte/dcte.ipynb`) zeroes everything outside the region, and
  the mask border then dominates the high-frequency energy: the WM mask filled with a
  constant scores 21.6 against 21.2 for the real image. Tissue-level legacy DCTE
  therefore hardly moves with undersampling. `05b_dcte_crop.py` takes the same
  ratio on a box around the eyes (lens + vitreous) and the brain with nothing
  zeroed; REGU eyeball −4/−17/−45 %, brain −10/−23/−53 %. DCTE also rises with
  noise, so the non-binned image (more readouts, less noise) scores lower than the
  clean bin; compare levels within one reconstruction type only.
- **MReyeQC label bug** (mreyeqc @ 8420cfb): `SEGM["LENS"] = 0` collides with
  background, so every lens IQM and `seg_cnr`/`seg_cjv` are computed on
  background + lens. `06_mreyeqc_iqms.py --labels fixed` patches it in-process;
  `--labels as-trained` reproduces it, which is what the shipped model saw.
- **MRIQC crashes and failures.** MRIQC aborts the whole run on the first crashed
  node (`config.nipype.stop_on_first_crash`, no CLI flag); `02_run_mriqc.sh` turns
  that off through a python entrypoint. Its Harmonize step normalises by the median
  of WM voxels with p(WM) >= 0.9 after erosion; LIBRE's weak WM/GM contrast leaves
  ~2k such voxels and erosion emptied the mask on 17/510 volumes (11 of sub-015),
  giving an all-zero image and a crashed IQM node. `patches/anatomical.py` (mounted
  by the run script) keeps the uneroded mask only when erosion empties it; the 17
  stale harmonize nodes were moved to `results/us_pulseq/stale_harmonize/` and
  recomputed. All 510 volumes now have MRIQC IQMs.
- **MReyeQC rating model** (`06 … predict`) needs scikit-learn 1.1.3 (the pickle's
  version): run it with the `mreyeqc_model` env.

## Results (n = 15, REGU unless stated)

| | 50 % | 75 % | 95 % |
| --- | --- | --- | --- |
| optic nerve–fat CNR | −1.6 % (≈) | −2.4 % (≈) | −13.8 % |
| rectus–fat CNR | −1.1 % (≈) | −3.1 % (≈) | −22.0 % |
| lens–vitreous CNR | −9.3 % | −23.3 % | −39.3 % |
| vitreous CV | +10.1 % | +19.5 % | +44.3 % |
| eyeball edge width | +3.6 % | +11.5 % | +26.9 % |
| WM–GM CNR | −3.5 % (≈) | −6.6 % | −16.9 % |
| MRIQC CJV / CNR | −0.2 / −0.1 % | −4.1 / +2.1 % | −1.6 / −3.6 % |
| MRIQC FWHM | +3.9 % | +12.5 % | +45.4 % |

≈ = equivalent within ±5 % (TOST p < 0.05). Readout budget (`analysis/knee.csv`):
nerve–fat CNR plateaus above 861 [776–1,060] readouts per bin, rectus–fat above
1,399 [1,176–2,272] (GRID: 1,305 and 3,104); eyeball edge width (+0.09 mm per
halving) and vitreous CV (+0.019) have no plateau. Clean bins hold
12,619 ± 3,925 readouts (2,126–17,343). Abstract draft: `abstract_ismrm2027.md`;
figures: `results/us_pulseq/figs/`.
