# Handoff — ROI-PCA compressed k-space for model training

Paste this as the opening prompt of a new chat that does not share this
project's memory. If the new chat runs in `/home/debi/jaime/repos/MR-EyeTrack`,
most of this loads automatically from the project memory and you only need the
last two sections.

---

## Context

MR-EyeTrack acquires 3D radial-spiral (LIBRE) k-space of the head with a
**52-channel** coil array, plus synchronised SR Research EyeLink eye tracking.
Recon uses the **monalisa** MATLAB library (`/home/debi/MatTechLab/monalisa`).

A side project compressed the raw k-space **on the channel axis only**, so a
model can be trained on a fraction of the data without losing the orbits. The
method is ROI-PCA: build the coil covariance over an eye ROI, keep the leading
eigenvectors, and linearly recombine 52 physical coils into `nv` virtual
channels. Trajectory, readout length and spoke count are untouched.

Code: `recon/ROI-PCA/` and `recon/ROVir/` (shared implementation, see
`recon/ROI-PCA/README.md`).

## The data

**15 subjects**, `sub-001` … `sub-015`, at
`data/study/sub-NNN/recon/ROI-PCA/kspace_ROI-PCA_8_woBin.mat` (~3.2 GiB each,
49 GB total). MATLAB v7.3 = HDF5, so `h5py` reads it directly.

```text
y     [38 349 120 x 8]   complex single   compressed k-space
t     [3 x 480 x 79894]  single           trajectory (kx, ky, kz), 1/mm
ve    [1 x 38 349 120]   single           volume elements (density weights)
meta                                      nv, nCh, retention, normalisation
```

- `38 349 120 = 480 x 79894` — 480 samples per readout, 79 894 readouts.
- `y` is `[nPt x nv]`; the channel axis is **last**.
- `t` and `ve` are **shared across channels** and identical for every `nv`.
- Reconstructed FoV 240 mm, matrix 240³ (1 mm isotropic).
- `nShot_off = 14` warm-up shots and the first (SI) segment of each shot are
  already removed.

**`y` is normalised** by mean |x0| over the eye ROI (recorded in
`meta.normalised`). If you reconstruct with monalisa's STEVA, note `delta` is
scale-dependent, so any comparison must use the same normalisation.

## Channel counts available

All 15 subjects also have `nv = 4` and `nv = 6` built
(`data/study/sub-NNN/recon/mitosius/ROI-PCA_{4,6,8}/woBin/`); only `nv = 8` is
exported as a single `kspace_*.mat`. Re-export any count with
`recon/ROVir/R5_export_compressed.m` (`variant='ROI-PCA'`, set `nv`).

Measured against each subject's own 52-channel **gridded** reconstruction
(no regularisation, so nothing hides the loss):

| nv | compression | ROI corr (min) | NRMSE (median) | SNR change (median / worst) |
|---:|---:|---:|---:|---:|
| 4 | 13.0x | 0.9992 | 0.0221 | -28.7% / -48.6% |
| 6 | 8.7x | 0.9994 | 0.0158 | -21.8% / -34.6% |
| **8** | **6.5x** | **0.9998** | **0.0083** | **-12.8% / -23.0%** |

**nv = 8 is the recommended operating point.** Fidelity never discriminates
(even 4 channels stays above 0.999 on every subject); SNR is the only axis that
does. Per-subject absolute SNR is in `recon/ROI-PCA/snr_by_channels.csv`.

## Things that will bite you

- **sub-005 has the lowest SNR at every channel count** (4.32 even at 52
  channels, vs 6.08 cohort median). It is the hardest subject, not a
  representative one — do not use it alone to sanity-check a model.
- **sub-015 is the most favourable subject** on both ROI energy and SNR cost.
  Early single-subject numbers from it (12% SNR loss at nv=4) became 29% median
  once all 15 were measured. Do not generalise from one subject here.
- **Only the unbinned (`woBin`) case is compressed.** That is fully sampled.
  The 4 gaze-direction bins are undersampled, where coil diversity matters more
  and these channel counts may not hold.
- **This is not anonymised.** The whole head is reconstructable from the
  compressed data. ROI-PCA keeps the eyes; it does not remove anything else.
  (The sibling `recon/ROVir/` variant does suppress the head, at much worse
  compression — see its README.)

## If you need to reconstruct from it

`recon/ROVir/R4_rovir_recon.m` (local) or
`recon/ROVir/hpc/S4_rovir_recon_chacha.m` (SLURM) do gridding (`bmMathilda`)
then compressed-sensing (`bmSteva`). They read the **mitosius** form, not the
packaged `.mat`:

```text
data/study/sub-NNN/recon/mitosius/ROI-PCA_8/woBin/     <- pass this to bmMitosius_load
data/study/sub-NNN/recon/ROI-PCA/C_rovir_8.mat         <- matching virtual coil maps
```

Same contents as `kspace_*.mat`, different packaging. Reconstructions already
exist for all 15 subjects at all three counts (`x0_*` gridded and `x_steva_*`),
with NIfTI versions for visual inspection.
