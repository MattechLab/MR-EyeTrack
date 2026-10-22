# ISMRM 2027 abstract draft — readout undersampling, pulseq cohort

**Title:** How many readouts does a gaze bin need? Structure-specific undersampling limits for motion-resolved eye MRI

**Authors (confirm):** Jaime Barranco (1,2), Yiwei Jia (1,2,3), Jean-Baptiste Ledoux (4,5), Eleonora Fornari (4,5),
Cédric Campos Carvalho (6), Elena Najdenovska (6), Laura Elena Railenau (6), Benedetta Franceschiello (1,2)
— affiliations as in the ROI-PCA abstract.

**Primary:** Acquisition & Reconstruction — Subsampling: Parallel imaging, compressed sensing
**Secondary:** Neuro — Head & Neck/ENT (alternative: Analysis Methods — Data Processing)
**General keyword:** Abbreviated/Rapid MRI (or own: "Image quality metrics")

## Impact (40 words max) — 40 words

Readout budgets per gaze position make motion-resolved eye MRI plannable: with our 3D radial T1-weighted LIBRE protocol (TR 8 ms), about 1,400 readouts (11 s of fixation) preserve orbital contrast with compressed sensing; eyeball focus needs about 7,000 (one minute).

## Synopsis (100 words max) — 100 words

Motivation: In motion-resolved eye MRI, each gaze bin keeps readouts from stable fixation, so image quality depends on fixation duration.

Goals: How many readouts per gaze bin preserve brain and eye-structure image quality?

Approach: Fifteen participants underwent eye-tracked 3D radial fat-suppressed (LIBRE) T1-weighted MRI (TR 8 ms); readouts within four gaze bins were retrospectively undersampled; images were reconstructed with gridding and compressed sensing and assessed with MRIQC and eye-structure metrics.

Results: With compressed sensing, optic nerve and rectus contrast stayed stable down to 900–1,400 readouts per gaze bin (7–11 s; gridding: 1,300–3,100), whereas eyeball focus needed about 6,900 (55 s).

## Body (750 words max) — 749 words

### Introduction

Motion-resolved MR-Eye reconstructs one image per gaze position from readouts acquired while the eye holds it, using synchronised eye tracking1–3, so the readouts per gaze bin depend on how long the participant fixates. In three subjects, quality was preserved up to 50% readout undersampling4, but a percentage does not transfer between bins of different size, and whole-brain masks cannot show which orbital structure fails first. We ask how many readouts a gaze bin needs, per eye structure, in 15 subjects.

### Methods

Fifteen healthy participants were scanned at 3T (MAGNETOM Prisma Fit, Siemens; 64-channel head/neck coil) with a Pulseq5 T1-weighted GRE with LIBRE water excitation6 and a 3D radial spiral-phyllotaxis trajectory7 (82,368 readouts, TR 8 ms, 11 min, 1 mm isotropic) with synchronised eye tracking (EyeLink 1000 Plus) during a four-position fixation task. Fixation readouts formed four "clean" gaze bins, each randomly undersampled by 50, 75 and 95% (240 binned datasets), plus one non-binned dataset per subject. All were reconstructed (240³) by gridding (GRID) and by compressed sensing with anisotropic spatial total variation (REGU; δ=1, 20 iterations) in Monalisa8.

Masks from each subject's non-binned REGU image served every reconstruction: MRIQC9 brain tissues, and orbital structures (lens, vitreous, optic nerve, fat, recti) segmented by A-eye10 on the same-session MPRAGE and rigidly registered to LIBRE (Fig1A). Two subjects without a usable MPRAGE received multi-atlas labels (leave-one-out Dice 0.94 eyeball, 0.74 optic nerve).

We computed MRIQC IQMs; 3D SSIM to the clean bin per structure; tissue CNR (optic nerve–fat, rectus–fat, lens–vitreous, WM–GM); vitreous coefficient of variation (CV); eyeball edge width (boundary-profile slope); DCT energy ratio (DCTE)11, a focus measure, on a box around each structure; and MReyeQC eye IQMs and rating model (FetMRQC12 adapted to the orbit13; trained on non-fat-suppressed T1-weighted images, so not directly applicable to LIBRE). Gaze bins were averaged per subject. Each level was compared with the clean bin by Wilcoxon tests (Holm-corrected), and a change was called negligible when an equivalence test placed it within ±5%. To express limits as a readout count N, each metric was fitted over all binned datasets with two line segments in log2N, flat above a threshold N* and falling below it; N* (95% CI: subject bootstrap) is reported only where this beat one straight line (F-test p<0.05).

### Results

Clean bins held 12,619 ± 3,925 readouts (101 ± 31 s of fixation; Fig1B).

With REGU, optic nerve–fat and rectus–fat CNR stayed within ±5% of the clean bin up to 75% undersampling (p≤0.025) and fell at 95% (−14%, −22%). Lens–vitreous CNR (−9%, −23%, −39%) and vitreous CV (+10%, +20%, +44%) degraded from 50%, and eyeball edges widened by 4%, 12% and 27%. Eyeball and brain DCTE fell by 4/17/45% and 10/23/53% (Fig5). At 50%, SSIM to the clean bin matched or exceeded that of the non-binned image in every structure (eyeball 0.73 vs 0.65); at 75% all fell below (Fig2).

In readouts (Fig3), REGU orbital contrast stayed stable down to 861 [776–1,060] readouts per bin for the optic nerve and 1,399 [1,176–2,272] for the recti (7 and 11 s of fixation); GRID needed 1,305 and 3,104. Eyeball edge width (+0.09 mm) and vitreous CV (+0.02) worsened with every halving of readouts without stabilising, and eyeball DCTE was stable only down to 6,889 [5,997–9,410] readouts (55 s). At equal readout count, edge width and CV were not significantly affected by whether a bin was large and heavily undersampled or small and lightly undersampled (p≥0.26). 

MRIQC's CJV and CNR missed undersampling in REGU (within ±5%) and in GRID until 95%, whereas FWHM tracked it (REGU +4/+13/+45%; Fig4). Applied anyway, the MReyeQC rating decreased with undersampling (REGU: 1.27 clean, 0.90 at 95%).

### Discussion

Readout count, not undersampling percentage, set image quality: at 75%, the 5 of 60 bins from shorter fixation fell below the rectus threshold; others stayed far above. Under total-variation regularisation, contrast measures (fixed-mask tissue CNR, MRIQC's CJV/CNR) missed the blurring that edge width and DCTE detected, so quality checks need sharpness, not only contrast. Seconds are specific to this protocol (TR 8 ms); other contrasts, resolutions and trajectories change SNR and k-space coverage and need their own readout calibration. Limitations: random removal thins each bin uniformly, whereas shorter or interrupted fixation (blinks, drowsiness) removes contiguous runs of spokes, leaving subject-specific angular coverage and point-spread functions; and the lens mask, defined on the MPRAGE at near-central gaze, is only approximate in eccentric gaze positions.

### Conclusion

About 1,400 readouts per gaze bin (11 s at TR 8 ms) preserve orbital contrast with compressed sensing; eyeball focus needs about 7,000 (55 s).

## Figure captions (500 characters max each)

Fig1. Design. A: one subject's 'left' gaze bin without binning, from the clean bin and after removing 50, 75 and 95% of its readouts, reconstructed with compressed sensing (REGU) and gridding (GRID); orbit and mid-sagittal slices, readouts per image on top. A-eye labels (lens, vitreous, optic nerve, recti) on the non-binned orbit. B: readouts per gaze bin for 15 subjects x 4 bins at each level (bar: median); right labels give fixation time at 8 ms per readout. (464 chars)

Fig2. SSIM to the clean gaze bin per structure after removing 50, 75 and 95% of its readouts (mean ± SD over 15 subjects of the 4-bin mean), REGU solid, GRID dashed. Grey band: SSIM of the non-binned image to the clean bin, i.e. the cost of ignoring eye motion. At 50% the undersampled bin stays at least as close to the clean bin as the non-binned image; at 75% it falls below for every structure. REGU preserves brain, nerve and muscle structure better than GRID. (465 chars)

Fig3. Readout budget. Each point is one subject x gaze bin x level (REGU, colour = level); lines are fits over all points with subject offsets (REGU black, GRID dashed). Optic nerve-fat and rectus-fat CNR and eyeball DCTE stay stable down to a threshold (blue line, 95% CI band): 861, 1,399 and 6,889 readouts for REGU (GRID: 1,305, 3,104, no DCTE threshold). Eyeball edge width and vitreous CV never stabilise. x axis: readouts per bin and fixation time. (455 chars)

Fig4. Top: MRIQC IQMs (whole image, MRIQC's own brain segmentation per image). Bottom: MReyeQC eye IQMs on fixed A-eye labels, and the pretrained MReyeQC rating. NB no binning, C clean bin, 50/75/95% undersampling; REGU solid, GRID dashed; mean ± SD over 15 subjects; arrows give the better direction. With REGU, CJV and CNR miss undersampling while FWHM, Dietrich SNR and EFC track it; globe and lens SNR fall from 50%, rectus SNR does not. (441 chars)

Fig5. REGU summary: value in the clean bin and % change from it after removing 50, 75 and 95% of its readouts (subject means). ≈: change shown to lie within ±5% (equivalence test); *: change shown to be non-zero (Wilcoxon, Holm); ≈*: small but consistent across subjects. Orbital tissue contrast is preserved to 75%; lens and vitreous metrics degrade from 50%, eyeball edges widen and DCTE drops progressively; MRIQC CJV and CNR miss undersampling. (448 chars)

Preview figure: `results/us_pulseq/figs/preview.png` (no caption).

## Acknowledgements

This work was supported by the Swiss National Science Foundation (SNSF) project 2.0 MR-Eye (grant #220433) and MR-Eye
Track, funded by the HES-SO University of Applied Sciences and Arts Western Switzerland (grant 129816/IA-RECHERCHE23-19).
The study protocol received positive approval from the local Ethics Committee on June 27, 2024 (BASEC 2024-00165).
Code: https://github.com/MattechLab/MR-EyeTrack (analysis/us_study_pulseq), MReyeQC https://github.com/jaimebarran/mreyeqc.

## References

1. Franceschiello B, et al. 3-Dimensional magnetic resonance imaging of the freely moving human eye. Prog Neurobiol. 2020;194:101885.
2. Jia Y, et al. MR-Eye: a high-resolution motion-resolved MRI protocol for anatomical imaging of the human eye. Proc Intl Soc Mag Reson Med. 2025;33:0639.
3. Barranco J, et al. Protocol for dynamic eye-brain imaging by integrating an eye-tracking system with magnetic resonance imaging. STAR Protoc. In press.
4. Barranco J, et al. Readout undersampling in motion-resolved eye MRI: quantitative effects on image quality. Proc Intl Soc Mag Reson Med. 2026 (abstract 04355).
5. Layton KJ, et al. Pulseq: a rapid and hardware-independent pulse sequence prototyping framework. Magn Reson Med. 2017;77(4):1544-1552.
6. Bastiaansen JAM, et al. Flexible water excitation for fat-free MRI at 3T using lipid insensitive binomial off-resonant RF excitation (LIBRE) pulses. Magn Reson Med. 2018;79(6):3007-3017.
7. Di Sopra L, et al. An automated approach to fully self-gated free-running cardiac and respiratory motion-resolved 5D whole-heart MRI. Magn Reson Med. 2019;82(6):2118-2132.
8. Leidi M, et al. Monalisa: an open source, documented MATLAB toolbox for magnetic resonance imaging reconstruction. bioRxiv. 2025. doi:10.1101/2025.06.06.658299
9. Esteban O, et al. MRIQC: advancing the automatic prediction of image quality in MRI from unseen sites. PLoS One. 2017;12(9):e0184661.
10. Barranco Hernandez J, et al. A-eye: automated 3D segmentation of healthy human eye and orbit structures and axial length extraction. 2025. doi:10.26039/TA7F-X088
11. Romanin L, et al. Exploring the limits of scan time reduction for ferumoxytol-enhanced whole-heart angiography in congenital heart disease patients. J Cardiovasc Magn Reson. 2025;27(1):101854. doi:10.1016/j.jocmr.2025.101854
12. Sanchez T, et al. FetMRQC: automated quality control for fetal brain MRI. PIPPI MICCAI Workshop. 2023. doi:10.1007/978-3-031-45544-5_1
13. Telley C. Integration of a quality control system for eye MRI images [master's thesis]. Lausanne: HES-SO Master of Science in Engineering; 2025.
