# Decoding eye state from raw k-space — what was tried and what it showed

**The eye moves and the binning works.** The optic nerve profile — a monotone
10–20× falloff from the optic nerve head to the orbital apex, with a null
control at 0.00 mm — is unambiguous. **The motion just does not reach individual
readouts at a level any model can use.** That is a property of the acquisition,
not of the labels, the compression, the architecture or the window length.

This directory holds the work behind that conclusion. `DESIGN.md` is the
chronological log, including the wrong turns and their corrections; this file is
the summary you can act on.

---

## 1. The conclusion, and the evidence for each half

### The eye really does rotate between gaze bins

Displacement measured in bands along each optic nerve (`nerve_profile.py`),
using A-eye labels carried into LIBRE space by a single rigid+affine transform:

| distance from globe | left L–R | left U–D | right L–R | right U–D |
|---|---:|---:|---:|---:|
| 10–14 mm (ONH) | +1.03 | +0.60 | +1.10 | +0.64 |
| 14–20 mm | +0.37 | +0.41 | +1.02 | +0.45 |
| 20–29 mm | +0.20 | +0.10 | +0.11 | +0.13 |
| 29–38 mm (apex) | +0.08 | +0.05 | +0.10 | +0.04 |

Monotone in all four series. This is the tethered-nerve signature and no
nuisance reproduces it: bulk head motion shifts every band equally, and the
brain control reads **0.00 mm**. The two nerves are independent structures and
agree. L–R exceeds U–D on both sides, matching the stimulus ratio.

Supporting: the eye tracker says the eye moved **18.0°** (L–R) and **7.90°**
(U–D); within-bin gaze scatter is 9–15 px IQR, so the bins are pure; and the
ET→readout alignment is correct.

### But single readouts carry nothing usable

Two independent measurements, neither involving a model.

**Gaze** (`pairtest_raw.py`). Take readout pairs pointing the same way — the
trajectory revisits directions ~80 s apart, hence usually in different gaze
states — and compare `||y_a − y_b||²` for same-bin against cross-bin pairs:

| sub | ROI-PCA nv=8 | uncompressed 52 ch |
|---|---|---|
| 015 | 0.9903 [0.9824, 0.9973] | 0.9861 [0.9803, 0.9932] |
| 013 | 0.9421 [0.9336, 0.9501] | 0.9501 [0.9440, 0.9567] |
| 010 | 0.9800 [0.9676, 0.9945] | 0.9608 [0.9499, 0.9731] |

Ratio 1.000 would mean cross-bin pairs differ by exactly the gaze signal. Every
CI sits **entirely below 1**, uncompressed as well as compressed.

**Motion** (`pairtest_motion.py`). Same method, still×motion pairs against
still×still. Point estimate ~0; one marginal hit (sub-001, ROI-PCA) failed to
replicate on that subject's own uncompressed data. Bound ≈ d′ < 0.12 per readout.

### Why both are true at once

An image is reconstructed from ~15 000 readouts, so it sees ~120× more evidence
than one readout does. A 1 mm displacement of a structure with lens/vitreous
CNR ≈ 2.7 shows up in the image and vanishes in the single measurement.

---

## 2. What was tested, and ruled out

| task | input | result | chance |
|---|---|---|---|
| 4-way gaze | raw k-space, ROI-PCA nv=8 | 0.2833 | 0.250 |
| 4-way gaze | raw k-space, SVD nv=8 | 0.2728 | 0.250 |
| 4-way gaze | raw k-space, uncompressed 52 ch | 0.2557 | 0.250 |
| 4-way gaze | per-readout statistics, 8 / 52 ch | 0.2518 | 0.250 |
| binary motion | per-readout statistics, 8 ch / 52 ch | 0.5018 / 0.5000 | 0.500 |
| binary motion | raw k-space + TCN, rebuilt labels | **0.5147 vs null 0.5206** | 0.500 |
| 4-way gaze | **unsupervised** k-means on window embeddings | z +0.52, p = 0.22 | — |
| 4-way gaze | **unsupervised** autoencoder + frozen probe | z +0.68, p = 0.035 | — |

The weak 0.27–0.29 gaze effect is **insensitive to everything**: spoke direction
(isotropic), `|k|` band (uniform), window length 0.5–5 s (flat), coil subspace,
and channel count 8/16/52. A real eye-position signal must care about at least
the first two.

The final motion run carries its own circular-shift null, and **the null scored
higher than the true labels** — so the ~0.52 both reach is what the architecture
produces from any temporally blocky label.

### Unsupervised learning finds the same weak effect, not a hidden one

Every model above is supervised. `unsupervised.py` asks the other question: does
the k-space hold gaze structure a classifier failed to *extract*? Labels are used
only to evaluate what was learned, never to fit it.

Both methods run on the per-direction-cluster detrended residual. That is not
optional: the dominant variance in raw readouts is the spoke direction (a probe
predicts it at 0.9865), so anything unsupervised on raw data rediscovers the
trajectory and nothing else.

| method | mean z | median z | subjects z>0 | p |
|---|---:|---:|---:|---:|
| k-means (k=4) on window embeddings | +0.64 | +0.52 | 6/11 | 0.224 |
| autoencoder (32-d code) + frozen linear probe | +1.03 | +0.68 | 8/11 | 0.035 |

The autoencoder trains properly (reconstruction MSE 0.104 -> 0.058), so its code
does capture the dominant structure. A linear probe on that frozen code reaches
**z ~ +1.0 -- statistically indistinguishable from the supervised linear probe**
(z ~ +1.0, p = 0.015-0.039). The same weak effect, reached from the other
direction; nothing new was uncovered.

k-means finds nothing (p = 0.22). Its raw agreement of ~0.31 looks above the
0.25 chance line, but its shift null sits at ~0.30 as well -- that offset is the
Hungarian cluster-to-class matching, not signal.

Two caveats. The probe evaluates on non-overlapping windows, so only ~94 test
windows per subject, which is why per-subject z is noisy (-0.64 to +3.77). And
p = 0.035 is uncorrected across the many tests in this document; it is at the
same marginal level as everything else here, not above it.

### Compression: the choice does not matter

| | eye-ROI energy | head energy | linear probe | trained |
|---|---:|---:|---:|---:|
| ROI-PCA nv=8 | 99.3% | 65.4% | 0.2895 | 0.2833 |
| SVD nv=8 | 16.6% | 15.2% | 0.2873 | 0.2728 |
| SVD nv=16 | 34.8% | 28.4% | 0.2873 | — |

A **6× difference in retained eye energy changes nothing downstream**. Combined
with the uncompressed pair-test bound, no ROI mask — including A-eye — can help,
because compression is a linear map and the signal is absent before it.
`recon/ROI-PCA/README.md` reached the same conclusion from the energy side.

**ROI-PCA nv=8 remains the right export.**

---

## 3. The label rebuild — worth having regardless

The original labels used a sliding window of `winLen=10, th=0.75`: a readout is
kept only if 75% of an 80 ms window agrees. That is a purity filter designed for
reconstruction binning, and it is the wrong rule for an ML label — a 40 ms
saccade tops out at 0.5 overlap and **can never be labelled at any position**.

`recon/2-Binning/S2_batch_all_winLen.m` rebuilds every mask type for every
subject in one pass (one raw read per subject rather than one per mask type),
with no prompts.

| class | expected RO | winLen10 | winLen3 | winLen3+sync |
|---|---:|---:|---:|---:|
| fixation | 1 006 682 | 94.2% | 95.2% | 95.2% |
| up / down / left / right | ~200 000 each | 93–94% | 95% | 95% |
| **saccade** | 54 974 | **37.9%** | **69.8%** | **70.0%** |
| blink | 149 895 | 90.2% | 95.5% | 95.5% |
| tracking-loss | 25 514 | 81.2% | 90.2% | 90.5% |

Per-subject saccade retention went from a **13× spread (4.5%–58%)** to under 2×
(44%–82%). The old saccade class was a duration-biased sample, and the bias
differed per subject.

The `_sync` correction changes retention by ~0.1% — it moves *which* readouts
carry a label, not how many — and applies

    et_index(k) = ratio × (offset_ms + TimeStamp_ms(k))

with per-subject values from `sync_table.py`. Two things it fixes:

- **a systematic ~61 ms offset every previous bin carried.** The first readout
  lands ~55 ms after the trigger that released it (`PMUTimeStamp` at readout 0);
  a trigger-only estimate misses this entirely. 61 ms is ~7.6 readouts.
- **per-subject EyeLink drift spanning −92 to +118 ppm**, up to 75 ms across a
  scan. It is not a constant and must be measured per subject.

Timing budget: start offset ~6 ms, scanner-vs-stimulus-PC drift +6 ppm (4 ms
per scan), EyeLink drift as above.

Three label generations coexist; nothing was overwritten:

    bins/<type>/eMask_th0.75_winLen10.mat          original
    bins/<type>/eMask_th0.75_winLen3.mat           window fix only
    bins/<type>/eMask_th0.75_winLen3_sync.mat      window + timing   <- use this

**Two subjects need care:** sub-007 and sub-009 have genuine −451/−456 ms ET
start offsets (their eye tracker began after the scan). sub-014's apparent
−2955 ms was a **stray trigger**, resolved to +63.5 ms — its data is fine.

---

## 4. Layout

| script | does |
|---|---|
| `dataset.py` | readout↔label mapping, verified and asserted |
| `precompute.py`, `extract_band_raw.py`, `svd_compress.py` | build k-space bands |
| `decode_gaze.py`, `compare_compression.py` | gaze decode, compression head-to-head |
| `train_small.py` | FiLM/TCN, `--task gaze\|motion`, `--variant`, `--shift` null |
| `pairtest_raw.py`, `pairtest_motion.py` | the model-free bounds |
| `nerve_profile.py`, `nerve_shift.py` | image-domain displacement |
| `sync_table.py`, `sync_audit.py` | timing corrections and cohort audit |
| `label_audit.py` | ET samples vs kept readouts, all three label generations |
| `features/` | per-readout statistics route (chance at every channel count) |

Derived data (~180 GB) under `data/derived/`: `kband` (ROI-PCA bands +
detrend caches), `kband_raw52` (uncompressed), `svd`, `features`,
`features_raw`.

Environment: conda env `mreye-ml` (torch 2.13 + cu130, finufft, nnUNet deps).
The vendored twixtools needs a shim — `scipy.integrate.cumtrapz` was removed in
modern scipy; alias it to `cumulative_trapezoid` before importing.

---

## 5. Corrections to earlier claims in this repo

**The CBMS 2025 result is not a target.** Its ~98% came from separating
redundant periodic MR structure (the SI projection), not motion; with correct
ET-derived classes it is ~50%. Our nulls corroborate it rather than fall short
of it. Related: the positive control in `features/evaluate_raw.py` predicts
*spoke angle* at 0.9865 from the same features — the same class of periodic
structure, and the coincidence of ~98% in both places went unnoticed at the
time. Our pipeline was never exposed to it only because segment 0 of each shot
is dropped first.

**Withdrawn:** the hypothesis that ROI-PCA destroyed phase information the
feature method relied on. It was invented to explain a gap that does not exist.

**Two traps worth not repeating.** Do not filter the trigger train for regular
2.5 s intervals — the *first* interval is systematically short (~1.99 s), so
that discards the genuine first trigger and reports a spurious ~2 s offset on
every subject. And trigger-train alignment cannot fix an offset: the train is
periodic, so any whole-period shift fits equally well. `PMUTimeStamp` at the
first readout is what breaks the ambiguity.

---

## 6. Where the lever actually is

Not in modelling. Six representations produced the same number and two
model-free bounds explain why.

The limit is that a ~1 mm displacement of a structure with lens/vitreous CNR
≈ 2.7 does not imprint on a single 8 ms readout. Higher lens/vitreous contrast,
shorter TR, or an orbit-targeted readout would each change that arithmetic.
Nothing downstream of the acquisition can.
