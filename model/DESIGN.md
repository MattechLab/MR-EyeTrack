# Eye state from raw k-space — data choice and model design

Decisions and measurements behind training a classifier on the compressed
k-space plus the eye tracker. Everything marked *measured* was computed from
the data on disk; the scripts are in the session scratchpad and the reusable
part is `model/dataset.py`.

## 1. Use ROI-PCA nv=8, not ROVir

The intuition for ROVir is sound — penalise the head, keep the orbits, and the
classifier has less nuisance to fight. It does not survive contact with the
data.

**Same probe, same subject, same labels** — a direction-conditioned linear
matched filter, cross-validated over time blocks, per-readout `d'`:

| contrast | ROI-PCA nv=8 (6.5x) | ROVir nv=20 (2.6x) |
|---|---:|---:|
| gaze up vs down | 0.080 | 0.082 |
| gaze left vs right | 0.102 | 0.103 |
| fixation vs blink | 0.042 | 0.032 |

ROVir needs 2.5x the channels to deliver the same information, and is worse on
the blink contrast. Whatever its recon-side merits, it is not a better input.

**ROVir's suppression is mostly already inside the ROI-PCA export.** Both are
52 -> nv linear maps, so any weight vector lying in `span(V_pca8)` is
recoverable from the export by a fixed 8-vector — something a first linear
layer learns for free. Solving the ROVir generalised eigenproblem *restricted*
to that span (sub-015):

| built from | best achievable ROI/head energy ratio |
|---|---:|
| an average physical coil | 0.061 |
| ROI-PCA nv=4 export | 1.49 |
| ROI-PCA nv=6 export | 2.04 |
| **ROI-PCA nv=8 export** | **2.44** |
| all 52 physical coils (true ROVir) | 5.61 |

So the export already permits ~40x an ordinary coil's head suppression, 44% of
what the full array could reach. The remaining gap costs 6.5x compression and a
3x SNR loss in the orbits (ROVir keeps 10.5% of ROI energy at nv=8 vs 99.3%).

**Practical asymmetry.** ROI-PCA is built for all 15 subjects at nv = 4, 6, 8.
ROVir exists for sub-015 only; matching it means re-running R1–R3 on 14
subjects (~40 min each against the 16 GB raw) for ~50 GB of extra data.

**On the visual impression of less background noise.** The two ITK-SNAP panes
were windowed differently (max 0.888 for ROI-PCA vs 1.498 for ROVir), so they
are not comparable by eye. The measured background sigma in the eye ROI goes
the other way: 6.26e-02 for ROI-PCA nv=8 against 1.26e-01 for ROVir nv=20
(`recon/ROI-PCA/README.md`). What ROVir genuinely does is suppress the *head*,
which reads as a cleaner field around the orbits.

**Keep ROVir for one thing:** it is the only variant that defaces. If the
trained model or the data ever leaves the group, ROVir is the anonymisable
path — ROI-PCA reconstructs the whole head.

Stay at **nv = 8**. The task is SNR-limited (below), and nv=4 costs both raw SNR
(-28.7% median vs -12.8%) and reachable suppression (2.44 -> 1.49).

## 2. What the k-space looks like to a model

- LIBRE spokes are **full-diameter projections**: `|k|` runs 1 -> 0 -> 1 over
  480 samples, so **sample 240 is k = 0** on every readout.
- The sample index is therefore *both* a k-space coordinate and an echo time.
  Antipodal spokes cover the same line but at mirrored echo times, so they are
  **not** interchangeable — folding them onto one hemisphere injects a
  systematic artefact (this bit me; it showed up as negative split-half
  reliability).
- Eye-scale structure sits at `|k|` 0.08–0.33 normalised, i.e. samples
  240 ± 20 to 240 ± 80. A 24 mm globe is ±20, a 6 mm lens ±80.
- Spoke direction changes every readout and is the dominant source of
  readout-to-readout variance. Angular neighbours exist (~12 spokes within 1°)
  but land a median 313 s away in time.

**The k=0 self-navigator alone is at chance** (*measured*: all AUCs 0.49–0.51
for blink, saccade and non-fixation on sub-015). DC integrates the whole FoV,
and ROI-PCA retains 65% of head energy against a head that carries 16.5x the
orbits' total energy, so eye motion moves DC by far less than the drift does.
A naive navigator will not work here; the model has to use `|k| > 0`, and
therefore has to be conditioned on the spoke direction.

## 3. Single readouts are hopeless; windows are mandatory

*Measured*: per-readout `d'` is at most ~0.1 for a sustained gaze contrast and
~0.04 for blink. Evidence pools as `d' * sqrt(W)`:

| window | duration | expected AUC at d' = 0.1 |
|---:|---:|---:|
| 16 | 0.13 s | 0.61 |
| 64 | 0.5 s | 0.72 |
| 256 | 2.0 s | 0.88 |
| 1024 | 8.2 s | 0.99 |

This is the physics, not a modelling shortfall: the binned reconstructions need
~19 000 readouts to render the eye, so one readout carries ~1/138 of that
evidence. Predict at readout resolution if you like, but the receptive field
has to span **hundreds of readouts (1–4 s)**.

## 4. Start with gaze state, not fixation vs non-fixation

The requested binary task is the *hardest* place to start, for three reasons.

1. **Event duration fights the integration window.** A saccade is 30–80 ms
   (4–10 readouts) and a blink 100–400 ms (12–50). Both are shorter than the
   window the SNR demands. A sustained gaze state lasts 5 s (≈625 readouts).
2. **The negative class is a grab-bag.** `fixation-ok` complement = saccades +
   blinks + tracking loss + a large "gap" class the 0.75 sliding window would
   not commit to (4.7% of readouts on sub-015, 13.5% on sub-005), which sits at
   event boundaries and is physiologically almost fixation.
3. *Measured*: the direction-conditioned template is reproducible across the
   first and second half of the scan for gaze contrasts (r ≈ +0.13, shuffled
   control ≈ 0) and not for fixation vs non-fixation.

**Recommended framing.** Train the model to output a **gaze state** (4-class,
or better a continuous 2D gaze estimate) at ~0.5 s resolution, then derive
fixation/saccade from its **temporal derivative**: a saccade is a change point
between two decodable states, which is a far stronger signal than a saccade's
own 60 ms of k-space. Blinks are the one non-fixation event with a large direct
signature (the lid sweeps the cornea) and can stay a direct output head.

## 5. Audit the labels before training

*Measured*, all 15 subjects, per-readout fractions:

| | min | median | max |
|---|---:|---:|---:|
| fixation-ok | 0.500 | 0.853 | 0.953 |
| blink-event | 0.006 | 0.074 | **0.415** |
| saccade-only | 0.001 | 0.005 | 0.078 |

A 66x spread in blink rate is not physiology. sub-002 (41.5% blink, 0.11%
saccade), sub-003, sub-005 and sub-007 (21–26% blink) almost certainly have
pupil-detection failure labelled as blink. Two consequences:

- **Subject identity is a shortcut.** With base rates from 0.50 to 0.95, a
  model that only learns which subject it is scores well above chance. Splits
  must be **leave-subjects-out**, and inputs normalised per subject.
- Decide explicitly whether the low-quality subjects are training data, a
  held-out robustness set, or excluded. sub-005 is both the lowest-SNR subject
  and among the worst-labelled — do not sanity-check on it, and do not
  sanity-check on sub-015 either (best on both axes).

## 6. Proposed architecture

Input, per readout `r` in a window of `W = 256` (2.0 s), stride 32:

| stream | shape | why |
|---|---|---|
| `y[r, 160:320, :]` re/im | 160 x 8 x 2 | eye-scale k-space, `\|k\|` <= 0.33 |
| residual vs angular neighbours | 160 x 8 x 2 | removes the static head; precompute once |
| spoke direction `d[r]` | 3 -> 32 | Fourier/SH features, the conditioning input |
| `|k|` axis | 160 | tells the conv where it is in k-space |

```
per-readout encoder   Conv1d over the sample axis (channels = 8 x 2 x 2 streams)
                      FiLM-modulated by the direction embedding      -> e_r in R^128
sequence model        dilated TCN or 4-layer transformer over W readouts
heads                 (a) gaze: 2D regression or 4-way, per readout
                      (b) blink: per-readout logit
                      (c) fixation: derived from d(gaze)/dt, or a third head
```

Non-negotiables, each traceable to a measurement above:

- **Condition on direction.** Without it the dominant input variance is
  trajectory rotation, not eye state.
- **Remove the global phase drift** per subject before training (B0 drift moves
  every readout's phase far more than the eye does).
- **Normalise per subject**, and split leave-subjects-out.
- Keep the **sample axis intact** — do not fold antipodal spokes, and do not
  reduce to magnitude only; the contrast is partly in phase.

## 7. Order of work

1. Cohort label audit — the blink/tracking-loss confusion above, against the
   raw EDF, on the four suspect subjects.
2. Precompute the training tensor: band 160–320, all 15 subjects
   (~800 MB each, ~12 GB total, fits in the 125 GB of RAM).
3. Positive control before any model: decode **gaze direction** from 2 s
   windows with a linear probe, leave-subjects-out. If that does not clear
   AUC 0.7, no deeper model will help, and the problem is upstream.
4. Then the encoder above, gaze head first, fixation/saccade derived.

`torch` is not installed in any conda env on this workstation; the GPU is a
10 GB RTX 3080, which is ample for this input size.

---

# Results of the positive control (section 7 step 3)

**The control did not come out positive.** Recording this properly because the
sections above were written before it ran and read more optimistically than the
evidence now supports.

## What was built

- `model/precompute.py` -> `data/derived/kband/`, 11 subjects (2, 3, 5, 7
  dropped for blink-rate), band 160:320 (`|k|` <= 0.33), 780 MiB each, 8.4 GB.
- `model/probe_gaze.py`, `probe_gaze2.py` — windowed 4-way gaze classifiers.
- `model/dprime.py`, `gaze_test.py` — bias-corrected separation with a
  circular-shift null.
- `model/effect_size.py` — image-domain effect size from the binned recons.

## What happened

| probe | result | why it is uninformative |
|---|---|---|
| matched filter, time-blocked CV | true 0.26 / null **0.33** | decodes time: drift makes a test block resemble nearby training blocks, so any temporally blocky label scores |
| + local detrend + 32 PCs | unchanged | drift was not the only route; 0.33 is the retained-window majority rate, i.e. the true chance level |
| unbiased d'^2, 8-shift null | \|z\| < 2 everywhere, median z ~ 0 | sensitivity floor is d' ~ 1, about 100x too coarse for the d' ~ 0.1 that matters |
| image domain, orbit wedge | ROI/control RMS 1.5–1.8 | **identical for all four gaze pairs** (1.53 on sub-015 for up/down, left/right, up/left, down/right) — that is a regional noise property, not a gaze effect |

A mask-orientation bug invalidated the first image-domain pass: `h5py` returns
MATLAB arrays with reversed axes, so `masks.mat` landed in the middle of the
head rather than on the orbits. Fixed in `effect_size.masks()`; the numbers
above are post-fix.

## What this does and does not establish

**Does not establish that the effect is absent.** The orbit `roiMask` is a
loose wedge of ~465 000 voxels at 240^3, while the moving structure — the lens
— occupies of order 500. Wedge-wide statistics dilute a real effect ~1000x, and
attempts to segment the globe automatically landed on the sinuses instead. The
largest smoothed between-bin difference in the wedge sits on the **scalp edge**,
which is where bulk head motion shows up, not eye motion.

**The stimulus geometry is favourable**, which is the encouraging part.
From `visual_stimuli/mreyetrack_4points.py`: screen 369.54 mm wide, 800x600,
viewed at 1020 mm, targets at norm offsets 0.5 (vertical) and 1.33*2/3
(horizontal). That is **+/-3.9 deg vertical and +/-6.9 deg horizontal**, so
7.8 deg up-to-down and 13.7 deg left-to-right. At a 12 mm globe radius the lens
travels **1.6 mm and 2.9 mm** respectively — well above the 1 mm voxel. The
effect should be resolvable in principle.

## The experiment that would actually settle it

Forward-simulate rather than probe. Take the four binned reconstructions as the
best available estimates of the gaze-state images, NUFFT them onto the real
trajectory to get noiseless per-readout k-space per state, and compute the
per-readout d' against the noise level measured from the data. No classifier,
no CV, no class priors, no drift — it returns the required window size as a
number. Do this before training anything.

Cheaper first step: segment the globes properly (the existing `roiMask` is not
a globe mask) and re-measure the between-bin difference there. If the lens
moves visibly between bins, the probes above were simply insensitive; if it
does not, that is worth knowing before investing in a model.

## Cross-subject transfer is a structural problem, not a sample-size one

The ROI-PCA transform is estimated per subject, so virtual channel 3 of
sub-001 has no relationship to virtual channel 3 of sub-004. A model trained
across subjects must be invariant to that: a per-subject input adapter, or a
first layer invariant to rotations of the channel basis. Adding subjects does
not help until a within-subject model works.

---

# Positive control, take 2 — it passes (weakly)

`model/decode_gaze.py`, all 11 retained subjects, both bin types.

Two fixes over `probe_gaze*.py` made the result interpretable:

- **balanced** 4-way accuracy, so chance is 0.25 by construction rather than the
  retained-window majority rate (~0.33) that made the earlier runs unreadable;
- a **distribution** of 8 circular shifts, giving a per-subject z rather than a
  single number to eyeball.

| window | clean: mean z | p (t-test) | subjects positive | filtered: mean z | p |
|---:|---:|---:|---:|---:|---:|
| 0.5 s | +1.13 | 0.034 | 9/11 | +0.87 | 0.064 |
| 1.0 s | +1.05 | 0.024 | 10/11 | +1.05 | 0.022 |
| 2.0 s | +1.02 | 0.018 | 10/11 | +0.98 | 0.039 |
| **3.1 s** | **+1.03** | **0.015** | **10/11** | **+1.08** | **0.022** |
| 5.1 s | +0.93 | 0.028 | 10/11 | +1.12 | 0.011 |
| 8.2 s | +0.74 | 0.151 | 7/11 | +0.61 | 0.225 |

Sign test at 1–5 s: 10/11 positive, p = 0.012. Balanced accuracy is 0.28–0.30
against 0.25 chance.

**Read this as "there is signal, and it is small".** A nearest-mean template on
32 PCs is close to the weakest possible decoder, so a learned nonlinear model
should do better — but nobody should expect this task to be easy.

Three things the sweep says about the design:

- **Window size barely matters between 0.5 s and 5 s.** Accuracy is flat, which
  means the per-readout scores are strongly correlated in time and pooling adds
  little independent evidence. Use ~3 s because it is the best-controlled point,
  not because longer is better.
- **8.2 s collapses** because gaze blocks last 5 s, so almost no window reaches
  70% purity. That is a labelling ceiling, not a signal ceiling.
- **`clean` and `filtered` perform the same.** `filtered` has ~20% more readouts
  per class (74 655 vs 61 404 on sub-015) and both are disjoint 4-class
  partitions, so prefer `filtered` for training data volume and keep `clean` as
  the cleaner evaluation set.

Subject spread is large: sub-013 is consistently strong (z +2.3 to +4.0),
sub-001 is consistently **negative** (z -1.1 to -2.2). Do not tune on one
subject.

## Answers to the data questions

- `kspace_ROI-PCA_8_woBin.mat` contains **`y`, `t`, `ve`, `meta` — no `C`.**
  The matching compressed coil maps are a separate file in the same folder,
  `C_rovir_8.mat` (48^3 x 8, 6 MiB). `C` is only needed to reconstruct, not to
  train.
- **`y` is already stripped of the redundant readouts.** The reader applied
  `nSeg-1` (SI projection) and `nShot-nShotOff` before export: 43 x 1858 =
  **79 894** readouts. Nothing further to remove.
- **`t` and `ve` are bit-identical across all 15 subjects** (same `.seq`), and
  both are stored float64: `t` is 0.86 GiB and `ve` 0.29 GiB of every 3.2 GiB
  file, i.e. a third of the payload duplicated 15 times. Store them once, in
  float32, and the dataset shrinks by ~17 GB.

---

# The architecture, and why each piece is there

`model/train_gaze.py`. Input is a sliding window of W readouts; output is
per-readout class probabilities, so a scan comes back as an array of readouts
labelled up/down/left/right. Extra classes (blink, saccade) drop in by widening
the final head.

```
[W, 160 samples, 8 ch] complex ──> 16 real planes
        │
        ├── per-subject channel adapter   (16x16, one per subject)
        │
   ReadoutEncoder : Conv1d over the 160 k-space samples
        │          + FiLM modulation from the spoke direction
        │                                          -> e_r in R^128, per readout
   TemporalTrunk  : dilated TCN over the W readouts
        │
   head           : Conv1d 1x1  -> 4 logits per readout
```

## FiLM — Feature-wise Linear Modulation

A **conditioning mechanism, not an architecture** (Perez et al., 2018). A side
input — here the spoke direction — goes through a small MLP that emits a scale
`gamma` and shift `beta` *per feature channel*, applied between conv layers as
`x -> gamma * x + beta`.

FiLM is not itself a CNN. The encoder around it is: a `Conv1d` stack running
along the 160 k-space samples of each readout. FiLM sits between those conv
layers and modulates their channels.

**Why it is needed here.** A readout's k-space profile is dominated by which way
its spoke points; the eye is a small perturbation on top. A plain CNN would
spend most of its capacity learning "what does a spoke pointing this way look
like" before reaching the part we care about. FiLM lets the direction
*reparameterise* the filters, giving a smoothly varying filter response per
direction.

It is the continuous version of what the linear probe does by hand: that fits
600 independent per-direction-cluster templates, while FiLM learns direction ->
filter as a smooth function and shares statistical strength across directions.
Cost is three small `Linear` layers.

## Dilated TCN — Temporal Convolutional Network

1D convolutions along the **readout (time)** axis. *Dilated* means the kernel
taps skip positions: kernel 3 at dilation 8 reads `t-8, t, t+8`. Doubling the
dilation each layer (1, 2, 4, 8, 16) grows the receptive field **exponentially
with depth** rather than linearly, so a few cheap layers span a long stretch.

**Why it is needed here.** One readout carries almost nothing (per-readout
d' <= 0.1), while a gaze state persists ~5 s = 625 readouts, so evidence has to
be integrated across many readouts. The linear probe integrated with a fixed sum
of per-readout scores; the TCN learns the integration instead — it can weight
readouts unequally, use context, and detect transitions (which is how
saccades/blinks will eventually be read off).

Chosen over an RNN or transformer because it is parallel over the sequence,
stable to train, cheap, and has an **exactly computable receptive field** —
which is how the bug below was caught. Given that decoding accuracy is flat from
0.5 s to 5 s, a transformer (global receptive field) is worth trying to test
whether long-range context helps at all.

## First training run — a negative result

30 epochs, 11 subjects, 3.1 s windows, `filtered` bins:

| | train loss | test balanced acc |
|---|---:|---:|
| epoch 0 | 1.42 | 0.250 |
| epoch 10 | 0.65 | 0.236 |
| epoch 29 | **0.105** | 0.256 |

Best test over all epochs 0.2695 against 0.25 chance — **worse than the linear
probe's 0.28–0.30**. The network memorised the training set and generalised
nothing. When a nearest-mean template beats a CNN, the answer is a *smaller*
model, not a bigger one.

Three causes, in fix order:

1. **Windows overlap 6x** (stride 64, width 384), so 9 548 "training windows"
   are more like ~1 500 independent ones across 11 subjects — far fewer than the
   parameter count supports.
2. **Receptive-field bug.** Dilations 1,2,4,8,16 give `1 + 2*(1+2+4+8+16) = 63`
   readouts = 0.50 s, but the window fed in is 384 readouts = 3.07 s. **The
   trunk sees only 16% of its own window.** Reaching 3 s needs dilations to
   ~128, i.e. three more layers.
3. **No validation split.** Test accuracy was printed every epoch, which is
   peeking; the "best epoch" figure is not trustworthy. Needs a 3-way split.

---

# Forward-simulation calibration — attempted, does not work as designed

`model/calibrate.py`. Plan was: push the four binned reconstructions through the
real forward model (multiply by each virtual coil map, NUFFT onto the actual
spokes), measure noise per k-space sample from the acquired data, and read off
d' per readout. k-space is the right domain because the noise there is white and
per-sample, unlike the spatially correlated noise of an undersampled radial
reconstruction.

## Two things it revealed

**1. The forward model reproduces magnitude but not phase.** Best orientation
gives magnitude correlation 0.85 against the acquired k-space, but complex
correlation is ~0 at every `|k|` except exactly at k = 0. Removing **one phase
per readout** lifts it 0.47 -> 0.63, which points at acquisition physics absent
from a plain NUFFT — LIBRE's off-resonant excitation, RF phase cycling, eddy
currents. This is *not* fatal for the measurement: a per-readout phase is common
to all four gaze states and cancels in `K[g1] - K[g2]`. Validating on
magnitude is therefore legitimate.

Note this also holds when validating against `x0_ROI-PCA_8_woBin.mat`, a
reconstruction made from exactly this compressed data with this trajectory — so
it is a genuine gap in the forward operator, not a file-convention mistake.
Redoing the projection with **monalisa's own NUFFT in MATLAB** would remove
every remaining convention question by construction.

**2. The result is impossible, and that is the finding.** Median d' per readout
came out at 107–152 across the six gaze pairs. A d' of 125 means a *single*
readout classifies gaze essentially perfectly. The direct measurement on real
data says d' ~ 0.1, and the best windowed decoder reaches balanced accuracy 0.29
using 384 readouts. The calibration is off by ~3 orders of magnitude.

The cause is not the NUFFT. It is the input: **the difference between two binned
reconstructions is dominated by their own reconstruction noise**, not by gaze.
Each bin uses ~19% of the readouts, so its noise is ~2.3x the full recon's, and
differencing two of them adds another sqrt(2). Forward-projecting that noise
produces a per-sample difference of ~11% of the signal, which is what generates
the absurd d'. The six pairs span only 107–152 (ratio 1.43) where geometry
predicts left/right (2.9 mm of lens travel) should clearly exceed up/down
(1.6 mm) — weakly consistent in ordering, far too flat in magnitude to be
signal-dominated.

**The binned reconstructions cannot serve as gaze-state templates.** That is a
real constraint on the project, not just on this script.

## The version that would work: a rotated-eye phantom

Avoid noisy templates entirely. Take the **`woBin` reconstruction** — the
highest-SNR image available, built from all 79 894 readouts — segment the globe,
and *synthesise* the gaze states by rotating the eye region through the known
angles (7.8 deg between up and down, 13.7 deg between left and right, from
`visual_stimuli/mreyetrack_4points.py`). Forward-project the original and the
rotated volumes.

The difference is then **pure signal with zero reconstruction noise**, and the
resulting d' is exactly the ceiling a perfect model could reach. That is the
number worth having, and nothing else in this document supplies it.

---

# Shrunken model, and the feature-based alternative

## Shrunken model (`model/train_small.py`) — matches the probe, does not beat it

Encoder 16/32/64 (was 32/64/128), d=64, dropout 0.3, 170k params; dilations now
auto-sized so the TCN receptive field (511 readouts) exceeds its window (256);
real 3-way 60/20/20 split with the model selected on val and test touched once.

| | value |
|---|---|
| best val balanced acc | 0.2717 (epoch 35) |
| **TEST balanced acc** | **0.2833** (chance 0.25) |
| linear probe reference | 0.28–0.30 |

An honest number this time, and it lands exactly on the linear probe. Two very
different estimators hitting the same 0.28–0.30 suggests that is the ceiling of
what this representation supports, not a tuning failure.

## Does binning-before-compression matter? No

The gaze bins are readout masks; ROI-PCA acts on the channel axis
(`y_virt = Vret.' y`). Selecting rows and mixing columns **commute**, so masking
the compressed `woBin` export gives bit-identical data to compressing the binned
mitosius. Only `mitosius/ROI-PCA_8/woBin/` exists on disk, so there is also no
way to accidentally feed uncompressed binned data.

## Feature-based approach (`model/features/`) — CBMS 2025 style

`extract.py` computes 23 statistics x {real, imag, magnitude, phase} x 8
channels = **736 features per readout** (the paper: 23 x 4 x 20 = 1840 on 20
physical coils). `evaluate.py` runs leave-one-subject-out with per-subject
standardisation and class balancing, on two tasks kept deliberately separate.

**Why the two tasks are not equivalent.** Motion corrupts a readout's phase
directly, and that mechanism does not depend on where the eye points — so
direction-blind summary statistics can see it. Gaze *direction* needs the eye's
position resolved, a sub-millimetre structural difference, and summary
statistics collapse the sample axis and know nothing about the spoke direction.
Expect this approach to work on motion and struggle on gaze.

**This updates the recommendation in section 4 above.** That section ranked
fixation-vs-motion as the hardest place to start. The CBMS result is evidence
the direct motion route works well, because motion has a signature that does not
require resolving position. Motion detection is plausibly the *easier* task.

### First look, sub-015, motion vs still (1260 vs 74905)

Best genuine single-feature AUC **0.529**, and the measure composition of the
top features runs re 30–40% / mag 20–32% / im 16–27% / **phase 12–24%** against
a 25% no-preference baseline. **Phase is the weakest measure here**, the
opposite of the paper's finding that 53 of 59 selected features were phase.

Two caveats before reading anything into that: this is univariate AUC on one
subject, not the paper's multivariate ANOVA/MI/RF/DT selection across subjects.

**A likely explanation, and it matters:** the paper used 20 **physical** coils
and reported that no single channel dominated — the information was spread
across the array. ROI-PCA collapses 52 physical coils into 8 virtual ones by
keeping only the eye-ROI signal subspace, mixing coil phase in the process. The
compression may destroy exactly what the CBMS features exploit. **To reproduce
that work, run it on the uncompressed data, not the ROI-PCA export.**

`mag_zero_cross` is degenerate by construction (a magnitude cannot cross zero)
and is dropped along with 7 other constant features.

---

# Direction-anisotropy test — the prediction fails

`model/anisotropy.py`, all 11 subjects, per-cluster discriminability averaged
across subjects (clusters are shared because the trajectory is bit-identical,
so this buys sqrt(11)).

A spoke measures the FT of the object's 1D projection onto its own direction,
and a projection is blind to displacement perpendicular to it. Real gaze
evidence must therefore concentrate: left/right in spokes along the L-R axis,
up/down in spokes along S-I, and the two must pick **different, orthogonal**
axes. The test measures dependence on all three array axes without assuming
which is which.

Result: **no concentration anywhere.** All z against the shifted-label null fall
between +0.4 and +1.7 — uniformly weakly positive, never significant, and with
no monotone rise across the quartiles of any axis. Left/right shows a mild trend
on axis2 (z 0.6 -> 1.6) but up/down shows nothing on any axis, and the
prediction requires both.

The test has some power: if the whole effect were concentrated in one quartile
it should read z ~ 2 there and ~0 elsewhere. Instead it is spread evenly.

**The 0.28-0.30 that both the linear probe and the CNN reach is isotropic in
spoke direction, which is not what eye position looks like.**

## The alternative that fits every observation

Small **head** motion accompanying gaze shifts. People move their head slightly
when they look left or right, and that motion is time-locked to the gaze label,
so a circular-shift null does not remove it. It is a whole-head effect, so it is
**isotropic in spoke direction** — exactly what was measured. It also explains
why accuracy is flat from 0.5 s to 5 s (bulk pose persists across a whole gaze
block) and why the effect survives every control aimed at drift.

Cheap next diagnostic: split the band by `|k|`. Bulk motion dominates **low**
`|k|` (whole-head scale); eye structure lives at `|k|` 0.08-0.33. If the
decodable evidence sits entirely at low `|k|`, it is the head, not the eye.

---

# Two decisive results: the |k| split, and the feature approach on compressed data

## |k| split — does NOT confirm the head-motion hypothesis

`model/kband_split.py`, 4-way gaze, 256-readout windows, 11 subjects, 4-shift null:

| band | median acc | median z | subjects z>0 |
|---|---:|---:|---:|
| low `\|k\|` < 0.06 (>16 mm, bulk head) | 0.2895 | 1.67 | 10/11 |
| mid `\|k\|` 0.06–0.15 (7–16 mm, globe) | 0.2763 | 1.41 | 10/11 |
| high `\|k\|` 0.15–0.33 (3–7 mm, lens/ONH) | 0.2777 | 1.56 | 9/11 |

Bulk head motion would put the evidence at **low** `|k|` and leave high `|k|` at
chance. High `|k|` is not at chance, and the three bands are within noise of
each other. So the head-motion hypothesis from the previous section is **not**
supported in its strong form.

Combined with the anisotropy result, the effect is uniform in **both** spoke
direction and `|k|` — it has no spatial structure at all. Neither an eye
displacement (high `|k|`, direction-dependent) nor a head displacement (low
`|k|`) looks like that. What does: a **global per-readout modulation** —
amplitude or phase applied to the whole readout — that happens to be time-locked
to the stimulus. Receiver/B0 drift, respiration, or postural muscle tension
entrained to the 5 s target cadence would all qualify.

Caveat on power: 1.67 vs 1.41 vs 1.56 is well inside noise, so this rules out
the strong version of the head hypothesis, not every version.

## CBMS features on ROI-PCA data — chance, and not because of transfer

`model/features/evaluate.py motion`, leave-one-subject-out, 11 subjects,
balanced, per-subject standardised:

| | median balanced accuracy (chance 0.500) |
|---|---:|
| HistGradientBoosting | 0.4984 |
| LinearSVC | 0.4995 |

Every subject lands between 0.487 and 0.506. The obvious explanation was that
ROI-PCA is a per-subject basis, so nothing transfers. **It is not that.**
`model/features/within_subject.py`, time-blocked 70/30 inside each subject:

| | median balanced accuracy |
|---|---:|
| HistGradientBoosting, within subject | 0.5018 (range 0.491–0.528) |

At chance within subject too. **The motion information is not present in the
ROI-PCA-compressed data**, as these features capture it.

This lines up with the earlier univariate result, where phase was the *weakest*
of the four measures (12–24% of top features against a 25% baseline) while CBMS
found 53 of its 59 selected features were phase. ROI-PCA mixes 52 physical coils
into 8 virtual ones by keeping only the eye-ROI signal subspace, and coil phase
is mixed in the process. **The compression appears to destroy exactly what the
method relies on.**

## Consequence: the compressed route is closed for motion detection

Everything needed for the uncompressed route is already on disk:

| | |
|---|---|
| raw 52-channel `.dat` | `data/study/sub-NNN/rawdata/*_T1wLIBRE.dat`, 16.6 GB each |
| Python twix reader | `old_study/code/twixtools/` |
| **the CBMS cohort itself** | `old_study/raw_data/` — 44 GB, the 2018–2019 subjects, with `old_study/et_data/` alongside |

That last row matters: the paper's own data is here. **Reproduce the 98% on it
first.** That validates the feature pipeline against a known result before it is
pointed at this cohort's uncompressed data, and it settles whether the gap is
the compression or something else in the implementation.

---

# Uncompressed 52-channel features — also chance. The compression was not the cause.

`model/features/extract_raw.py` reads the raw `.dat` with twixtools (**exactly
82 368 image scans** = nSeg x nShot, acquisition order, 52 x 480 complex64) and
computes 23 stats x {re, im, mag, phase} x 52 channels = **4784 features per
readout**, 1.42 GiB per subject, ~20 min each.

## The pipeline is verified, so the null result is trustworthy

| control | result |
|---|---:|
| predict segment group (= spoke angle) from the same features, same split | **0.9865** (chance 0.250) |
| mean \|corr\| between putative SI navigators (index % 44 == 0) | **0.9997** |
| same-segment non-navigator readouts (index % 44 == 7) | 0.8349 |
| random readouts | 0.6205 |

0.9997 only happens if those really are the same k-space line every shot, which
confirms the (shot, segment) order with segment fastest and confirms the drop
removes exactly the SI navigators. And a pipeline that predicts spoke angle at
98.65% is not broken.

## Results, three subjects (015, 013, 009)

| test | compressed nv=8 | uncompressed 52 ch |
|---|---:|---:|
| motion, within subject | 0.5018 | 0.5062 (0.4705 / 0.5062 / 0.5070) |
| motion, leave-one-subject-out | 0.4984 | 0.4997 |
| gaze, within subject | — | 0.2488 (0.2488 / 0.2501 / 0.2419) |

sub-009 contributes 10 804 motion readouts (~3200 in the test split, SE ~0.009),
so this is properly powered, not a small-sample artefact.

**The ROI-PCA hypothesis is dead.** Uncompressed data gives chance too, so the
compression was not what stood between this cohort and the CBMS result.

## What is left

1. **Implementation** — my 23 statistics may differ materially from theirs.
2. **Labels** — theirs were ET labels *refined* by a PCA+GMM data-driven step
   that excluded 37.9% of the original motion samples; mine are the ET masks at
   `th0.75_winLen10` with transitional readouts dropped. That sliding window
   smooths over a 62 ms span, which can blur a 30–80 ms saccade.
3. **Cohort/acquisition** — theirs: 2018–19, 20 coils, FoV 192, TR 6.4 ms,
   TE 2.94 ms, 3723 x 22. Ours: 2026, 52 coils, FoV 240, TR 8.0 ms, TE 3.6 ms,
   1872 x 44.

Only one experiment separates these: **run this pipeline on
`old_study/raw_data/` (the paper's own cohort, 44 GB, on disk).** 98% there means
the implementation is right and the difference is cohort or labels; chance there
means the implementation differs from theirs. Nothing else distinguishes them,
and every further model built on this cohort is a guess until it is resolved.

---

# Compression comparison — the channel map is not the lever

4-way gaze, 256-readout windows, 11 subjects, 4-shift null, everything
downstream of the channel map held identical:

| variant | median balanced acc | median z | subjects z>0 |
|---|---:|---:|---:|
| ROI-PCA nv=8 | 0.2895 | 1.69 | 10/11 |
| SVD nv=8 | 0.2873 | 1.80 | 10/11 |
| SVD nv=16 | 0.2873 | 1.89 | 9/11 |

Trained model (`train_small.py`, identical architecture and protocol):

| variant | best val | TEST balanced acc |
|---|---:|---:|
| ROI-PCA nv=8 | 0.2717 | **0.2833** |
| SVD nv=8 | 0.2697 | **0.2728** |

**Every representation lands at 0.27–0.29.** Eye-ROI-targeted (ROI-PCA),
whole-FoV (SVD), 8 channels or 16 — no difference, and doubling the channel
count buys nothing. Combined with the isotropy and `|k|`-uniformity results,
the weak effect is insensitive to *everything* about the spatial representation.

`train_small.py` now takes `--variant`, so all three run through one code path
rather than forked copies.

## Feature route: closed at every channel count

| features | motion within-subject | motion LOSO | gaze within-subject |
|---|---:|---:|---:|
| ROI-PCA nv=8 (736 feats) | 0.5018 | 0.4984 | — |
| uncompressed 52 ch (4784 feats) | 0.4992 | **0.5000** | 0.2518 |

All 11 subjects, chance throughout, while the same features predict spoke angle
at 0.9865. Per-readout summary statistics do not capture eye state in this
cohort at any channel count.

---

# Step 4: FiLM/TCN on uncompressed 52 channels — no better

296k params (adapter widened to 104x104), same protocol, 42 s/epoch.

| representation | linear probe | trained model (TEST) |
|---|---:|---:|
| ROI-PCA nv=8 | 0.2895 | 0.2833 |
| SVD nv=8 | 0.2873 | 0.2728 |
| SVD nv=16 | 0.2873 | — |
| **uncompressed 52 ch** | — | **0.2557** |

More channels did not help; it was slightly worse. Train loss fell 1.45 -> 0.33
while val peaked at 0.2978 (epoch 17) and decayed — the extra capacity bought
overfitting against the same 4103 training windows, nothing else.

## Read the model numbers with their error bars

The test split is 682 windows over 4 classes, so macro balanced accuracy carries
SE ~0.017. **0.2833, 0.2728 and 0.2557 are all within noise of each other and
none is individually distinguishable from 0.25.** The only real evidence that
any effect exists is the 11-subject linear probe against a circular-shift null
(median z ~1.0-1.2, sign test p = 0.012), and that is weak.

## The representation axis is exhausted

Across everything tried, the effect is unchanged by:

| varied | result |
|---|---|
| spoke direction | isotropic (anisotropy test, all z +0.4..+1.7) |
| `\|k\|` band | uniform (0.2895 / 0.2763 / 0.2777) |
| window length 0.5–5 s | flat |
| coil subspace (eye-ROI vs whole-FoV) | identical |
| channel count 8 / 16 / 52 | identical |
| per-readout summary features | chance at every channel count |

A genuine eye-position signal should care about at least the first two. This one
cares about none of them, which is the signature of a **non-spatial, weakly
label-locked nuisance** rather than gaze.

## What is worth doing instead

Stop adding representational capacity — six variations have now produced the
same number. The unanswered question is not *how* to model the signal but
*whether there is one to model*, and one experiment settles it:

**the rotated-eye phantom.** Take the `woBin` reconstruction (highest SNR, all
79 894 readouts), segment the globe, and synthesise the gaze states by rotating
the eye region through the known angles (7.8 deg up-down, 13.7 deg left-right).
Forward-project the original and rotated volumes onto the real trajectory. The
difference is pure signal with zero reconstruction noise, so the resulting d'
per readout is the exact ceiling a perfect model could reach.

The earlier calibration attempt failed only because it used the *binned*
reconstructions as templates and their differences are dominated by their own
noise. The phantom removes that, and it is the only measurement that can
distinguish "the task is infeasible at this SNR" from "we have not found the
right model".

---

# Where the signal is lost: everything upstream is correct, the contrast is not

Landmark coordinates clicked by the user in ITK-SNAP for all four gaze bins
(lens and optic nerve head, both eyes), converted with the empirically derived
mapping `axis0 = 239-y, axis1 = 239-x, axis2 = z` (correlation 1.000000).

## The eye tracker: the eye really did move, and the bins are clean

| | measured from ET | stimulus geometry |
|---|---:|---:|
| left vs right excursion | **18.00 deg** | 13.75 deg |
| up vs down excursion | **7.90 deg** | 7.78 deg |

Within-bin gaze scatter is only 9–15 px IQR (0.2–0.4 deg), so the bins are
**pure**, not contaminated. Implied displacement at an 11 mm globe radius
(confirmed by the clicked lens-to-ONH distance of 22.0 / 22.8 mm):
lens 3.46 / 1.52 mm, **lens-minus-ONH differential 6.9 / 3.0 mm**.

## The temporal alignment is correct

ET target changes map onto the readout bins with the expected ~0.5 s saccade
latency, and with the left-right mirror the scanner's mirror introduces:

| ET message | time | readout bin that switches on | time |
|---|---:|---|---:|
| dot moved up | 30.5 s | up | 31.0 s |
| dot moved left (screen) | 35.5 s | **right** | 36.0 s |
| dot moved down | 40.5 s | down | 40.9 s |
| dot moved right (screen) | 45.5 s | **left** | 45.9 s |

ET recording 663.1 s, MRI 82 368 x 8 ms = 658.9 s. Nothing is misaligned.

## The images show none of it

lens-minus-ONH differential, measured by sub-voxel cross-correlation:

| recon | L: L-R | L: U-D | R: L-R | R: U-D | ET implies |
|---|---:|---:|---:|---:|---:|
| STEVA (regularised) | -0.61 | +0.41 | -1.00 | +0.81 | 6.9 / 3.0 |
| **gridded x0 (none)** | **-0.14** | **+0.06** | **+0.03** | **-0.28** | 6.9 / 3.0 |
| measurement noise floor | ±0.3–0.9 | | | | |

The noise floor comes from repeating the measurement on `clean_0.50/0.75/0.95`,
which have identical bin purity and therefore identical true displacement.
Regularisation is **not** the cause — the unregularised gridded recons show even
less. Everything is within noise of zero, against a 3–7 mm prediction.

## The likely reason, and its caveat

Lens-to-vitreous contrast-to-noise in the gridded recon is **1.6 (left eye) and
0.4 (right eye)** — at or below the noise. A structure at CNR < 2 cannot have
its 3 mm displacement tracked, by cross-correlation or by a model.

*Caveat:* the clicked lens coordinate sits at axis0 = 41 with the globe centre
at ~52 and radius ~11, i.e. exactly the anterior pole — so that box may be on
the corneal surface rather than the lens body, and the measured CNR would then
understate the true lens contrast. Worth re-measuring 3–5 voxels posterior.

## What this explains

Every null in this document follows from one fact: **the moving structures do
not survive into these reconstructions at usable SNR.** That accounts for the
chance-level features at every channel count, the 0.28 ceiling that was
insensitive to spoke direction, `|k|` band, window length, coil subspace and
channel count, and the absence of any visible lens displacement between bins.

The upstream pipeline — eye tracking, bin purity, ET-to-readout alignment — is
sound. The loss is between the readouts and the image, and the question is now
an acquisition/contrast question, not a machine-learning one.

---

# The eye does rotate between bins: the optic nerve settles it

`model/nerve_profile.py`. Displacement measured in bands along each optic nerve,
gradient-weighted, using A-eye labels carried into LIBRE space by one
MPRAGE->woBin transform. Nerves split by connected component (454 / 419 vox) --
a midline split is wrong because the nerves converge medially toward the chiasm.

| distance from globe | left L-R | left U-D | right L-R | right U-D |
|---|---:|---:|---:|---:|
| 10–14 mm (ONH) | **+1.03** | **+0.60** | **+1.10** | **+0.64** |
| 14–20 mm | +0.37 | +0.41 | +1.02 | +0.45 |
| 20–29 mm | +0.20 | +0.10 | +0.11 | +0.13 |
| 29–38 mm (apex) | **+0.08** | **+0.05** | **+0.10** | **+0.04** |

**Monotone decreasing in all four series**, falling 10–20x from the optic nerve
head to the orbital apex, with correlation 0.75–0.89 throughout.

This is the tethered-nerve signature and it cannot be produced by any nuisance:

- **bulk head motion shifts every band equally** — the brain control reads
  0.00 mm and the apex bands read 0.04–0.10 mm, so there is no bulk motion to
  confuse it;
- the two nerves are independent structures and agree closely;
- **L-R exceeds U-D** on both sides (1.03/1.10 vs 0.60/0.64, ratio ~1.7)
  matching the stimulus ratio 18.0/7.9 = 2.3.

**The binning works and the eye rotates.** Everything upstream — eye tracking,
bin purity, ET-to-readout alignment, reconstruction — is doing its job.

## The remaining gap is magnitude, not existence

Measured 1.0–1.1 mm at the ONH against ~3.5 mm predicted from the ET excursion
at an 11 mm radius: about 30%, consistent across every estimator tried (clicks,
boxes, label ROIs, gradient weighting, corneal cap, nerve bands).

So the signal a decoder has to work with is ~1 mm of structure displacement, not
the ~3 mm the stimulus geometry implies — 3x smaller in amplitude, ~10x in
energy. That is the quantitative reason the k-space decoding sat at chance, and
it is a property of the data rather than of any model.

Note also the ET measured 18.0 deg where the stimulus geometry gives 13.75 deg,
a 31% overestimate suggesting a calibration gain issue; correcting it narrows
the prediction to ~2.6 mm but does not close the gap.

---

# The per-readout signal, measured directly — and it is not there

`model/empirical_dprime.py`. The forward-simulation ceiling (`ceiling.py`) kept
returning d' of 8-25 per readout, which would make a single readout classify
gaze perfectly. Every decode we ran says otherwise, so the simulation is wrong,
not the data. This measures the same quantity with no forward model at all.

The trajectory revisits directions, and those revisits land ~80 s apart, hence
usually in different gaze states. So take readout pairs pointing the same way
(within 1.0 deg; the median nearest-neighbour spacing is 0.71 deg) and compare
`||y_a - y_b||^2` for same-bin against cross-bin pairs. Same-bin pairs carry
noise, drift and everything that is not gaze; cross-bin pairs carry that plus
the gaze difference.

sub-015, 135 097 labelled pairs (34 520 same-bin, 100 577 cross-bin):

| | median `\|\|y_a - y_b\|\|^2` |
|---|---:|
| same-bin | 3.3739e11 |
| cross-bin | 3.3412e11 |
| **ratio** | **0.9903** |

| contrast | pairs | ratio |
|---|---:|---:|
| up/down | 17 190 | 0.9834 |
| left/right | 16 598 | 0.9414 |
| up/left | 12 676 | 0.9902 |
| down/right | 19 214 | 0.9941 |

**Cross-bin pairs are no more different than same-bin pairs** — every ratio is
at or below 1. There is no measurable per-readout gaze difference in the
k-space. The bound is roughly d' < 0.1 per readout, and probably far below it.

## This reconciles everything

The eye **does** rotate — the nerve profile establishes that beyond doubt in the
image domain. But an image is reconstructed from ~15 000 readouts, so per-readout
the evidence is ~1/122 of what the image shows. At d' < 0.1, reaching AUC 0.9
needs W > 330 readouts and, at the effect sizes the decodes actually exhibited,
far more than the 5 s (625-readout) gaze block allows.

So both facts are true at once: **the binning works and the eye moves, and the
per-readout k-space signal is still too small to classify.** That is a property
of the acquisition, not of any model or representation, and it explains every
null in this document.

## Compression comparison, reconstructed

`model/recon_compare.py`: adjoint NUFFT per virtual channel then RSS, identical
treatment for both 6.5x compressions.

| | eye-ROI mean | head mean | background sigma | eye-ROI SNR |
|---|---:|---:|---:|---:|
| ROI-PCA nv=8 | 0.440 | 0.758 | 0.0896 | 4.9 |
| SVD nv=8 | 0.448 | 0.791 | 0.0798 | 5.6 |

Visually indistinguishable; the difference map is a smooth low-frequency gain
pattern, not structure. Note SVD gives slightly *better* RSS SNR despite
retaining only 16.6% of eye-ROI energy against ROI-PCA's 99.3% — RSS ignores
coil sensitivities, so it favours the channels carrying the most total signal.
With a proper sensitivity-weighted recon ROI-PCA should win in the orbit.

---

# CORRECTION: the 98% we were chasing was an artifact

**Everything in this document that treats the CBMS 2025 result as a target to
reproduce is wrong, and the sections below supersede it** — specifically
"Feature-based approach (CBMS 2025 style)", "CBMS features on ROI-PCA data",
"Uncompressed 52-channel features", and every mention of reproducing 98%.

Per the author (2026-09-02): that classifier was not separating motion from
no-motion. It was separating **redundant periodic MR structure — the SI
projection, repeating across shots and subjects**. Re-run with the correct
ET-derived classes, its balanced accuracy is **~50%**.

## What this changes

**Our null results are corroboration, not failure.** Chance at 8 channels,
chance at 52, chance with per-readout statistics, chance with raw k-space and a
temporal model — all agree with the corrected CBMS number. There was never a 98%
to reach.

**It also explains a result of ours that should have raised an eyebrow.** The
positive control in `evaluate_raw.py` predicted the *segment group* — i.e. the
spoke angle — at **0.9865** from the same features. That is the same class of
periodic structure the CBMS classifier was picking up. It was recorded here as
proof the feature pipeline worked; it was equally a demonstration of the
artifact, and the coincidence of ~98% in both places was not noticed at the time.

Our pipeline was never exposed to it: segment 0 of every shot (the SI navigator)
is dropped before anything else, so the periodicity could not leak into the
labels. That is luck of construction, not a safeguard that was designed in.

**Withdrawn as a consequence:** the recommendation to run the feature pipeline on
`old_study/raw_data/` to "reproduce the 98% first". There is nothing to
reproduce. Likewise the hypothesis that ROI-PCA compression destroyed the phase
information the method relied on — that was invented to explain a gap that does
not exist.

## What still stands

Independent of CBMS, and unaffected by this correction:

- the eye genuinely rotates between bins (optic nerve profile, monotone 10–20x
  falloff from ONH to apex, null control);
- the per-readout k-space difference between gaze states is consistent with
  zero (95% CI on the cross/same pair ratio entirely below 1);
- the 0.27–0.29 gaze decode ceiling, insensitive to spoke direction, `|k|`,
  window length, coil subspace and channel count.

## The one defect this leaves worth fixing

The motion labels come from a sliding window of `winLen=10, th=0.75`, i.e. a
readout is labelled only if **75% of an 80 ms window** agrees. A saccade lasts
30–80 ms, so a 40 ms saccade reaches at most 0.5 overlap and **can never be
labelled at any position**. Blinks (100–400 ms) pass easily. That matches the
counts: sub-015 has 350 saccade readouts against 910 blink, and the cohort
saccade fraction ranges 0.1%–7.8%, implausible as physiology.

`winLen=10, th=0.75` is the right rule for **reconstruction binning**, where a
conservative purity filter is what you want. Reusing it as an **ML label** is the
error; the two requirements are opposite.

| winLen | window | th | saccade readouts (40 ms event) | blink (200 ms) |
|---:|---:|---:|---:|---:|
| 10 (current) | 80 ms | 0.75 | **0** | ~15 |
| 3 | 24 ms | 0.50 | ~2–4 | ~22 |
| 1 | 8 ms | any overlap | ~5 | ~25 |

Prerequisite before rebuilding: **the ET-MRI synchronisation must be known to
~10 ms**. The block-level check (target change 30.5 s -> bin onset 31.0 s)
bounds it at 0.5 s including saccade latency, which is not a measurement of sync
error. An 8 ms label is meaningless if the sync is uncertain at 100 ms.

---

# ET-MRI synchronisation: start offset is 6 ms, drift is unresolved

The PsychoPy log carries the scanner trigger, which the ET log does not: the
scanner emits `Keypress: s` every ~2.5 s and PsychoPy waits on it.

| event | PsychoPy clock |
|---|---:|
| "waiting for trigger" displayed | 420.208 s |
| `eyetracker.setRecordingState(True)` | 727.279 s |
| **first scanner trigger** | **727.285 s** |

**The ET recording started 6.0 ms before the first scanner trigger.** PsychoPy
had been listening for 307 s beforehand, so that trigger is the scan start, not
a mid-scan pulse. The binning code assumes MRI t=0 == ET t=0; at the start of
the scan that assumption is good to **6 ms, well under one TR (8 ms)**.

## Drift is the open question, and it is the binding constraint

Fitting trigger index against time gives a period of **2.502568 +/- 0.000029 s**
(271 triggers, 675.2 s span, residual rms 37.6 ms consistent with 60 Hz keyboard
polling). If the scanner's nominal period is exactly 2.5000 s, that is a
**+1027 ppm clock error = +677 ms accumulated over a 659 s scan**. The nominal
period could not be confirmed from the available files.

An independent check using ET target changes against readout-bin onsets gives
the *opposite* sign (-1.5 ms/s, -995 ms over the scan), but that measure is
dominated by saccade latency and by the sliding window's own filtering -- early
lags run 2.9-4.4 s, which is behaviour, not timing. It does not settle the
question.

So drift is bounded loosely somewhere within +/-1 s over the scan, i.e. up to
~1.5 ms per second. **At winLen=1 (8 ms) even 0.5 ms/s accumulates to 40 TRs of
error by the end of the acquisition**, which would scramble per-readout labels
in the second half of every scan while leaving 80 ms labels roughly intact.

## How to settle it without new acquisitions

The ET pupil-size trace carries a cardiac oscillation, and the twix header
carries PMU timestamps from the scanner's pulse recording. Cross-correlating the
two in sliding windows across the scan measures the relative timing directly, at
several points in time, in the two clocks that actually matter. That is the
measurement to make before choosing a new `winLen`.

Until then: **6 ms at the start is established, drift is not**, and shortening
the label window is only safe to the extent that drift is known.

## RESOLVED: the drift was my arithmetic, not the hardware

The scanner records the external trigger itself. `PMUTimeStamp` in the twix mdh
is a sawtooth (262 resets over the scan), and its resets are the trigger times
**in the scanner's own clock** -- which makes the comparison direct instead of
inferred.

| clock | fitted trigger period |
|---|---:|
| scanner (`PMUTimeStamp` resets vs `TimeStamp`) | 2.502583 s (residual rms 3.0 ms) |
| stimulus PC (PsychoPy `Keypress: s`) | 2.502568 s |
| **ratio** | **1.00000609 = +6 ppm** |

**Drift over a 659 s scan: +4 ms.** The earlier +1027 ppm figure came entirely
from assuming the nominal period was exactly 2.5000 s. It is not -- the scanner
says 2.5026 s, and both clocks agree on that to 6 ppm.

### Total synchronisation budget

| term | value |
|---|---:|
| start offset (ET recording vs first trigger) | +6 ms |
| drift, scanner vs stimulus PC, whole scan | +4 ms |
| **total** | **~10 ms, about one TR** |

**Per-readout (8 ms) labels are therefore viable.** The timing objection to
shortening `winLen` is removed.

### The one link still inferred

The ET mask timeline is EyeLink *samples*, so the EyeLink's own 1 kHz clock is
what ultimately matters, and it was measured only against PsychoPy. The
`dot moved` messages fit an interval of 5.01367 s in EyeLink clock against a
PsychoPy-commanded 5 s rendered at a measured 59.77 Hz -- which is 299, 300 or
301 frames, i.e. 5.0025 / 5.0192 / 5.0360 s. The frame count is not recorded, so
this leg is bounded only to roughly +/-2000 ppm (up to ~1 s over the scan) rather
than measured.

To close it: count EyeLink samples between two `dot moved` messages directly in
the `.tsv.gz` rather than trusting message timestamps, or read the EDF's own
sample clock. Worth doing before committing to `winLen = 1`; `winLen = 3`
(24 ms) already recovers saccades and tolerates this residual.

## The EyeLink leg, measured: +33 ppm

The `.tsv.gz` rows *are* the 1 kHz EyeLink samples (the notebook uses
`len(recording)` as milliseconds directly), so the sample count is an exact
EyeLink-clock measurement of the recording interval. Anchoring it against
PsychoPy at both ends:

| | |
|---|---:|
| EyeLink samples recorded | 663 051 |
| PsychoPy `setRecordingState` True -> False | 663 028.9 ms |
| **ratio** | **+33 ppm** -> **+22 ms over a 659 s scan** |

(Do **not** use the `eye tracker stopped` message for this -- it fires 16 s after
recording actually ends and gives a nonsensical +22 619 ppm.)

Uncertainty is dominated by the ET start/stop command latency at each end, which
is ~7 ms judging by `setRecordingState(True)` at PsychoPy 727.2791 against the
EyeLink's own `start recording` 7 ms later. At +/-5 ms per end that is +/-15 ppm,
so the honest figure is **33 +/- ~15 ppm, i.e. 22 +/- 10 ms across the scan**.

### Full synchronisation budget

| term | value | how established |
|---|---:|---|
| start offset, ET vs first trigger | +6 ms | PsychoPy log, direct |
| drift, scanner vs stimulus PC | +4 ms | PMU resets vs keypresses, +6 ppm |
| drift, EyeLink vs stimulus PC | +22 +/- 10 ms | sample count vs interval, +33 ppm |
| **net EyeLink vs scanner, end of scan** | **~24 ms (~3 TRs)** | |

## Verdict on winLen: use 3, not 1

`winLen = 1` is an 8 ms window, and the residual misalignment reaches ~24 ms by
the end of the acquisition -- three readouts. A 40 ms saccade is only 5 readouts
long, so labels would be right early in the scan and smeared by more than half
an event late in it, in a way that varies systematically with scan time. That is
a worse failure mode than a slightly wide window, because it is *structured*.

`winLen = 3` (24 ms, th 0.50) recovers **2-4 readouts per saccade against zero
at the current setting**, and its window is the same size as the worst-case
timing error, so the error costs at most one window rather than several.

`winLen = 1` becomes defensible only after applying a drift correction --
`et_index = readout_time_ms * (1 + 27e-6)`, the EyeLink-vs-scanner rate implied
above -- and even then the +/-10 ms residual is comparable to one TR.

---

# Cohort synchronisation audit — winLen must be 3, and two subjects need checking

`model/sync_audit.py`, all 15 subjects, from the PsychoPy log (`Keypress: s`
triggers, `setRecordingState`) and the `.tsv.gz` sample count.

| sub | offset | EyeLink ppm | EL drift | net | |
|---|---:|---:|---:|---:|---|
| 001 | +11.1 ms | +33 | +21.6 | 32.7 | |
| 002 | +6.6 | +36 | +23.8 | 30.4 | |
| 003 | +4.8 | -6 | -4.2 | **9.0** | best |
| 004 | +7.6 | -6 | -4.2 | 11.8 | |
| 005 | +7.0 | -61 | -40.3 | 47.3 | |
| 006 | +7.7 | -55 | -36.3 | 44.0 | |
| 007 | **-506** | +118 | +77.5 | 583.6 | offset suspect |
| 008 | +12.9 | +113 | +74.4 | **87.3** | worst drift |
| 009 | **-509** | +23 | +14.9 | 523.9 | offset suspect |
| 010 | +6.9 | +29 | +19.1 | 26.0 | |
| 011 | +5.9 | -92 | -60.6 | 66.5 | |
| 012 | +11.1 | -5 | -3.3 | 14.4 | |
| 013 | +7.4 | -3 | -1.8 | **9.2** | best |
| 014 | **-2955** | +34 | +22.5 | 2977.8 | offset suspect |
| 015 | +6.0 | +33 | +22.0 | 28.0 | |

Median net 32.7 ms; excluding the three suspects the range is **9-87 ms**.

## Two traps found while measuring this

**Do not "clean" the trigger train by requiring regular 2.5 s intervals.** The
*first* interval is systematically short (~1.99 s), so that filter discards the
genuine first trigger and reports a spurious ~2 s offset on every subject.

**Trigger-train alignment cannot fix the offset either** (`model/sync_align.py`).
The train is periodic at 2.5 s, so any whole-period shift aligns equally well:
it placed sub-015 at +2002 ms where the direct reading gives +6 ms, a one-pulse
error. Residuals of ~4 ms confirm the *rates* agree; the absolute pulse identity
does not follow.

## The three suspects

sub-007, 009 and 014 each have exactly **one trigger before**
`setRecordingState(True)`. Either that pulse started the scan (so the ET began
0.5-3.0 s late and those bins are badly misaligned), or it is a stray from a
preceding sequence and the true offset is ~8 ms. The periodic ambiguity above
means the trigger train cannot decide. **sub-009 and sub-014 are in the 11-subject
training set**, so this is not academic: if the offsets are real, ~18% of that
data carried scrambled gaze labels.

Resolving it needs a non-periodic anchor -- e.g. `PMUTimeStamp` at the very
first readout gives the phase since the preceding trigger, which pins the scan
start within one period rather than modulo one.

## Verdict

`winLen = 3` (24 ms) covers the median 33 ms only after the per-subject drift
correction, which is now measurable from files every subject already has:
`et_index = readout_ms * (1 + ppm*1e-6)` with the per-subject ppm above. With
that correction the residual is the start offset alone, ~5-13 ms, and winLen 3
has comfortable margin.

**Without the correction, winLen = 5 (40 ms)** is the honest choice, since four
subjects (005, 006, 008, 011) exceed 40 ms of drift on their own. `winLen = 1`
is not defensible for any subject except 003 and 013.

---

# Would a better ROI mask (A-eye) help the compression? No — and here is the bound

`model/pairtest_raw.py`. The question was whether feeding the A-eye segmentation
into ROI-PCA, instead of the hand-drawn wedge, would retain gaze information the
current compression throws away.

Compression is a **linear map**, so it can only lose information. That means the
compressed null result cannot by itself prove the information was never present.
Running the same pair test on the **uncompressed 52 channels** bounds what *any*
channel recombination — any ROI mask, any ranking criterion — could deliver.

| sub | representation | features | ratio | 95% CI | |
|---|---|---:|---:|---|---|
| 015 | ROI-PCA nv=8 | 1280 | 0.9903 | [0.9824, 0.9973] | no excess |
| 015 | **uncompressed 52** | 8320 | **0.9861** | [0.9803, 0.9932] | **no excess** |
| 013 | ROI-PCA nv=8 | 1280 | 0.9421 | [0.9336, 0.9501] | no excess |
| 013 | **uncompressed 52** | 8320 | **0.9501** | [0.9440, 0.9567] | **no excess** |
| 010 | ROI-PCA nv=8 | 1280 | 0.9800 | [0.9676, 0.9945] | no excess |
| 010 | **uncompressed 52** | 8320 | **0.9608** | [0.9499, 0.9731] | **no excess** |

Ratio 1.000 would mean cross-bin readout pairs differ from same-bin pairs by
exactly the gaze signal. Every confidence interval sits **entirely below 1** —
uncompressed as well as compressed. There is no per-readout gaze excess to
recover, so **no ROI mask can help**.

This closes the question from the decoding side, and matches what
`recon/ROVir/README.md` already found from the compression side: shrinking the
ROI 29x (wedge 3432 voxels -> one globe, 120) moved retained energy by ~2 points.

## Three independent reasons, now consistent

1. ROI size was already measured as a weak lever (ROVir README).
2. **SVD retains 16.6% of eye-ROI energy against ROI-PCA's 99.3% — a 6x
   difference — and decodes identically** (0.2873 vs 0.2895 probe; 0.2728 vs
   0.2833 trained).
3. The uncompressed bound above: the signal is absent before any compression.

Note the ratios sit slightly *below* 1 rather than at 1. Cross-bin pairs are
marginally more similar than same-bin pairs, most likely a time-gap confound
(same-bin pairs span longer intervals and accrue more drift). It does not affect
the conclusion, which only requires the excess not to be positive.
