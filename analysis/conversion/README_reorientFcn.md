# Deriving the `reorientFcn` for a new dataset

`mat2nii_twix` needs a `reorientFcn` (7th argument) whenever a new sequence is
added. This document is the procedure for deriving it.

Budget ~15 minutes per new sequence. The result is stable for every subject
acquired with that sequence, so it is a one-time cost per sequence, not per scan.

---

## 1. What the function has to achieve

`mat2nii_twix` takes the direction cosines from the **reference** NIfTI and
writes the `reorientFcn` output *verbatim*:

```matlab
refAffine = niftiinfo(refNifti).Transform.T';
dirCos    = refAffine(1:3,1:3) ./ vecnorm(refAffine(1:3,1:3), 2, 1);
affine(1:3,1:3) = dirCos .* voxelSize;      % ← reference's axes
niftiwrite(vol, ...)                         % ← your array, unmodified
```

So the job of `reorientFcn` is:

> **Permute and flip the raw reconstruction array until its dim1/dim2/dim3 run
> along the same physical directions as the reference's dim1/dim2/dim3.**

**The target is the reference's storage order — not RAS.** These coincide only
when the reference is RAS-stored. Getting this wrong is what produced the
left-right mirrored T2w volume: the fix was assumed to be on dim 1 (the L-R axis
in RAS) when for that reference L-R is dim 3.

### Find the target frame first

```bash
python3 -c "import nibabel as nib; print(nib.aff2axcodes(nib.load('REF.nii.gz').affine))"
```

or in MATLAB, read the normalised direction cosines and label each column by its
dominant component:

```matlab
refInfo = niftiinfo(refNifti);
A = refInfo.Transform.T';
R = A(1:3,1:3) ./ vecnorm(A(1:3,1:3), 2, 1);
disp(R)   % column k = direction that reference dim k runs along, in RAS
```

Known references in this project:

| Reference | Storage order | Target for `reorientFcn` |
|---|---|---|
| MR-EyeTrack MPRAGE | `('R','A','S')` | (+R, +A, +S) — plain RAS |
| Yiwei 2.0 MPRAGE (dataset 0005) | `('P','I','L')` | (−A, −S, −R) |

Write the target down before touching the data.

---

## 2. Produce the unmodified volume

Convert **with no `reorientFcn`** so the array is written exactly as the
reconstruction produced it:

```matlab
mat2nii_twix(matFile, twixFile, seqFile, refNifti, outNii, metaFile);
```

Do not start from a header-derived guess. A non-identity starting transform makes
every rotation you observe a *composition* of your guess and the true error,
which is materially harder to invert than starting from identity. (The
twix-header/trajectory route was tested in
`analysis/test_orientation/mat2nii_auto_orient_test.m`; for isotropic 3D radial
it carries no axis-assignment information — `dRowDir`/`dColDir` are absent from
the headers and the trajectory covariance is uniform at ≈0.33 on the diagonal.)

---

## 3. Work out which raw dim is which

Open the file in Mango (RAS, neurological, axial main view) or ITK-SNAP.

The viewer chooses which plane to draw from the **affine**, and right now the
affine is the reference's while the data is raw. So the slot that slices raw
dim *k* is fixed by the reference, and is known before you look at anything:

- the slot showing **sagittal** slices the reference dim whose direction is ±R
- the slot showing **coronal** slices the reference dim whose direction is ±A
- the slot showing **axial** slices the reference dim whose direction is ±S

| Reference | sagittal slot slices | coronal slot slices | axial slot slices |
|---|---|---|---|
| `('R','A','S')` | raw dim1 | raw dim2 | raw dim3 |
| `('P','I','L')` | raw dim3 | raw dim1 | raw dim2 |

### 3a. Identify the axes

For each slot, note **what anatomy actually appears in it**. A slot slices along
one raw dim, so the anatomy tells you that dim's true identity:

| Anatomy seen in the slot | ⇒ that raw dim is |
|---|---|
| sagittal (profile, corpus callosum, brainstem) | an R axis |
| coronal (two hemispheres, temporal lobes low) | an A axis |
| axial (ventricles, eyes at front, symmetric) | an S axis |

### 3b. Identify the signs

The in-plane rotation of each slot gives the signs. For a slot displaying
correctly you expect:

| Slot | right of image | top of image |
|---|---|---|
| axial | +R | +A |
| coronal | +R | +S |
| sagittal | +A | +S |

Read off the rotation and apply it. "Rotated 90° to the left" means the content
appears 90° **counter-clockwise** from the above, so what should be at the top is
now at the left:

- 180° → right = −(expected right), top = −(expected top)
- 90° CCW → right = −(expected top), top = +(expected right)
- 90° CW → right = +(expected top), top = −(expected right)

Combine 3a and 3b into a signed triple for the raw array, e.g. `(−A, −R, +S)`
meaning: raw dim1 increases toward posterior, dim2 toward left, dim3 toward
superior.

---

## 4. Build the function

Given raw axes and the target frame from step 1:

1. **Permute** — for each target position, find which raw dim holds that axis.
   `permute(v, [p1 p2 p3 (4:ndims(v))])` where `p_k` is the raw dim that belongs
   in position `k`. Always carry the `(4:ndims(v))` tail so 4D (multi-bin)
   volumes survive.
2. **Flip** — after permuting, flip every position whose sign disagrees with the
   target's sign at that position.

```matlab
reorientFcn = @(v) flip(flip(permute(v, [p1 p2 p3 (4:ndims(v))]), f1), f2);
```

### Worked example — MR-EyeTrack

Raw `(−A, −R, +S)`, target `(+R, +A, +S)`.

- +R lives in raw dim2, +A in raw dim1, +S in raw dim3 → `permute(v,[2 1 3 ...])`
- after permute: `(−R, −A, +S)` → signs wrong at positions 1 and 2 → `flip(...,1)`, `flip(...,2)`

```matlab
@(v) flip(flip(permute(v, [2 1 3 (4:ndims(v))]), 1), 2)
```

### Worked example — Yiwei T2w

Raw `(−A, −R, +S)`, target `(−A, −S, −R)` (reference is P,I,L).

- −A is raw dim1, −S is raw dim3 (as +S, sign to fix), −R is raw dim2 (already −R)
  → `permute(v,[1 3 2 ...])`
- after permute: `(−A, +S, −R)` → only position 2 disagrees → `flip(...,2)`

```matlab
@(v) flip(permute(v, [1 3 2 (4:ndims(v))]), 2)
```

---

## 5. Verify

Re-run `mat2nii_twix` with the function and confirm every slot shows the correct
plane with no residual rotation or flip. `check_nifti_orientation.m` reorients by
the affine for display, so it should now agree with the reference.

---

## 6. Check left-right — this step is not optional

**Visual inspection cannot detect a left-right mirror.** A mirrored brain looks
entirely plausible, and neither axis identity (3a) nor rotation (3b) changes
under a mirror — so a mirrored result passes steps 3–5 cleanly. Mask-overlap
metrics such as Dice are blind to it too.

Only intensity comparison against a reference of known handedness resolves it:

```bash
python3 analysis/test_orientation/tissue_check/lr_flip_test.py
```

It registers the volume and a mirrored copy to a reference and compares MI/CC.
Prefer a same-session reference in the same contrast family over the MPRAGE — for
the T2w the MPRAGE margins were −0.006/−0.027 (ambiguous) while the same-session
T1w LIBRE gave −0.164/−0.046 (decisive). See
`analysis/test_orientation/tissue_check/README.md`.

---

## Confirmed conventions

Raw axes are expressed **relative to each dataset's reference storage order**.

| Dataset | Reference | Raw axes | `reorientFcn` |
|---|---|---|---|
| MR-EyeTrack LIBRE | R,A,S | (−A, −R, +S) | `@(v) flip(flip(permute(v,[2 1 3 (4:ndims(v))]),1),2)` |
| Yannick AudioBOLD | R,A,S | (−R, −A, +S) | `@(v) flip(flip(v,1),2)` |
| Yiwei 2.0 T1w LIBRE (MID00030) | P,I,L | (+R, −A, +S) | `@(v) flip(flip(permute(v,[2 3 1 (4:ndims(v))]),2),3)` |
| Yiwei 2.0 T2w LIBRE (MID00025) | P,I,L | (−A, −R, +S) | `@(v) flip(permute(v,[1 3 2 (4:ndims(v))]),2)` |

Sanity check on these: the physical gradient frame predicts `(−R, −A, +S)`.
MR-EyeTrack and Yiwei T2w both differ from it by a Gx↔Gy swap, which is exactly
what the `swap1` in both sequence filenames denotes; Yiwei T1w differs by a
single dim1 sign. Earlier revisions of this table listed the Yiwei raw axes
against an assumed RAS target and did not reconcile with each other.
