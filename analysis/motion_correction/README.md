# Rigid Motion Correction for Monalisa Mathilda Recon Pipeline

This repository provides two helper functions to perform **rigid-body motion correction** in a Monalisa reconstruction pipeline, using spm motion estimation.  

---

## Overview

The motion correction workflow is designed to be integrated into any existing reconstruction pipeline. 
It focuses on **applying motion corrections directly to k-space data and trajectories** after estimating motion from a sequential 4D reconstruction.

The main steps are:

1. **Estimate rigid motion parameters** from a sequential 4D reconstruction using SPM:
    ```matlab
    motionParams = estimate_rigid_motion_params_spm( ... );
    ```
2. **Apply rigid-body motion to raw k-space and trajectory**:
    ```matlab
    [y_corr, t_corr, meta] = apply_rigid_motion_to_kspace( ... );
    ```

---

## Assumptions

- You have already computed a **4D sequential reconstruction** (`x_cs`).
- Motion parameters are estimated using **SPM**.
- Motion is applied directly to **k-space data (`y`)** and **trajectory (`t`)**.
- The repository does **not include a full reconstruction pipeline**; these functions are helpers to integrate into your own pipeline.
- I provide an example script, that does not run, just to showcase where in the reconstruction pipeline the helper function are designed to be used.

---

## Files

- `estimate_rigid_motion_params_spm.m` – Computes rigid-body motion parameters from a 4D reconstruction using SPM.
- `apply_rigid_motion_to_kspace.m` – Applies motion corrections to raw k-space data and trajectories.
- `example_script.m` – Illustrates where and how to use the two functions in a standard reconstruction pipeline. **This script does not run as-is**; it is for illustration purposes only.
- `extract_seq_params.m` – Utility to extract parameters from Pulseq sequence files. **Probably will not apply to your .seq files, you might need to change this file.**

---

## Example Usage

```matlab
% Path to your 4D reconstruction
pathTo4DImage = '/path/to/x_cs.mat';

% Estimate motion parameters
motionParams = estimate_rigid_motion_params_spm( ...
    pathTo4DImage, ...
    'x_cs', ...
    [2 2 2 5], ...
    'timeseries');

% Apply rigid motion to raw k-space
[y_moco, t_moco] = apply_rigid_motion_to_kspace( ...
    y, t, timestampMs, motionParams, motionTimeMs );
```
## Installations needed
You need to have Monalisa, (Pulseq if you are using it for trajectory generation) and SPM installed


## Example Results
Already with a 3 minutes acquisition, where the participant was ask to stand still, we can observe an improvement of the anatomical details:
<img width="2360" height="969" alt="MOCO" src="https://github.com/user-attachments/assets/c2dc3888-a1a8-4c97-a353-d1a538aee59f" />



