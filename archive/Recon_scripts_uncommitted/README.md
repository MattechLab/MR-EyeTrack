# Uncommitted `recon/Recon_scripts` changes (rescued)

Verbatim copies of three files that existed **only as uncommitted working-tree
changes** inside the `recon/Recon_scripts` submodule. They matched no blob in
any commit of that submodule — local or upstream — so they were one `git
checkout` away from being lost permanently. Saved here on 2026-09-23.

Provenance:

| | |
|---|---|
| submodule | `recon/Recon_scripts` |
| branch | `mac` (50 commits behind `origin/mac`) |
| checked-out commit | `e5cd8b7a13194626a1bc266f37ba3077e61d32cb` |
| file mtimes | 2026-01-19 → 2026-01-22 |

Paths below mirror the submodule layout.

| file | vs. submodule HEAD |
|---|---|
| `func/et_binning_related/eyeGenerateBinningWin.m` | ~150 lines changed |
| `archive/sop_check_mrtrack_local_movement/s1_checkCoil.m` | +45 lines |
| `archive/sop_check_mrtrack_local_movement/s2_extract_k_center.m` | +1 line |

## Why `eyeGenerateBinningWin.m` matters

Its signature was changed, dropping the interactive `user_input`/`input_info`
batch path:

```matlab
% submodule HEAD (7+2 args)
function cMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg, th_ratio, ETDir, winLen, display_binning, ...
    user_input, input_info)

% this copy (9 args)
function cMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg, th_ratio, ETDir, winLen, display_binning, save_binning, mask_type)
```

**This 9-argument form is the one the active pipeline calls.**
`recon/2-Binning/S2_eyeMask_t1_binning_pulseq.m:74` and
`S2_eyeMask_t1_binning_pulseq_nomo.m:68` both invoke it as:

```matlab
eMask = eyeGenerateBinningWin(rawDir, nShotOff, nSeg, th_ratio, ETDir, winLen, true, true, mask_type);
```

Note that `recon/Binning/eyeGenerateBinningWin.m` is a *different, older* copy
carrying the 7-arg signature. S2 does `addpath(genpath(...recon))` then
`addpath(genpath(...recon/Recon_scripts))`, and the second call prepends — so
the submodule working-tree copy shadows `recon/Binning/` and is what actually
runs. Restoring the submodule to its HEAD would therefore break S2 with a
"Too many input arguments" error.

These copies are deliberately placed **outside `recon/`** so `genpath` does not
add them to the MATLAB path and create a third shadowing candidate. They are a
preservation snapshot, not an importable source.
