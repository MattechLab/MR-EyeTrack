#!/usr/bin/env bash
# A-eye (nnU-Net Task313_Eye, 5-fold ensemble, CPU) on every subject's MPRAGE,
# both orbits via quadrant cropping, then the distance gate that drops the
# mouth/jaw fat components. Reuses the tissue_check tooling unchanged.
#
#   bash analysis/us_study_pulseq/03_run_aeye.sh [sub-001 ...]
set -euo pipefail
REPO=/home/debi/jaime/repos/MR-EyeTrack
TC=$REPO/analysis/test_orientation/tissue_check
OUT=$REPO/results/us_pulseq/aeye          # must stay under /home/debi (Docker mounts)
PY=/home/debi/miniconda3/envs/mreyetrack/bin/python
SUBS=${@:-$(seq -f "sub-%03g" 1 15)}

for sub in $SUBS; do
  mpr=$REPO/data/study/$sub/dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/$sub.nii.gz
  if [[ -f $OUT/$sub/${sub}_aeye_both_gated.nii.gz ]]; then echo "skip $sub"; continue; fi
  $PY "$TC/run_aeye.py" --input "$mpr" --out "$OUT/$sub"
  $PY "$TC/aeye_gate.py" --labels "$OUT/$sub/${sub}_aeye_both.nii.gz"
done
