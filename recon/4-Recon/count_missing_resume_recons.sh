#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/home/debi/jaime/repos/MR-EyeTrack/data/study"

missing=0
for n in $(seq 1 15); do
  s=$(printf '%03d' "$n")
  for mask in clean clean_0.50 clean_0.75 clean_0.95; do
    for region in 0 1 2 3; do
      x_path="$BASE_DIR/sub-${s}/recon/${mask}/x/x_steva_regionidx_${region}_nIter_20_delta_1.000.mat"
      if [[ ! -f "$x_path" ]]; then
        missing=$((missing + 1))
      fi
    done
  done
done

printf '%s\n' "$missing"
