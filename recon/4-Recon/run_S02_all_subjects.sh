#!/bin/bash
# Run S02_chuv_woC.m for subjects 1-15 (per-coil gridding, no C).
# Outputs: data/study/sub-NNN/recon/woBin/x0_noC_48.mat
#
# Usage: bash recon/4-Recon/run_S02_all_subjects.sh

set -euo pipefail

SCRIPT="$(realpath "$(dirname "$0")/S02_chuv_woC.m")"
BASEDIR="$(realpath "$(dirname "$0")/../../data/study")"
LOGFILE=/tmp/run_S02_all_subjects.log
MATRIX_SIZE=48

echo "=== S02 batch started at $(date) ===" | tee -a "$LOGFILE"

for sub in $(seq 1 15); do

    SUBJECT_DIR="${BASEDIR}/$(printf 'sub-%03d' "$sub")"
    OUT="${SUBJECT_DIR}/recon/woBin/x0_noC_${MATRIX_SIZE}.mat"

    if [ -f "$OUT" ]; then
        echo "[skip] Subject ${sub}: output already exists" | tee -a "$LOGFILE"
        continue
    fi

    echo "--- Subject ${sub} started at $(date) ---" | tee -a "$LOGFILE"
    matlab -nodisplay -nosplash -batch \
        "subject_num=${sub}; saveflag=1; matrix_size=${MATRIX_SIZE}; run('${SCRIPT}')" \
        2>&1 | tee -a "$LOGFILE"
    echo "--- Subject ${sub} done at $(date) ---" | tee -a "$LOGFILE"

done

echo "=== S02 batch finished at $(date) ===" | tee -a "$LOGFILE"
