#!/usr/bin/env bash
set -euo pipefail

BINNING_DIR="/home/debi/jaime/repos/MR-EyeTrack/recon/2-Binning"
LOG_FILE="/tmp/S2_batch_fine_masks.log"

: > "$LOG_FILE"

echo "Log: $LOG_FILE"

for mask in fixation_ok saccade_only blink_event tracking_loss; do
    SCRIPT="$BINNING_DIR/S2_eyeMask_t1_binning_pulseq_batch_${mask}.m"
    echo ""
    echo "=========================================="
    echo "Running: $mask"
    echo "=========================================="
    matlab -batch "diary('$LOG_FILE'); diary on; run('$SCRIPT'); diary off"
    echo "Done: $mask"
done

echo ""
echo "All four masks complete. Log: $LOG_FILE"
