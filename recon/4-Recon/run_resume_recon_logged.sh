#!/usr/bin/env bash
set -euo pipefail

LOG_FILE="/tmp/S4_resume_recon_debi_4fr_pulseq.log"
SCRIPT_FILE="/home/debi/jaime/repos/MR-EyeTrack/recon/4-Recon/S4_resume_recon_debi_4fr_pulseq.m"

: > "$LOG_FILE"

matlab -batch "diary('$LOG_FILE'); diary on; run('$SCRIPT_FILE'); diary off"
