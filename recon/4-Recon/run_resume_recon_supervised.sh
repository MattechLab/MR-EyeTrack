#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/debi/jaime/repos/MR-EyeTrack"
RUNNER="$REPO_ROOT/recon/4-Recon/run_resume_recon_logged.sh"
COUNT_SCRIPT="$REPO_ROOT/recon/4-Recon/count_missing_resume_recons.sh"

SUP_LOG="/tmp/S4_resume_recon_supervisor.log"
SUP_PID="/tmp/S4_resume_recon_supervisor.pid"
CHILD_PID="/tmp/S4_resume_recon_worker.pid"
LOCK_DIR="/tmp/S4_resume_recon_supervisor.lock"

mkdir -p /tmp

timestamp() {
  date '+%Y-%m-%d %H:%M:%S %z'
}

log() {
  printf '[%s] %s\n' "$(timestamp)" "$*" | tee -a "$SUP_LOG"
}

cleanup() {
  rm -f "$SUP_PID"
  rm -f "$CHILD_PID"
  rmdir "$LOCK_DIR" 2>/dev/null || true
}

on_exit() {
  local exit_code=$?
  log "Supervisor exiting with status $exit_code"
  cleanup
  exit "$exit_code"
}

trap on_exit EXIT
trap 'log "Supervisor interrupted"; exit 130' INT TERM

if ! mkdir "$LOCK_DIR" 2>/dev/null; then
  echo "Another supervisor instance appears to be running: $LOCK_DIR" >&2
  exit 1
fi

printf '%s\n' "$$" > "$SUP_PID"
log "Supervisor started with pid $$"

while true; do
  missing="$("$COUNT_SCRIPT")"
  if [[ "$missing" == "0" ]]; then
    log "All recon outputs are present. Nothing left to run."
    break
  fi

  log "Missing recon outputs remaining: $missing"
  log "Launching resume runner"

  bash "$RUNNER" &
  child_pid=$!
  printf '%s\n' "$child_pid" > "$CHILD_PID"
  log "Runner pid: $child_pid"

  set +e
  wait "$child_pid"
  child_status=$?
  set -e

  rm -f "$CHILD_PID"

  missing_after="$("$COUNT_SCRIPT")"
  if [[ "$child_status" -eq 0 ]]; then
    log "Runner exited cleanly"
  else
    log "Runner exited with status $child_status"
  fi

  if [[ "$missing_after" == "0" ]]; then
    log "All recon outputs are now complete."
    break
  fi

  log "Still missing $missing_after outputs; restarting in 15 seconds"
  sleep 15
done
