#!/usr/bin/env bash
# Stage ROVir inputs on the filer, submit the recon on chacha, pull results back.
#
#   ./push_rovir.sh -s 15 -n 20 -b woBin --submit    # stage + sbatch
#   ./push_rovir.sh -s 15 -n 20 -b woBin --pull      # copy results into the repo
#   ./push_rovir.sh -s 15 -n 20 -b woBin --submit --no-stage
#   ./push_rovir.sh --scripts-only                   # just refresh the .m/.sh on chacha
#
# The filer is one share seen from both machines, so the mitosius is written
# ONCE and stays there - no per-job transfer:
#
#   this workstation  /mnt/filer01/MatTechLab/jaime.barranco/MR-EyeTrack/data/study
#   chacha            ~/mnt/jaime.barranco/MR-EyeTrack/data/study
#   inside apptainer  /usr/src/app/data/study      (recon_rovir.sh binds it)
#
# Both are //filer01.hevs.ch CIFS mounts of fs_projets/MatTechLab. Only the
# scripts go over ssh, because chacha keeps those on its own disk under
# shared_datasets/ rather than on the filer.
set -euo pipefail

HOST=chacha
SUBJECT=15
NV=20
BIN=woBin
VARIANT=ROVir
DO_SUBMIT=0
DO_PULL=0
SCRIPTS_ONLY=0
NO_STAGE=0
DRYRUN=""

LOCAL_REPO="/home/debi/jaime/repos/MR-EyeTrack"
FILER_DATA="/mnt/filer01/MatTechLab/jaime.barranco/MR-EyeTrack/data/study"
REMOTE_SCRIPTS="/home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon"

usage() { sed -n '2,9p' "$0"; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--subject)   SUBJECT="$2"; shift 2 ;;
    -n|--nv)        NV="$2";      shift 2 ;;
    -b|--bin)       BIN="$2";     shift 2 ;;
    -v|--variant)   VARIANT="$2"; shift 2 ;;
    --submit)       DO_SUBMIT=1;  shift ;;
    --pull)         DO_PULL=1;    shift ;;
    --scripts-only) SCRIPTS_ONLY=1; shift ;;
    --no-stage)     NO_STAGE=1;   shift ;;
    --dry-run)      DRYRUN="--dry-run"; shift ;;
    -h|--help)      usage ;;
    *) echo "unknown arg: $1"; usage ;;
  esac
done

SUBJ=$(printf "sub-%03d" "$SUBJECT")
echo "host=$HOST subject=$SUBJ nv=$NV bin=$BIN variant=$VARIANT"

if [[ ! -d "$FILER_DATA" ]]; then
  echo "ERROR: filer not mounted at $FILER_DATA" >&2
  echo "It is an autofs mount - 'ls /mnt/filer01' usually wakes it." >&2
  exit 1
fi

# --- scripts (ssh; chacha keeps these off the filer) ------------------------
echo "==> refreshing ROVir scripts on $HOST"
ssh "$HOST" "mkdir -p $REMOTE_SCRIPTS/ROVir/hpc"
rsync -a $DRYRUN "$LOCAL_REPO/recon/ROVir/" "$HOST:$REMOTE_SCRIPTS/ROVir/"
[[ $SCRIPTS_ONLY -eq 1 ]] && { echo "scripts only, done."; exit 0; }

RD_LOCAL="$LOCAL_REPO/data/study/$SUBJ/recon"
RD_FILER="$FILER_DATA/$SUBJ/recon"
VDIR="$VARIANT"

# --- results back into the repo --------------------------------------------
if [[ $DO_PULL -eq 1 ]]; then
  echo "==> pulling results from the filer"
  mkdir -p "$RD_LOCAL/$VDIR"
  rsync -av $DRYRUN \
    --include="x_steva_${VARIANT}_${NV}_*" --include="x0_${VARIANT}_${NV}_*" --exclude="*" \
    "$RD_FILER/$VDIR/" "$RD_LOCAL/$VDIR/"
  exit 0
fi

# --- stage inputs on the filer ---------------------------------------------
if [[ $NO_STAGE -eq 0 ]]; then
  echo "==> staging coil maps, masks, transform"
  mkdir -p "$RD_FILER/$VDIR"
  rsync -av $DRYRUN \
    "$RD_LOCAL/$VDIR/C_rovir_${NV}.mat" \
    "$RD_LOCAL/$VDIR/masks.mat" \
    "$RD_LOCAL/$VDIR/rovir_transform.mat" \
    "$RD_FILER/$VDIR/"

  # Ship x0 when it exists so the cluster skips the Mathilda step
  X0="$RD_LOCAL/$VDIR/x0_${VARIANT}_${NV}_${BIN}.mat"
  [[ -f "$X0" ]] && { echo "==> staging precomputed x0"; rsync -av $DRYRUN "$X0" "$RD_FILER/$VDIR/"; }

  # rsync -c would reread both sides over CIFS; size+mtime is enough here and
  # skips the 6.5 GB copy entirely on repeat runs.
  echo "==> staging mitosius (skipped if already identical on the filer)"
  mkdir -p "$RD_FILER/mitosius/${VARIANT}_${NV}/${BIN}"
  rsync -av $DRYRUN --info=progress2 \
    "$RD_LOCAL/mitosius/${VARIANT}_${NV}/${BIN}/" \
    "$RD_FILER/mitosius/${VARIANT}_${NV}/${BIN}/"
fi

if [[ $DO_SUBMIT -eq 1 ]]; then
  echo "==> submitting"
  ssh "$HOST" "cd $REMOTE_SCRIPTS/ROVir/hpc && \
    sbatch --export=ALL,SUBJECT=$SUBJECT,NV=$NV,BIN=$BIN,VARIANT=$VARIANT recon_rovir.sh"
  echo "watch with:  ssh $HOST squeue -u \$USER"
fi
