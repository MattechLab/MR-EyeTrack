#!/usr/bin/env bash
# MRIQC (nipreps/mriqc:latest, v25.1.0.dev45) on the pulseq undersampling BIDS set.
# REGU (steva*) and GRID (mathilda*) run as two separate MRIQC datasets so their
# group tables stay apart. Work dirs are kept on /home on purpose: the brain
# tissue segmentation and the MNI transforms inside them are reused for the
# region masks (04_register_masks.py), and /tmp was wiped last year.
#
#   bash analysis/us_study_pulseq/02_run_mriqc.sh [wobin|regu|grid|all]
#
# `wobin` is a small extra pass on the 15 woBin REGU volumes only; its brain
# segmentation and MNI transform are the masks used by 04_register_masks.py,
# available long before the 255-volume REGU run finishes.
set -euo pipefail
ROOT=/home/debi/jaime/repos/MR-EyeTrack/results/us_pulseq
BIDS=$ROOT/bids/derivatives/reconstructions
CFG=$(cd "$(dirname "$0")" && pwd)/cfg
WHICH=${1:-all}
# MRIQC's SynthStrip node peaks at 12.5 GB RSS on a 240^3 volume and MRIQC does
# not budget for it: it runs nprocs/omp of them at once, and 6-8 concurrent were
# OOM-killed (exit 137) inside the 62.8 GB Docker VM. 24/8 caps it at 3.
NPROCS=${NPROCS:-24}
OMP=${OMP:-8}
MEMGB=${MEMGB:-50}
SUBS=$(seq -f "%03g" 1 15 | tr '\n' ' ')

# MRIQC aborts the whole run on the first crashed node (config.nipype.stop_on_first_crash
# defaults to True and has no CLI flag): one heavily undersampled volume failing its
# "inhomogeneity-corrected data seem empty" check killed a 255-volume run. The
# entrypoint below turns that off, so a failing volume just gets no IQMs.
# patches/anatomical.py: Harmonize falls back to the uneroded WM mask when erosion
# empties it. Without it, 17/510 LIBRE volumes (weak WM/GM contrast: only ~2k voxels
# reach p(WM) >= 0.9) got an all-zero harmonized image and no IQMs.
PATCH=$(cd "$(dirname "$0")" && pwd)/patches/anatomical.py
MRIQC_PY='import sys; from mriqc import config; config.nipype.stop_on_first_crash = False; from mriqc.cli.run import main; sys.argv[0] = "mriqc"; main()'

run() {  # $1 = wobin|regu|grid (cfg/filter_$1.json)
  local tag=$1
  mkdir -p "$ROOT/mriqc_$tag" "$ROOT/work_$tag"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp --entrypoint python \
    -v "$BIDS":/data:ro -v "$ROOT/mriqc_$tag":/out -v "$ROOT/work_$tag":/work -v "$CFG":/cfg:ro \
    -v "$PATCH":/opt/conda/lib/python3.11/site-packages/mriqc/interfaces/anatomical.py:ro \
    nipreps/mriqc:latest -c "$MRIQC_PY" /data /out participant --participant-label $SUBS -m T1w \
    --nprocs "$NPROCS" --omp-nthreads "$OMP" --mem-gb "$MEMGB" --no-sub --notrack \
    --bids-filter-file "/cfg/filter_$tag.json" -w /work
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$BIDS":/data:ro -v "$ROOT/mriqc_$tag":/out \
    nipreps/mriqc:latest /data /out group -m T1w --no-sub --notrack
}

case $WHICH in
  wobin) run wobin ;;
  regu) run regu ;;
  grid) run grid ;;
  all)  run wobin; run regu; run grid ;;
esac
