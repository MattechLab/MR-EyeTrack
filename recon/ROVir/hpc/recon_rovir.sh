#!/bin/bash
#SBATCH --job-name=rovir_recon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=220G
#SBATCH --time=03:00:00
#SBATCH --output=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon/4-Recon/HPC/logs/%x_%j.out
#SBATCH --error=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon/4-Recon/HPC/logs/%x_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=jaime.barrancohernandez@hevs.ch

# ROVir reconstruction on chacha/disco.
#
# Mirrors recon/4-Recon/HPC/recon_woBin_steva.sh. Unlike that one, this needs no
# raw .dat in the container: S4_rovir_recon_chacha.m takes FoV as a parameter,
# so only the ROVir mitosius and C_rovir_<nv>.mat have to be pushed.
#
# Usage:
#   sbatch recon_rovir.sh                       # sub-015, nv=20, woBin
#   sbatch --export=ALL,SUBJECT=15,NV=20,BIN=mask_0 recon_rovir.sh
#
# Push the inputs first with push_rovir.sh.

echo "Job ID: $SLURM_JOB_ID"

SUBJECT="${SUBJECT:-15}"
NV="${NV:-20}"
BIN="${BIN:-woBin}"
VARIANT="${VARIANT:-ROVir}"
echo "subject=$SUBJECT nv=$NV bin=$BIN variant=$VARIANT"

SIF="/home/jaime.barrancohernandez/shared_datasets/monalisa/monalisa_251215.sif"
MATLAB_CMD="addpath(genpath('/usr/src/app')); compile_mex_for_monalisa; \
subject_num=${SUBJECT}; nv=${NV}; binName='${BIN}'; variant='${VARIANT}'; \
run('/usr/src/app/scripts/ROVir/hpc/S4_rovir_recon_chacha.m');"

export OMP_NUM_THREADS=1

apptainer exec \
  --bind /home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon:/usr/src/app/scripts \
  --bind /home/jaime.barrancohernandez/mnt/jaime.barranco/MR-EyeTrack/data/study:/usr/src/app/data/study \
  --bind /home/jaime.barrancohernandez/shared_datasets/pulseq:/usr/src/app/pulseq \
  --writable-tmpfs \
  --env MLM_LICENSE_FILE=27000@matlablm.hevs.ch \
  $SIF \
  bash -c "chmod -R ugo+x /usr/src/app/scripts && matlab -batch \"$MATLAB_CMD\""
