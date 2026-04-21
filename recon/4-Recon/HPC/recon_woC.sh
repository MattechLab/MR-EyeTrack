#!/bin/bash
#SBATCH --job-name=mreyetrack_recon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon/4-Recon/HPC/logs/%x_%j.out
#SBATCH --error=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon/4-Recon/HPC/logs/%x_%j.err
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jaime.barrancohernandez@hevs.ch

echo "Job ID: $SLURM_JOB_ID" 

SIF="/home/jaime.barrancohernandez/shared_datasets/monalisa/monalisa_251215.sif"
MATLAB_CMD="addpath(genpath('/usr/src/app')); compile_mex_for_monalisa; run('/usr/src/app/scripts/4-Recon/HPC/S02_chuv_woC_chacha.m'); exit;"

export OMP_NUM_THREADS=1

apptainer exec \
  --bind /home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon:/usr/src/app/scripts \
  --bind /home/jaime.barrancohernandez/mnt/jaime.barranco/MR-EyeTrack/data/study:/usr/src/app/data/study \
  --bind /home/jaime.barrancohernandez/shared_datasets/pulseq:/usr/src/app/pulseq \
  --writable-tmpfs \
  --env MLM_LICENSE_FILE=27000@matlablm.hevs.ch \
  $SIF \
  bash -c "chmod -R ugo+x /usr/src/app/scripts && matlab -nodisplay -nosplash -r \"$MATLAB_CMD\""
