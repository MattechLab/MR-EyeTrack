#!/bin/sh
#SBATCH --job-name=mreyetrack_recon
#SBATCH --nodelist=dance        # dance, disco, chacha
#SBATCH --ntasks=1               # total number of tasks across all nodes 
#SBATCH --cpus-per-task=20       # cpu-cores per task (>1 if multi-threaded tasks) 
#SBATCH --mem=480G               # total memory per node (4 GB per cpu-core is default) 
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS) 
#SBATCH --mail-type=begin        # send email when job begins 
#SBATCH --mail-type=end          # send email when job ends 
#SBATCH --mail-user=jaime.barrancohernandez@hevs.ch 
#SBATCH --output=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/data/study/logs/%x_%j.out 
#SBATCH --error=/home/jaime.barrancohernandez/shared_datasets/mreyetrack/data/study/logs/%x_%j.err 
#module purge 

srun apptainer exec \ 
--bind /home/jaime.barrancohernandez/shared_datasets/mreyetrack/recon:/usr/src/app/scripts \
--bind /home/jaime.barrancohernandez/shared_datasets/mreyetrack/data/study:/usr/src/app/data/study \
--writable-tmpfs --fakeroot \ 
--env MLM_LICENSE_FILE=27000@matlablm.hevs.ch \ 
/home/jaime.barrancohernandez/shared_datasets/monalisa/monalisa_251215.sif \ 
bash -c "chmod -R ugo+x /usr/src/app/scripts && matlab -batch 'addpath(genpath(\"/usr/src/app/\")); compile_mex_for_monalisa; S4_recon_chacha_4fr_woBin_pulseq; exit;'" 
