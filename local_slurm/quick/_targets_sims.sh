#!/bin/bash

#SBATCH --mem=4g
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/_targets_sims_quick.out
#SBATCH --error=logs/_targets_sims_quick.err

cd /home/naimrashid/Downloads/DeSurv-paper

Rscript --max-connections=1012 local_slurm/quick/submit_targets_sims.R
