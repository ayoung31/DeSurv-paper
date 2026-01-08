#!/bin/bash

#SBATCH --mem=16g
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/_targets_quick.out
#SBATCH --error=logs/_targets_quick.err

Rscript --max-connections=1012 local_slurm/submit_targets_quick.R
