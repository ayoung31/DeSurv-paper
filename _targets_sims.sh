#!/bin/bash

#SBATCH --mem=64g
#SBATCH -t 04:00:00
#SBATCH --output=logs/_targets_sims.out
#SBATCH --error=logs/_targets_sims.err

module load r/4.4.0

Rscript --max-connections=1012 submit_targets_sims.R