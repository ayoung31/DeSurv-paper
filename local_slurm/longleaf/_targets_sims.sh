#!/bin/bash

#SBATCH --mem=32g
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/_targets_sims.out
#SBATCH --error=logs/_targets_sims.err

module load r/4.5.0

export DESURV_CPU_LIMIT=200
export R_LIBS_USER=/proj/rashidlab/nur2/R_libs/4.5

Rscript --max-connections=1012 submit_targets_sims.R
