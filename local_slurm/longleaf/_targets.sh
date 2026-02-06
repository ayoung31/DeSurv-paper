#!/bin/bash

#SBATCH --mem=32g
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/_targets.out
#SBATCH --error=logs/_targets.err

module load r/4.4.0

export DESURV_CPU_LIMIT=200

Rscript --max-connections=1012 submit_targets.R
