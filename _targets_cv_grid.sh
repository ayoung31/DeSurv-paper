#!/bin/bash

#SBATCH --mem=64g
#SBATCH -t 12:00:00
#SBATCH --output=logs/_targets_cv_grid.out
#SBATCH --error=logs/_targets_cv_grid.err

module load r/4.4.0

Rscript --max-connections=1012 submit_targets_cv_grid.R
