#!/bin/bash

#SBATCH --mem=16g
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/_targets_bladder.out
#SBATCH --error=logs/_targets_bladder.err

Rscript --max-connections=1012 submit_targets_bladder.R
