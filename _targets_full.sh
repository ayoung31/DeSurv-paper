#!/bin/bash
#SBATCH --job-name=desurv_full
#SBATCH --mem=8g
#SBATCH -t 2-00:00:00
#SBATCH --output=logs/_targets_full.out
#SBATCH --error=logs/_targets_full.err

echo "Starting DeSurv FULL pipeline at $(date)"
echo "Using store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"
Rscript --max-connections=1012 submit_targets_full.R
echo "Pipeline completed at $(date)"
