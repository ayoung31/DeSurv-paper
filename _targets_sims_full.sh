#!/bin/bash
#SBATCH --job-name=desurv_sims_full
#SBATCH --mem=8g
#SBATCH -t 3-00:00:00
#SBATCH --output=logs/_targets_sims_full.out
#SBATCH --error=logs/_targets_sims_full.err

echo "Starting DeSurv SIMULATION pipeline (FULL) at $(date)"
echo "Using store: store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full"
echo "Analysis methods: bo, bo_alpha0, bo_tune_ntop, bo_tune_ntop_alpha0, fixed, fixed_alpha0"
echo "Scenarios: R0_easy, R00_null, R0k6, R_mixed (100 replicates each)"
Rscript --max-connections=1012 submit_targets_sims_full.R
echo "Simulation pipeline completed at $(date)"
