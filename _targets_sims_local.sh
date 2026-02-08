#!/bin/bash
#SBATCH --job-name=desurv_sims_full
#SBATCH --mem=28G
#SBATCH --cpus-per-task=19
#SBATCH --time=96:00:00
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

# Full simulation reproduction: 6 specs × 4 scenarios × 100 replicates = 2,400 runs
# Estimated time: ~60-80 hours with 2 concurrent workers

# Use SLURM_SUBMIT_DIR if available, otherwise script directory
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"

# Wait for system load to drop if needed
echo "Starting simulation pipeline at $(date)"
echo "Full reproduction: 2,400 simulation runs"
echo ""

# Load R
module load r/4.5.0 2>/dev/null || true
export R_LIBS_USER=/proj/rashidlab/nur2/R_libs/4.5

# Run the simulation pipeline
Rscript --max-connections=1024 -e '
library(targets)

# Use the same store as the main pipeline
tar_config_set(store = "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full")

cat("Starting full simulation pipeline\n")
cat("Store:", tar_config_get("store"), "\n")
cat("Time:", as.character(Sys.time()), "\n\n")

# Run the simulation targets
tar_make(script = "_targets_sims.R")

cat("\nSimulation pipeline completed at:", as.character(Sys.time()), "\n")
'

echo "Simulation pipeline finished at $(date)"
