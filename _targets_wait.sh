#!/bin/bash
#SBATCH --job-name=desurv_pipeline
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00

# Use SLURM_SUBMIT_DIR if available, otherwise script directory
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"

# Load R 4.5 and set user library path
module load r/4.5.0
export R_LIBS_USER=/proj/rashidlab/nur2/R_libs/4.5
export DESURV_CPU_LIMIT=200

Rscript -e 'targets::tar_make()'
