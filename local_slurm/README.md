# Local Slurm Configuration

This folder contains configuration files optimized for running the DeSurv pipeline on a local desktop with Slurm (20 CPUs, 31GB RAM).

## Quick Start

### Option 1: Easy/Quick Runs (Functionality Testing)

For quick functionality tests with minimal compute time:

1. **Copy or symlink the config file:**
   ```bash
   # Either copy to use local configs
   cp local_slurm/targets_configs.R targets_configs.R

   # Or use both by sourcing local after main (configs merge by name)
   ```

2. **Submit the job:**
   ```bash
   sbatch local_slurm/_targets.sh
   ```

### Option 2: Full Runs with Local Resources

For production runs using local Slurm resources:

1. **Use the local targets_setup.R:**
   - Copy `local_slurm/targets_setup.R` to overwrite the main `targets_setup.R`, OR
   - Modify `submit_targets.R` to source `local_slurm/targets_setup.R` instead

2. **Submit:**
   ```bash
   sbatch local_slurm/_targets.sh
   ```

## Files in this Folder

| File | Description |
|------|-------------|
| `_targets.sh` | Batch script for main pipeline (16GB mem, no module load) |
| `_targets_bladder.sh` | Batch script for bladder cancer pipeline |
| `_targets_sims.sh` | Batch script for simulation pipeline |
| `_targets_quick.sh` | Batch script for quick/easy mode testing |
| `submit_targets_quick.R` | R script that swaps configs, runs, then restores |
| `targets_setup.R` | Slurm controllers configured for local resources |
| `_targets_sims.R` | Simulation config with decoupled init/CPU counts |
| `targets_configs.R` | Easy configs with reduced initializations |
| `paper_targets.yaml` | Local store path for paper rendering |

## Key Differences from HPC Versions

### Resource Limits

| Parameter | HPC | Local |
|-----------|-----|-------|
| Workers | 202 | 20 |
| CPUs per task | 30-100 | 19 |
| Memory | Varies | 16GB max |
| Module load | Required | Not needed |

### Easy Config Parameters

The `targets_configs.R` in this folder uses "easy" settings for quick runs:

| Parameter | Easy | Full | Effect |
|-----------|------|------|--------|
| `ninit` | 2 | 30 | Random inits during BO |
| `bo_n_init` | 4 | 20 | Initial BO samples |
| `bo_n_iter` | 4 | 50 | BO iterations |
| `bo_candidate_pool` | 200 | 4000 | Candidates per step |
| `bo_max_refinements` | 0 | 1-2 | Refinement rounds |
| `ninit_full` | 10 | 100 | Final model inits |
| `desurv_ncores_grid` | 19 | 30 | Parallel cores |

### Simulation Pipeline Note

In `_targets_sims.R`, `SIM_CV_NSTARTS` (30) is **decoupled** from `cpus_per_task` (19):
- `SIM_CV_NSTARTS=30`: Number of algorithm initializations (statistical parameter)
- `cpus_per_task=19`: Slurm resource allocation (infrastructure parameter)

The parallel code schedules 30 initializations across 19 available CPUs.

## Switching Between Configs

To switch from easy to full configs while staying local:

1. Edit `local_slurm/targets_configs.R`
2. Change parameter values (e.g., `ninit = 30`, `bo_n_iter = 50`)
3. Keep `desurv_ncores_grid = 19` to respect local CPU limits

Or source the original `targets_configs.R` and only override `targets_setup.R`.

## Paper Rendering with Local Store

To render the paper using a local store path (avoiding conflicts with collaborators):

```bash
# Copy the local yaml to paper/ before rendering
cp local_slurm/paper_targets.yaml paper/_targets.yaml

# Render the paper
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'

# Restore original yaml after (or add paper/_targets.yaml to .gitignore)
git checkout paper/_targets.yaml
```
