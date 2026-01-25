# Local Slurm Configuration

This folder contains configuration files for running the DeSurv pipeline on a local desktop with Slurm (20 CPUs, 31GB RAM).

**Key feature:** These configs don't modify the student's original files - they use symlinks that can be easily switched.

## Quick Start

### 1. Set Up Local Configuration

```bash
# From the project root directory:

# Quick mode (fast testing - default)
./local_slurm/setup_local.sh

# OR full quality mode (HPC-quality iterations)
./local_slurm/setup_local.sh --full
```

This creates symlinks from `local_slurm/` configs to the main directory.

### 2. Install DeSurv Package

```bash
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'
```

### 3. Run the Pipeline

```bash
# Main PDAC analysis
sbatch local_slurm/_targets.sh

# Or for simulations
sbatch local_slurm/_targets_sims.sh

# Or for bladder cancer
sbatch local_slurm/_targets_bladder.sh
```

### 4. Switch Modes or Restore

```bash
# Switch to full quality mode
./local_slurm/setup_local.sh --full

# Switch back to quick mode
./local_slurm/setup_local.sh --quick

# Check current mode
./local_slurm/setup_local.sh --status

# Restore original HPC configs
./local_slurm/setup_local.sh --restore
```

## Two Modes: Quick vs Full

| Parameter | QUICK (testing) | FULL (publication) | HPC (original) |
|-----------|-----------------|-------------------|----------------|
| `ninit` | 4 | 30 | 30 |
| `bo_n_init` | 4 | 20 | 20 |
| `bo_n_iter` | 4 | 50 | 50 |
| `bo_candidate_pool` | 200 | 4000 | 4000 |
| `ninit_full` | 19 | 100 | 100 |
| `desurv_ncores_grid` | 19 | 19 | 30-100 |
| **Purpose** | Fast testing | Publication quality | Full HPC |

### When to Use Each Mode

- **QUICK mode**: Initial testing, debugging, verifying pipeline runs
- **FULL mode**: Final runs for publication-quality results (slower but statistically sound)

## What the Setup Script Does

1. **Backs up** original config files to `.config_backup/`
2. **Creates symlinks** from main directory to `local_slurm/<mode>/` versions
3. **Creates** the `logs/` directory if it doesn't exist
4. **Tracks** current mode in `.current_mode` file

```
targets_setup.R       -> local_slurm/<mode>/targets_setup.R
targets_bo_configs.R  -> local_slurm/<mode>/targets_bo_configs.R
targets_run_configs.R -> local_slurm/<mode>/targets_run_configs.R
targets_val_configs.R -> local_slurm/<mode>/targets_val_configs.R
targets_figure_configs.R -> local_slurm/<mode>/targets_figure_configs.R
```

## Directory Structure

```
local_slurm/
├── setup_local.sh          # Setup script (run this!)
├── README.md               # This file
├── _targets.sh             # Batch script for main pipeline
├── _targets_bladder.sh     # Batch script for bladder cancer
├── _targets_sims.sh        # Batch script for simulations
├── _targets_quick.sh       # Batch script for quick testing
├── quick/                  # Quick testing configs
│   ├── targets_setup.R
│   ├── targets_bo_configs.R
│   ├── targets_run_configs.R
│   ├── targets_val_configs.R
│   └── targets_figure_configs.R
└── full/                   # Full quality configs
    ├── targets_setup.R
    ├── targets_bo_configs.R
    ├── targets_run_configs.R
    ├── targets_val_configs.R
    └── targets_figure_configs.R
```

## Key Differences from HPC Versions

### Resource Limits

| Parameter | UNC Longleaf | Local Desktop |
|-----------|--------------|---------------|
| Workers | 202 | 20 |
| CPUs per task | 30-100 | 19 |
| Memory | Varies | 16GB max |
| `module load` | Required | Not needed |

### Parallelization Note

In **full mode**, `ninit_full=100` means 100 random initializations for the final model. These will be scheduled across the 19 available CPUs (the parallel code handles this automatically). Similarly, `ninit=30` during BO cross-validation will parallelize across 19 CPUs.

## Troubleshooting

### "module: command not found"

This is expected on local systems without module system. The local configs don't use `module load`.

### Jobs not starting

Check Slurm status:
```bash
sinfo           # Cluster status
squeue          # Job queue
```

### Permission denied on setup script

```bash
chmod +x local_slurm/setup_local.sh
```

### Check current configuration status

```bash
./local_slurm/setup_local.sh --status
```

### View job logs

```bash
tail -f logs/_targets.out
tail -f logs/_targets.err
```

## Simulation Pipeline Note

In `_targets_sims.R`, `SIM_CV_NSTARTS` (30) is **decoupled** from `cpus_per_task` (19):
- `SIM_CV_NSTARTS=30`: Number of algorithm initializations (statistical parameter)
- `cpus_per_task=19`: Slurm resource allocation (infrastructure parameter)

The parallel code schedules 30 initializations across 19 available CPUs.

## Paper Rendering with Local Store

To render the paper using a local store path:

```bash
# Copy the local yaml to paper/ before rendering
cp local_slurm/paper_targets.yaml paper/_targets.yaml

# Render the paper
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'

# Restore original yaml after
git checkout paper/_targets.yaml
```

## Verifying the Pipeline Works

After setup, run a quick verification:

```bash
# Check that R can load all required packages
Rscript -e 'library(targets); library(DeSurv); print("OK")'

# Check Slurm is working
sinfo

# Check current mode
./local_slurm/setup_local.sh --status

# Run a minimal test (optional)
Rscript -e 'targets::tar_glimpse()'
```
