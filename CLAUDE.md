# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains code to reproduce findings in the DeSurv paper: "A survival-driven deconvolution tool to prognostic subtype discovery". It uses the R `targets` pipeline framework to run DeSurv from original data to final figures.

## Dependencies

The pipeline depends on the `DeSurv` package located at `../DeSurv`. Install before running:
```r
devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)
```

ORA/KEGG enrichment requires `clusterProfiler` and `org.Hs.eg.db` packages.

## Common Commands

### Running Pipelines

Run on Slurm cluster (via batch script):
```bash
sbatch _targets.sh           # Main TCGA/CPTAC pipeline
sbatch _targets_bladder.sh   # Bladder cancer pipeline
sbatch _targets_sims.R       # Simulation pipeline
```

**Note:** These pipelines require a Slurm cluster environment. The targets use Slurm-based crew controllers for distributed execution.

Specify a different pipeline file:
```bash
Rscript -e 'targets::tar_make(script = "_targets_bladder.R")'
```

### Testing
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Pipeline Inspection
```r
targets::tar_manifest()       # View all targets
targets::tar_visnetwork()     # Visualize dependency graph
targets::tar_read(target_name) # Read a target's cached value
targets::tar_invalidate(target_name) # Force target to rerun
```

## Architecture

### Pipeline Files

- `_targets.R` - Main pipeline for TCGA/CPTAC pancreatic cancer data
- `_targets_bladder.R` - Pipeline for IMVigor210 bladder cancer data
- `_targets_sims.R` - Simulation pipeline with scenario-based data generation
- `targets_setup.R` - Shared configuration: Slurm controllers, package loading, global parameters
- `targets_common_pipeline.R` - Shared target definitions (`COMMON_DESURV_TARGETS` list) used across pipelines

### Pipeline Workflow Stages

The common pipeline (`targets_common_pipeline.R`) implements this workflow:

1. **Bayesian Optimization** - `desurv_bo_results` runs coarse BO then iterative refinement via `desurv_cv_bayesopt_refine()`
2. **Data Preprocessing** - `data_filtered` applies optimal gene filter from BO
3. **Seed-based Fitting** - `desurv_seed_fits` runs 100 fits with different random seeds to evaluate initialization robustness
4. **Consensus Initialization** - `desurv_consensus_init` synthesizes best features from seed fits to create W0, H0, beta0
5. **Final Model Fit** - `fit_desurv` uses consensus initialization for the production model
6. **Downstream Analysis** - Top genes extraction, ORA enrichment, clustering on validation datasets

Three model variants are fitted in parallel: DeSurv (full), DeSurv (alpha=0), and standard NMF baseline.

### Slurm Controllers

Defined in `targets_setup.R`. Select via `resources = tar_resources(crew = tar_resources_crew(controller = "name"))`:

| Controller | Workers | Memory | CPUs | Time | Use Case |
|------------|---------|--------|------|------|----------|
| `low_mem` | 202 | 1GB/cpu | 1 | 120min | Light tasks, data loading |
| `cv` | 202 | 2GB/cpu | NINIT | 720min | Cross-validation, BO |
| `full` | 202 | 32GB total | NINIT_FULL | 600min | Final model fitting |
| `med_mem` | 120 | 8GB total | 1 | 200min | Medium memory tasks |

### Store Path Convention

Stores are named `store_PKG_VERSION={version}_GIT_BRANCH={branch}` to enable parallel runs of different package versions or git branches. The path is set in `submit_targets*.R`.

### Key Configuration Constants (in `_targets*.R`)

- `NGENE_CONFIG` - Gene filter size; single value fixes it, range (e.g., `c(2000, 5000)`) enables BO tuning
- `NTOP_CONFIG` - Genes per factor for downstream analysis; same tuning behavior
- `LAMBDAW_CONFIG`/`LAMBDAH_CONFIG` - Factor penalty terms; singleton fixes, range enables BO tuning
- `DESURV_BO_BOUNDS` - Bayesian optimization search space for k, alpha, lambda, nu

The helper `maybe_add_numeric_bound()` in `R/bo_helpers.R` handles conditional BO bound construction based on whether CONFIG values are fixed or ranges.

### R/ Directory

Helper functions sourced by `targets_setup.R`. Key modules:
- `bo_helpers.R` - `maybe_add_numeric_bound()` for dynamic BO bound construction
- `load_data.R` / `load_data_internal.R` - Dataset loading with multi-dataset gene union support
- `ora_analysis.R` - ORA/KEGG pathway enrichment via clusterProfiler
- `run_clustering*.R` - Consensus clustering workflows
- `simulation_functions/` - Scenario generation (`simulate_*.R`), precision/recall metrics

### Data Flow

1. Raw data in `data/original/` (RDS files with expression + survival)
2. Preprocessing via `DeSurv::preprocess_data()`
3. Hyperparameter tuning via `DeSurv::desurv_cv_bayesopt_refine()`
4. Final model fitting via `DeSurv::desurv_fit()` with consensus initialization
5. Downstream: ORA enrichment, clustering, paper rendering

### Simulation Pipeline

`_targets_sims.R` generates synthetic data to benchmark DeSurv:
- `SIM_DATASETS_PER_SCENARIO` controls replicates per scenario (default: 100)
- Scenarios are defined in `SIMULATION_SCENARIOS` list (e.g., R0_easy, R2_correlated)
- Analysis modes: "fixed" (predetermined params), "fixed_alpha0", "bayesopt" (full BO)
- Uses sequential controller (no Slurm) by default

### Paper Generation

R Markdown files in `paper/` are rendered via `tarchetypes::tar_render()`. The main document is `paper/paper.Rmd` which includes child documents for each section (`02_introduction*.Rmd`, `03_methods.Rmd`, `04_results.Rmd`, `05_discussion.Rmd`).
