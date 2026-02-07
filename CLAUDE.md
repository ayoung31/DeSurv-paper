# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Reproducible research repository for **DeSurv**: a survival-driven deconvolution tool for discovering prognostic cancer subtypes. Contains targets pipelines to reproduce all analyses, simulations, and manuscript figures from raw data.

**See [CHANGELOG.md](CHANGELOG.md) for recent changes and development history.**

**See [CONSISTENCY_STANDARDS.md](CONSISTENCY_STANDARDS.md) for code-documentation-paper alignment guidelines.**

## Essential Commands

```bash
# Install DeSurv package (required before running pipelines)
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'

# Run all tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test_bo_helpers.R")'

# Pre-submission check (REQUIRED before sbatch)
./scripts/preflight_check.sh

# Submit pipelines to Slurm
sbatch _targets_wait.sh         # Main analysis (tcgacptac + bladder via config)
sbatch _targets_sims_local.sh   # Simulation studies

# Run specific pipeline directly
Rscript -e 'targets::tar_make(script = "_targets.R")'

# Rebuild simulation figures from pre-computed results (see below)
Rscript -e 'source("sim_figs.R"); sim_results_table <- readRDS("store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/objects/sim_results_table"); save_sim_figs_by_scenario(build_sim_figs_by_scenario(sim_results_table), sim_dir="figures/sim", figure_configs=list())'

# Render manuscript
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'
```

**Note:** Bladder analysis runs via the same `_targets.R` pipeline with bladder configs in `targets_bo_configs.R`. There is no separate `_targets_bladder.R` file.

## Architecture

### Pipeline System

Uses `targets` R package for declarative workflow management with Slurm distribution via `crew.cluster`.

**Pipeline files:**
- `_targets.R` - Main TCGA/CPTAC pancreatic cancer analysis (bladder also runs via config here)
- `_targets_sims.R` - Method validation simulations

**Configuration system:** All hyperparameters flow through separate config files (symlinks to `local_slurm/`):
- `targets_bo_configs.R` → `targets_bo_configs()` - Bayesian optimization bounds and settings
- `targets_run_configs.R` → `targets_run_configs()` - Full-model run parameters (references BO via `bo_key`)
- `targets_val_configs.R` → `targets_val_configs()` - External validation datasets and modes
- `targets_figure_configs.R` - Figure generation settings

Each config gets a hash-based `config_id` and readable `path_tag`. Validation runs at pipeline startup via `validate_desurv_configs()`.

**IMPORTANT:** There is NO `targets_configs.R` (monolithic file). The config files are split into separate modules. Old scripts referencing `targets_configs.R` are obsolete.

### Three-Phase Optimization Flow

1. **Bayesian Optimization** - Hyperparameter search (k, alpha, lambda, nu, ngene, ntop, penalties)
2. **Full-Model Run** - Best params with larger initialization count
3. **External Validation** - Evaluation on held-out datasets

### Key Modules

| File | Purpose |
|------|---------|
| `R/targets_config.R` | Config hashing, validation, diff utilities |
| `R/bo_helpers.R` | BO k-selection logic, bound management |
| `R/figure_targets.R` | All manuscript figures |
| `R/cluster_alignment.R` | Align discovered subtypes to references |
| `R/enrichment_map.R` | Pathway enrichment analysis |
| `targets_setup.R` | Library loading, Slurm controller definitions |
| `targets_common_pipeline.R` | Shared target definitions |

### Crew Resource Controllers

Defined in `targets_setup.R`. Five controllers (`cv`, `default`, `med_mem`, `full`, `low_mem`) tuned for 20 CPUs / 30GB RAM. See `targets_setup.R` for worker counts, memory budgets, and CPU constraints.

**IMPORTANT — avoid OOM crashes:**
- Changing `desurv_parallel_grid` or `desurv_ncores_grid` in `targets_bo_configs.R` **invalidates all BO targets** (they are inside `bo_config`, which is hashed). Only change these if you're willing to rerun BO from scratch.
- Worker counts and `seconds_idle` in `targets_setup.R` do NOT invalidate targets — safe to adjust anytime.

### Data Modes

- `data_mode = "external"` - Separate train/validation files
- `data_mode = "split"` - Train/test split from single dataset
- `mode = "external"` / `mode = "train_split"` - Validation source

## DeSurv Package Dependency

The pipeline requires the DeSurv package to be **installed** (not just checked out). The pipeline uses `library(DeSurv)`, not `pkgload::load_all()`.

```bash
# Install before running pipelines
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'
```

**Note:** If you modify the DeSurv package, you must reinstall it for changes to take effect in the pipeline.

**Core DeSurv functions used:**
- `desurv_data()` - Validate/encapsulate expression matrix + survival
- `desurv_cv_bayesopt_refine()` - BO with iterative bound refinement
- `desurv_fit()` - Main model fitting
- `predict.desurv_fit()` - Risk score predictions
- `desurv_get_top_genes()` - Factor interpretation

## Configuration Patterns

**Fixed vs tuned parameters:**
```r
# Single value = fixed
ngene_config = c(2000)

# Range = BO will tune
ngene_config = c(2000, 5000)
```

The helper `maybe_add_numeric_bound()` handles this logic.

**Compare configs:**
```r
desurv_config_diff(old_config, new_config, type = "bo")
```

## Paper & Manuscript

### File Structure

```
paper/
├── paper.Rmd                 # Main document (combines child sections)
├── 02_introduction_30102025.Rmd  # Introduction
├── 03_methods.Rmd            # Methods section
├── 04_results.Rmd            # Results section (loads most figures)
├── 05_discussion.Rmd         # Discussion section
├── supplement.Rmd            # Supplementary document
├── supp_methods.Rmd          # Supplementary methods (detailed algorithms)
├── references_30102025.bib   # Bibliography
├── pnas.csl                  # Citation style
├── pnas-new.cls              # PNAS LaTeX class
├── paper.pdf                 # Rendered output
└── supp_methods.pdf          # Rendered supplement
```

### Dependency Chain

Pipeline targets → `tar_load()` in `.Rmd` files → `figures/panels/*.pdf` → `paper/paper.Rmd` (child includes) → `paper/paper.pdf`.

### Key Targets Loaded in Results

| Target | Figure | Content |
|--------|--------|---------|
| `sim_figs_by_scenario` | Fig 2E, 3 | Simulation results (loaded via `readRDS`, not `tar_load`) |
| `fig_bo_heat_tcgacptac` | Fig 2D | BO heatmap |
| `fig_residuals_tcgacptac` | Fig 2A | NMF residuals |
| `fig_cophenetic_tcgacptac` | Fig 2B | Cophenetic correlation |
| `fig_silhouette_tcgacptac` | Fig 2C | Silhouette width |
| `fig_gene_overlap_heatmap_*` | Fig 4A-B | Gene program correlations |
| `fig_variation_explained_*` | Fig 4C, 6A | Variance vs survival |
| `fig_hr_forest_tcgacptac` | Fig 5A | Forest plot |
| `fig_median_survival_*` | Fig 5B-C | KM curves |
| `tar_fit_desurv_tcgacptac` | Fig 6B | Fitted model for projection |

### Figure Generation

Figure targets are defined in `R/figure_targets.R` with this pattern:
- `build_fig_*_panels()` - Creates individual panel plots
- `combine_fig_*_panels()` - Arranges panels into composite figure
- `save_fig_*()` - Saves to `figures/` and `figures/panels/`

### Simulation Figures (HPC Store Workaround)

Simulation analyses (`sim_analysis_result`, 2400 branches) are run on HPC by the student. The aggregated `sim_results_table` and `sim_figs_by_scenario` objects are downloaded into the local store's `objects/` directory. However, the store **metadata** (`meta/`) is not transferred, so `targets` does not recognize these objects as up-to-date.

**DO NOT run `tar_make(names = "sim_figs_by_scenario", script = "_targets_sims.R")` locally** — it will attempt to recompute all 2,400 analysis branches from scratch and appear to freeze.

Instead, rebuild simulation figures directly from the pre-computed results table:

```r
source("sim_figs.R")
store <- "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full"
sim_results_table <- readRDS(file.path(store, "objects", "sim_results_table"))
sim_figs_by_scenario <- build_sim_figs_by_scenario(sim_results_table)

# Save rebuilt object back to store (REQUIRED for paper rendering)
saveRDS(sim_figs_by_scenario, file.path(store, "objects", "sim_figs_by_scenario"))

# Also save standalone PDFs
save_sim_figs_by_scenario(sim_figs_by_scenario,
                          sim_dir = "figures/sim",
                          figure_configs = list())
```

**IMPORTANT:** You must save the rebuilt object back to the store. The paper (`04_results.Rmd:25`) loads `sim_figs_by_scenario` via `readRDS()` from the store. The HPC-built ggplot objects are **not renderable locally** due to ggplot2 version mismatch (`ggplot_build` fails on deserialized objects from different versions). Always rebuild from `sim_results_table` (raw data, version-agnostic) rather than reusing serialized ggplot objects.

### Rendering Commands

```bash
# Render main paper (requires targets store)
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'

# Render supplementary methods
Rscript -e 'rmarkdown::render("paper/supp_methods.Rmd")'

# The paper uses this targets store (set in paper/_targets.yaml):
# store: store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full
```

### Targets Store for Paper

The paper reads from a specific targets store. **All four store references must match:**

| File | Setting |
|------|---------|
| `_targets.yaml` (root) | `main: store:` |
| `paper/_targets.yaml` | `main: store:` |
| `paper/paper.Rmd` | `params: tar_store:` |
| `_targets_sims_local.sh` | `tar_config_set(store = ...)` |

**Current store:** `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`

To render with a different store:
```r
# In paper.Rmd, the store is set via params$tar_store
rmarkdown::render("paper/paper.Rmd",
  params = list(tar_store = "your_store_name"))
```

**Run `scripts/verify_consistency.R` to check all store paths match.**

## Code Style

- 2-space indentation, snake_case for functions/files
- Configuration edits go in `targets_*_configs.R` files, not hardcoded
- Keep helpers small and composable
- Commit messages: short, imperative, lowercase
- Shell scripting patterns: see [docs/shell-best-practices.md](docs/shell-best-practices.md)

## Consistency Checking

Before committing changes, run `Rscript scripts/verify_consistency.R`. Key docs: [CONSISTENCY_STANDARDS.md](CONSISTENCY_STANDARDS.md), [DEFAULTS.md](DEFAULTS.md), [CHANGELOG.md](CHANGELOG.md).

## Optional Dependencies

ORA/KEGG enrichment requires `clusterProfiler` and `org.Hs.eg.db`.

## Local Development & Config Architecture

### Config Symlink System

The root config files are **symlinks** to files in `local_slurm/`:
```
targets_bo_configs.R     → local_slurm/targets_bo_configs.R
targets_run_configs.R    → local_slurm/targets_run_configs.R
targets_val_configs.R    → local_slurm/targets_val_configs.R
targets_figure_configs.R → local_slurm/targets_figure_configs.R
```

The `local_slurm/` directory contains:
- Root level: Current active configs (symlink targets)
- `quick/`: Reduced iteration configs for testing
- `full/`: Publication-quality configs

**WARNING:** The `local_slurm/full/` configs may diverge from root configs. Always verify which config is active before running pipelines.

### Obsolete Files (Do Not Use)

| File | Status |
|------|--------|
| `local_slurm/targets_configs.R` | Obsolete monolithic file, never sourced |
| `local_slurm/submit_targets_quick.R` | References old file structure |
| `local_slurm/_targets_sims.R` | Incomplete stub, needs main file |

### Switching Configs

```bash
# Check current mode
cat local_slurm/.current_mode

# To use quick mode (for testing only):
# 1. Backup current symlinks
# 2. Point symlinks to local_slurm/quick/
# 3. Run pipeline
# 4. Restore symlinks
```

## Current Branch Configuration

**Branch:** `naimedits0125` (or check with `git branch --show-current`)
**Store:** `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`

### Local Desktop Mode

The pipeline is configured for local desktop execution (bypassing HPC Slurm issues):
- Main pipeline uses `crew_controller_local` with 2 workers
- Simulation pipeline uses `crew_controller_local` with 2 workers
- Resource limits: 19 CPUs, ~30GB RAM available

### Key Configuration Files

| File | Purpose |
|------|---------|
| `targets_bo_configs.R` | BO hyperparameter bounds (k, alpha, lambda, nu, ngene, ntop) |
| `targets_setup.R` | Controller definitions, package loading |
| `paper/_targets.yaml` | Store path for paper rendering |
| `paper/paper.Rmd` | Paper params including `tar_store` |

For monitoring, diagnostics, pitfall avoidance, and Slurm lifecycle: see [docs/pipeline-operations.md](docs/pipeline-operations.md).

## Open Issues

### Pipeline

- **`targets_setup.R:110`** - `DEFAULT_NINIT` referenced before definition

### DeSurv Package

- **`R/cv_helpers.R:231`** - Missing gene subset validation in CV
- **`R/predict_methods.R:71-106`** - Non-finite validation inconsistency

See [CODE_REVIEW.md](CODE_REVIEW.md) for full findings and resolved items.

## Reference Documents

| Document | Content |
|----------|---------|
| [CODE_REVIEW.md](CODE_REVIEW.md) | Bug findings, architecture issues, open items |
| [NARRATIVE_ARC.md](NARRATIVE_ARC.md) | Paper storytelling structure analysis |
| [PNAS_REVIEW.md](PNAS_REVIEW.md) | Publication readiness, compliance checklist |
| [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md) | Panel-by-panel visualization review |
| [SUGGESTED_TEXT.md](SUGGESTED_TEXT.md) | Drop-in replacement language for claims |
| [MANUSCRIPT_PIPELINE_ANALYSIS.md](MANUSCRIPT_PIPELINE_ANALYSIS.md) | Target-to-figure mapping, config impact |
| [docs/pipeline-operations.md](docs/pipeline-operations.md) | Monitoring, diagnostics, pitfalls, Slurm lifecycle |
| [docs/shell-best-practices.md](docs/shell-best-practices.md) | Shell scripting patterns for this project |
