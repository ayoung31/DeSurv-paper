# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Reproducible research repository for **DeSurv**: a survival-driven deconvolution tool for discovering prognostic cancer subtypes. Contains targets pipelines to reproduce all analyses, simulations, and manuscript figures from raw data.

## Essential Commands

```bash
# Install DeSurv package (required before running pipelines)
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'

# Run tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Submit pipelines to Slurm
sbatch _targets.sh              # Main TCGA/CPTAC analysis
sbatch _targets_bladder.sh      # Bladder cancer analysis
sbatch _targets_sims.sh         # Simulation studies

# Run specific pipeline directly
Rscript -e 'targets::tar_make(script = "_targets.R")'

# Render manuscript
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'
```

## Architecture

### Pipeline System

Uses `targets` R package for declarative workflow management with Slurm distribution via `crew.cluster`.

**Pipeline files:**
- `_targets.R` - Main TCGA/CPTAC pancreatic cancer analysis
- `_targets_bladder.R` - Bladder cancer analysis
- `_targets_sims.R` - Method validation simulations

**Configuration system:** All hyperparameters flow through `targets_configs.R`:
- `targets_bo_configs()` - Bayesian optimization bounds and settings
- `targets_run_configs()` - Full-model run parameters (references BO via `bo_key`)
- `targets_val_configs()` - External validation datasets and modes

Each config gets a hash-based `config_id` and readable `path_tag`. Validation runs at pipeline startup via `validate_desurv_configs()`.

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

### Slurm Resource Controllers

Defined in `targets_setup.R`:
- `cv_comp_controller` - BO phase (202 workers, 2GB/CPU)
- `full_run_controller` - Full runs (202 workers, 32GB total)
- `low_mem_controller` - Light tasks (1GB/CPU)

### Data Modes

- `data_mode = "external"` - Separate train/validation files
- `data_mode = "split"` - Train/test split from single dataset
- `mode = "external"` / `mode = "train_split"` - Validation source

## DeSurv Package Dependency

The pipeline requires a local checkout of `../DeSurv`. Each `_targets*.R` script auto-loads it via `pkgload::load_all("../DeSurv")`.

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

```
Pipeline Targets (targets store)
       ↓
   tar_load() calls in .Rmd files
       ↓
paper/04_results.Rmd ──────────────────────────────┐
paper/supplement.Rmd                               │
       ↓                                           │
figures/panels/*.pdf (555 panel files)             │
figures/*.pdf (combined figures)                   │
       ↓                                           │
paper/paper.Rmd ← child includes ─────────────────┤
  ├── 02_introduction_30102025.Rmd                 │
  ├── 04_results.Rmd                               │
  ├── 05_discussion.Rmd                            │
  └── 03_methods.Rmd                               │
       ↓                                           │
paper/paper.pdf                                    │
paper/supp_methods.pdf                             │
```

### Key Targets Loaded in Results

The results section (`04_results.Rmd`) loads these targets via `tar_load()`:

| Target | Figure | Content |
|--------|--------|---------|
| `sim_figs_by_scenario` | Fig 2E, 3 | Simulation results |
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

### Rendering Commands

```bash
# Render main paper (requires targets store)
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'

# Render supplementary methods
Rscript -e 'rmarkdown::render("paper/supp_methods.Rmd")'

# The paper uses this targets store (set in paper/_targets.yaml):
# store: store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main
```

### Targets Store for Paper

The paper reads from a specific targets store defined in `paper/_targets.yaml`. To render with a different store:
```r
# In paper.Rmd, the store is set via params$tar_store
rmarkdown::render("paper/paper.Rmd",
  params = list(tar_store = "your_store_name"))
```

## Code Style

- 2-space indentation, snake_case for functions/files
- Configuration edits go in `targets_configs.R`, not hardcoded
- Keep helpers small and composable
- Commit messages: short, imperative, lowercase

## Optional Dependencies

ORA/KEGG enrichment requires `clusterProfiler` and `org.Hs.eg.db`.

## Local Development

For quick testing on reduced resources, use configs from `local_slurm/`:
```bash
cp local_slurm/targets_configs.R targets_configs.R
sbatch local_slurm/_targets.sh
```

## Known Issues & Code Review

See [CODE_REVIEW.md](CODE_REVIEW.md) for comprehensive code review findings including:
- **Critical bugs** requiring immediate attention (e.g., `compute_metrics.R:9` variable assignment)
- **Architectural issues** in pipeline configuration and resource management
- **DeSurv package issues** in C++ numerical stability and validation
- **Test coverage gaps** for critical functions

### Critical Items to Address

1. **`R/compute_metrics.R:9`** - `sdZ` incorrectly assigned from `meanZ`
2. **`R/select_best_init.R:5`** - Filtered dataframe `keep` never used
3. **`targets_common_pipeline.R`** - Remove `browser()` debug statements
4. **`targets_setup.R:110`** - `DEFAULT_NINIT` referenced before definition

### DeSurv Package Critical Items

1. **`src/functions.cpp:472`** - Beta backtracking logic inverted
2. **`R/cv_helpers.R:231`** - Missing gene subset validation in CV
3. **`R/predict_methods.R:71-106`** - Non-finite validation inconsistency
