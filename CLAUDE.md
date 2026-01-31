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

# Run tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Submit pipelines to Slurm
sbatch _targets.sh              # Main TCGA/CPTAC analysis
sbatch _targets_bladder.sh      # Bladder cancer analysis
sbatch _targets_sims.sh         # Simulation studies (HPC mode)
sbatch _targets_sims_local.sh   # Simulation studies (local desktop mode)

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

**Configuration system:** All hyperparameters flow through separate config files:
- `targets_bo_configs.R` → `targets_bo_configs()` - Bayesian optimization bounds and settings
- `targets_run_configs.R` → `targets_run_configs()` - Full-model run parameters (references BO via `bo_key`)
- `targets_val_configs.R` → `targets_val_configs()` - External validation datasets and modes
- `targets_figure_configs.R` - Figure generation settings

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

Defined in `targets_setup.R` (HPC mode - not used in local desktop mode):
- `cv_comp_controller` - BO phase (202 workers, 2GB/CPU)
- `full_run_controller` - Full runs (202 workers, 32GB total)
- `low_mem_controller` - Light tasks (1GB/CPU)

**Note:** Local desktop mode uses `crew_controller_local` instead. See "Current Branch Configuration" section.

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
# store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full
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
- Configuration edits go in `targets_*_configs.R` files, not hardcoded
- Keep helpers small and composable
- Commit messages: short, imperative, lowercase

## Consistency Checking

Before committing changes, run:
```bash
Rscript scripts/verify_consistency.R
```

This checks:
- Store path consistency across config files
- Config file validation
- Pipeline validation
- Debug statement detection
- Figure target status
- Documentation completeness

**Key documents for consistency:**
- [CONSISTENCY_STANDARDS.md](CONSISTENCY_STANDARDS.md) - Full guidelines
- [DEFAULTS.md](DEFAULTS.md) - Parameter defaults registry
- [CHANGELOG.md](CHANGELOG.md) - Change history

## Optional Dependencies

ORA/KEGG enrichment requires `clusterProfiler` and `org.Hs.eg.db`.

## Local Development

For quick testing on reduced resources, use configs from `local_slurm/`:
```bash
cp local_slurm/targets_configs.R targets_configs.R
sbatch local_slurm/_targets.sh
```

## Current Branch Configuration

**Branch:** `naimedits0125`
**Store:** `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full`

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

### Monitoring Jobs

```bash
squeue -u $USER                    # View job queue
tail -f logs/slurm-<JOBID>.out     # Monitor job output
sacct -j <JOBID>                   # Job accounting info
```

## Known Issues & Code Review

See [CODE_REVIEW.md](CODE_REVIEW.md) for comprehensive code review findings including:
- **Critical bugs** - mostly resolved (see status below)
- **Architectural issues** in pipeline configuration and resource management
- **DeSurv package issues** in C++ numerical stability and validation
- **Test coverage gaps** for critical functions

### Critical Items (Status)

1. ~~**`R/compute_metrics.R:9`** - `sdZ` incorrectly assigned from `meanZ`~~ **RESOLVED**: Bug fixed and file removed (was dead code - see BUG1_IMPACT_ANALYSIS.md)
2. ~~**`R/select_best_init.R:5`** - Filtered dataframe `keep` never used~~ **RESOLVED**: File removed (was dead code - never called in pipeline)
3. ~~**`targets_common_pipeline.R`** - Remove `browser()` debug statements~~ **RESOLVED**: All 6 browser() statements removed
4. **`targets_setup.R:110`** - `DEFAULT_NINIT` referenced before definition

### DeSurv Package Critical Items

1. ~~**`src/functions.cpp:472`** - Beta backtracking logic inverted~~ **NOT A BUG**: Analysis confirmed the logic is correct for maximization (see BUG4_ANALYSIS.md)
2. **`R/cv_helpers.R:231`** - Missing gene subset validation in CV
3. **`R/predict_methods.R:71-106`** - Non-finite validation inconsistency

## Narrative Arc Analysis

See [NARRATIVE_ARC.md](NARRATIVE_ARC.md) for detailed analysis of the paper's storytelling structure:
- Five-act story structure (Problem → Solution → Proof → Insight → Generalization)
- Figure-by-figure narrative role analysis
- The central insight: "Variance ≠ Prognosis"
- Rhetorical strategies (problem-solution framing, escalating claims, consistent baseline)
- Information disclosure hierarchy (main text vs supplement)

## PNAS Submission Review

See [PNAS_REVIEW.md](PNAS_REVIEW.md) for comprehensive publication readiness assessment including:
- **PNAS compliance checklist** (word count, figures, references, required sections)
- **Critical code bugs** that may affect reported results
- **Paper-code inconsistencies** in mathematical formulations
- **Mathematical notation issues** (symbol overloading, typos, incomplete equations)
- **Scientific content assessment** (strengths, weaknesses, PNAS suitability)
- **Supplement completeness** review
- **Action checklist** organized by priority

### Key PNAS Issues

1. **Figures:** 6 figures (PNAS max: 4) - move 2 to supplement
2. **References:** 71 references (PNAS max: 50) - reduce by 21+
3. **Placeholders:** Keywords, author contributions, conflict of interest incomplete
4. **Methods section:** Too brief (~30 lines) for PNAS
5. **Discussion section:** Too brief (6 lines) for PNAS

## Figure Analysis

See [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md) for panel-by-panel evaluation of all manuscript figures:
- **Visualization effectiveness** - whether each panel optimally supports its narrative point
- **Proposed improvements** - alternative plotting approaches and annotations
- **Aesthetic assessment** - structure, proportions, color consistency, labels
- **Narrative arc impact** - how figure improvements strengthen the paper's story

### High-Impact Figure Changes (Minimal Effort)

1. **Figure 2:** Reorder panels B-C-D above A (problem before solution)
2. **Figure 4:** ADD FACTOR LABELS to heatmaps (dramatically improves interpretability)
3. **Figures 5-6:** Add HR, CI, p-value to all KM plots (standard practice)
4. **Figure 3:** Add "true k=3, n=100 replicates" annotation
5. **All figures:** Standardize legend terminology ("Standard NMF" vs "DeSurv")

## Suggested Replacement Text

See [SUGGESTED_TEXT.md](SUGGESTED_TEXT.md) for drop-in replacement language implementing editorial recommendations:
- **Revised Introduction** - 5-paragraph replacement with proper prior work acknowledgment
- **Related Work paragraph** - Explicitly acknowledges DECODER and prior PDAC work
- **6 targeted claim fixes** - Specific sentence replacements for W-vs-H, cross-cancer, overfitting claims
- **Revised Significance Statement** - More concrete alternative

### Key Text Changes

1. **Acknowledge Huang et al. and Le Goff et al.** - prior survival-NMF work cited fairly
2. **Acknowledge DECODER** - UNC prior work explicitly recognized
3. **Temper cross-cancer claim** - "consistent with prior reports" not novel discovery
4. **Add W-vs-H gradient explanation** - preempts mathematical reviewers
5. **Add explicit out-of-sample statement** - preempts double-dipping concerns

## Manuscript Pipeline Analysis

See [MANUSCRIPT_PIPELINE_ANALYSIS.md](MANUSCRIPT_PIPELINE_ANALYSIS.md) for comprehensive documentation of:
- **Complete dependency chain** from paper.pdf back to raw data
- **Figure-to-target mapping** showing which targets produce which figures
- **Configuration system** explaining BO, Run, and Validation configs
- **Simulation pipeline** with scenarios, analysis methods, and metrics
- **Real data pipeline** with datasets, preprocessing, and validation
- **Impact of configuration choices on paper claims**
- **Quick vs Full mode differences** and their effects on results

### Key Insights for Reproducibility

1. **Paper claims require both DeSurv AND NMF baselines** - `bo` and `bo_alpha0` analysis methods must both run
2. **Quick mode may produce different results** - fewer initializations, shorter BO, missing comparisons
3. **Figure 3 boxplots need alpha=0 runs** - only shows one box if NMF baseline is missing
4. **External validation requires all cohorts** - Dijk, Moffitt, PACA_AU_array, PACA_AU_seq, Puleo

### Critical Targets for Paper

| Target | Figure | Must Have |
|--------|--------|-----------|
| `sim_figs_by_scenario` | Fig 2E, 3 | Both bo and bo_alpha0 results |
| `fig_bo_heat_tcgacptac` | Fig 2D | Sufficient bo_n_iter |
| `tar_fit_desurv_tcgacptac` | Multiple | Stable factors (high ninit_full) |
| `fig_extval_*` | Fig 6 | All validation datasets |
