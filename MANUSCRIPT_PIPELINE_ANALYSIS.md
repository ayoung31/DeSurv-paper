# DeSurv Manuscript Pipeline Analysis

This document provides a comprehensive analysis of the DeSurv paper's computational pipeline, tracing all dependencies from the manuscript PDF back through targets, figures, and source data.

---

## Table of Contents

1. [Manuscript Generation Overview](#1-manuscript-generation-overview)
2. [Complete Dependency Chain](#2-complete-dependency-chain)
3. [Figure-to-Target Mapping](#3-figure-to-target-mapping)
4. [Configuration System](#4-configuration-system)
5. [Simulation Pipeline](#5-simulation-pipeline)
6. [Real Data Pipeline](#6-real-data-pipeline)
7. [Impact of Configuration Choices on Claims](#7-impact-of-configuration-choices-on-claims)
8. [Quick vs Full Mode Differences](#8-quick-vs-full-mode-differences)
9. [Key File Locations](#9-key-file-locations)

---

## 1. Manuscript Generation Overview

### Document Structure

The manuscript is rendered from `paper/paper.Rmd` which includes four child documents:

```
paper/paper.Rmd (main)
├── paper/02_introduction_30102025.Rmd (no data dependencies)
├── paper/04_results.Rmd (21 target dependencies)
├── paper/05_discussion.Rmd (no data dependencies)
└── paper/03_methods.Rmd (no data dependencies)
```

**Supplementary documents** (rendered separately):
- `paper/supplement.Rmd` - Additional analyses
- `paper/supp_methods.Rmd` - Detailed algorithms

### Targets Store Configuration

The paper reads from a targets store specified in `paper/_targets.yaml`:
```yaml
main:
  store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125
```

This can be overridden via render parameters:
```r
rmarkdown::render("paper/paper.Rmd",
  params = list(tar_store = "your_store_name"))
```

---

## 2. Complete Dependency Chain

### Targets Loaded by Manuscript

**paper/04_results.Rmd** (primary data consumer):

| Target | Purpose | Figure |
|--------|---------|--------|
| `sim_figs_by_scenario` | Simulation results | Fig 2E, 3 |
| `fig_bo_heat_tcgacptac` | BO heatmap | Fig 2D |
| `fig_residuals_tcgacptac` | NMF residuals | Fig 2A |
| `fig_cophenetic_tcgacptac` | Cophenetic correlation | Fig 2B |
| `fig_silhouette_tcgacptac` | Silhouette width | Fig 2C |
| `fig_gene_overlap_heatmap_desurv_tcgacptac` | DeSurv gene overlap | Fig 4A |
| `fig_gene_overlap_heatmap_std_desurvk_tcgacptac` | NMF gene overlap | Fig 4B |
| `fig_variation_explained_tcgacptac` | Variance vs survival | Fig 4C |
| `fig_desurv_std_correlation_tcgacptac` | Factor correspondence | Fig 4D |
| `fig_hr_forest_tcgacptac` | HR forest plot | Fig 5A |
| `fig_median_survival_desurv_tcgacptac` | DeSurv KM curves | Fig 5B |
| `fig_median_survival_std_desurvk_tcgacptac` | NMF KM curves | Fig 5C |
| `fig_extval_panel_a_tcgacptac` | Validation forest | Fig 6A |
| `fig_extval_panel_b_tcgacptac` | Validation DeSurv KM | Fig 6B |
| `fig_extval_panel_c_tcgacptac` | Validation NMF KM | Fig 6C |
| `fig_variation_explained_bladder` | Bladder variance plot | Fig 7 |
| `tar_fit_desurv_tcgacptac` | Fitted model object | Projections |
| `tar_data_filtered_tcgacptac` | Processed training data | Analysis |
| `raw_data_bladder` | Bladder cancer data | Fig 7 |

**paper/paper.Rmd**:
- `version_info` - Package/git metadata for reproducibility

**paper/supplement.Rmd**:
- `full_model` - Convergence analysis
- `fig_sc_plot` - Single-cell validation (dynamically resolved)

### Upstream Dependencies

Each figure target depends on model fitting and data processing:

```
Raw Data (data/original/*.rds)
    ↓
tar_data (combined datasets)
    ↓
tar_data_filtered (gene selection, preprocessing)
    ↓
desurv_bo_results (Bayesian optimization)
    ↓
tar_params_best (optimal k, alpha, lambda, nu)
    ↓
tar_fit_desurv (trained model: W, H, beta)
    ↓
fig_* targets (visualizations)
    ↓
paper/04_results.Rmd (tar_load)
    ↓
paper/paper.pdf
```

---

## 3. Figure-to-Target Mapping

### Main Paper Figures

| Figure | Content | Targets Required | Key Claims |
|--------|---------|------------------|------------|
| **Fig 1** | Model schematic | None (static) | Method overview |
| **Fig 2** | Rank selection | `fig_bo_heat`, `fig_residuals`, `fig_cophenetic`, `fig_silhouette`, `sim_figs_by_scenario` | DeSurv resolves rank ambiguity |
| **Fig 3** | Simulation results | `sim_figs_by_scenario` | DeSurv outperforms NMF |
| **Fig 4** | Biological factors | `fig_gene_overlap_heatmap_*`, `fig_variation_explained`, `fig_desurv_std_correlation` | Factors capture known programs |
| **Fig 5** | Training survival | `fig_hr_forest`, `fig_median_survival_*` | Prognostic stratification |
| **Fig 6** | External validation | `fig_extval_panel_*` | Generalization to new cohorts |
| **Fig 7** | Bladder cancer | `fig_variation_explained_bladder` | Cross-cancer applicability |

### Figure Generation Functions

Defined in `R/figure_targets.R` and `targets_common_pipeline.R`:

```r
# Pattern: build_fig_*_panels() → combine_fig_*_panels() → save_fig_*()

# Example for BO figure:
build_fig_bo_panels()      # Creates 4 individual panels
combine_fig_bo_panels()    # Arranges into grid
save_fig_bo()              # Saves to figures/
```

### Figure Output Locations

```
figures/
├── panels/                    # Individual panels
│   ├── fig_bo_cvk_tcgacptac.pdf
│   ├── fig_bo_cvalpha_tcgacptac.pdf
│   └── ...
├── fig_bo_full_tcgacptac.pdf  # Combined figures
├── fig_bio_full_tcgacptac.pdf
└── sim/                       # Simulation figures
    └── ...
```

---

## 4. Configuration System

### Configuration Types

The pipeline uses three interconnected configuration types:

#### 1. BO Configuration (`targets_bo_configs.R`)

Controls data loading and Bayesian optimization:

```r
list(
  # Data settings
  data_mode = "external",           # "external" or "split"
  data_loader = "load_data",
  train_datasets = c("TCGA_PAAD", "CPTAC"),

  # BO bounds
  desurv_bo_bounds = list(
    k_grid = list(lower = 2L, upper = 10L, type = "integer"),
    alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
    lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
    nu_grid = list(lower = 0, upper = 1, type = "continuous")
  ),

  # Fixed parameters
  ngene_config = c(2000),
  ntop_config = c(100, 200),
  lambdaw_config = c(0),
  lambdah_config = c(0),

  # BO settings
  ninit = 30,          # CV random initializations
  bo_n_init = 20,      # Initial BO sample points
  bo_n_iter = 50,      # BO iterations
  bo_candidate_pool = 4000,
  bo_max_refinements = 1,

  # CV settings
  nfold = 5
)
```

#### 2. Run Configuration (`targets_run_configs.R`)

Controls final model training:

```r
list(
  ninit_full = 100,    # Full model initializations
  run_tol = 1e-5,
  run_maxit = 4000,
  std_nmf_k_grid = 2:12,
  coxnet_alpha_grid = seq(0, 1, by = 0.1)
)
```

#### 3. Validation Configuration (`targets_val_configs.R`)

Controls external validation:

```r
list(
  mode = "external",   # "external" or "train_split"
  val_datasets = c("Dijk", "Moffitt_GEO_array", "PACA_AU_array",
                   "PACA_AU_seq", "Puleo_array")
)
```

### Configuration Resolution

Configs are resolved and validated at pipeline startup:

```r
# In _targets.R:
RESOLVED_BO_CONFIGS <- resolve_desurv_bo_configs(BO_CONFIGS_RAW)
RESOLVED_RUN_CONFIGS <- resolve_desurv_run_configs_by_bo(...)
RESOLVED_VAL_CONFIGS <- resolve_desurv_val_configs_by_bo(...)
validate_desurv_configs(BO, RUN, VAL)
```

Each resolved config gets:
- `config_id`: SHA1 hash for cache invalidation
- `path_tag`: Human-readable identifier for output paths
- Merged defaults for missing parameters

---

## 5. Simulation Pipeline

### Scenarios

| Scenario ID | Description | k | Beta | Purpose |
|-------------|-------------|---|------|---------|
| `R0_easy` | Default test | 3 | [2.0, 0, 0] | Main comparison |
| `R0k6` | Higher complexity | 6 | [2.0, -1.5, 1.0, 0, 0, 0] | Scalability test |
| `R00_null` | No survival signal | 3 | [0, 0, 0] | Specificity test |
| `R_mixed` | Mixed markers | 3 | [2.0, 0, 0] | Robustness test |

Each scenario: 100 replicates, N=200 samples, G=3000 genes

### Analysis Methods

| Method ID | Alpha | Description | Used For |
|-----------|-------|-------------|----------|
| `bo` | Tuned (0-0.95) | DeSurv with BO | **Main method** |
| `bo_alpha0` | Fixed at 0 | Standard NMF with BO | **NMF baseline** |
| `bo_tune_ntop` | Tuned | DeSurv + ntop tuning | Full pipeline |
| `bo_tune_ntop_alpha0` | Fixed at 0 | NMF + ntop tuning | Full NMF baseline |
| `fixed` | 0.6 | Fixed DeSurv params | Debugging |
| `fixed_alpha0` | 0 | Fixed NMF params | Debugging |

### Metrics Computed

- **C-index**: Survival prediction accuracy (test set)
- **k selection**: Recovered rank vs true rank
- **Precision**: Gene signature recovery (learned vs true markers)
- **Recall/F1**: Additional signature metrics

### Paper Claims from Simulations

**Figure 2E**: "DeSurv more consistently recovers the true rank"
- Data: `sim_figs_by_scenario$R0_easy$k_hist`
- Comparison: `bo` vs `bo_alpha0` k distributions

**Figure 3A**: "Higher survival prediction performance"
- Data: `sim_figs_by_scenario$R0_easy$cindex_box`
- Metric: Test C-index comparison

**Figure 3B**: "Substantially higher precision for gene recovery"
- Data: `sim_figs_by_scenario$R0_easy$precision_box`
- Key insight: NMF shows near-zero precision

---

## 6. Real Data Pipeline

### Datasets

**Training (PDAC)**:
- TCGA_PAAD: The Cancer Genome Atlas
- CPTAC: Clinical Proteomic Tumor Analysis Consortium

**External Validation (PDAC)**:
- Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, Puleo_array

**Bladder Cancer** (separate analysis):
- IMVigor210: 70/30 train/validation split

### Data Flow

```
1. Load raw data (load_data_internal)
   - Expression matrix (.rds)
   - Survival data (.survival_data.rds)
   - Subtype annotations (_subtype.csv)

2. Combine datasets (load_data)
   - Union of genes across datasets
   - Aligned expression matrices
   - Combined sample metadata

3. Preprocess (preprocess_training_data)
   - Gene selection (ngene = 2000)
   - Optional transformation (method_trans_train)

4. Bayesian optimization (desurv_cv_bayesopt_refine)
   - 5-fold CV
   - Gaussian process surrogate
   - k, alpha, lambda, nu tuned

5. Full model training (desurv_fit)
   - 100 random initializations
   - Consensus initialization option
   - Final W, H, beta matrices

6. External validation
   - Project onto learned factors
   - Compute risk scores
   - Calculate C-index per cohort
```

### Model Outputs

| Output | Contents | Usage |
|--------|----------|-------|
| `tar_fit_desurv` | W, H, beta, convergence | Factor interpretation, predictions |
| `tar_tops_desurv` | Top genes per factor | Biological analysis |
| `fit_std_desurvk` | NMF at DeSurv's k | Comparison baseline |
| `fit_std_elbowk` | NMF at elbow-selected k | Alternative baseline |
| `desurv_bo_history` | All BO evaluations | Figure 2D heatmap |

---

## 7. Impact of Configuration Choices on Claims

### Critical Parameters and Their Effects

#### 1. `ninit` (CV Initializations)

**Affects**: Reliability of BO parameter selection

| Value | Effect on Claims |
|-------|------------------|
| 4 (quick) | High variance in selected parameters; results may not be reproducible |
| 30 (full) | Stable parameter estimates; supports "robust optimization" claim |

**Paper claim at risk if too low**: "Bayesian optimization identifies optimal hyperparameters"

#### 2. `bo_n_iter` (BO Iterations)

**Affects**: Quality of hyperparameter search

| Value | Effect on Claims |
|-------|------------------|
| 4 (quick) | May not find global optimum; suboptimal k, alpha |
| 50 (full) | Thorough exploration; supports "comprehensive search" |

**Paper claim at risk if too low**: "Selected k=3 is optimal for PDAC"

#### 3. `ninit_full` (Final Model Initializations)

**Affects**: Model stability and reproducibility

| Value | Effect on Claims |
|-------|------------------|
| 19 (quick) | Factors may vary across runs |
| 100 (full) | Consensus factors; supports "reproducible subtypes" |

**Paper claim at risk if too low**: "Discovered factors correspond to known PDAC programs"

#### 4. `alpha_grid` Range

**Affects**: DeSurv vs NMF comparison validity

| Setting | Effect |
|---------|--------|
| `lower = 0, upper = 0` | Forces alpha=0 (NMF baseline) |
| `lower = 0, upper = 0.95` | Full DeSurv optimization |

**Paper claim requires**: Both `bo` (alpha tuned) AND `bo_alpha0` (alpha=0) runs

#### 5. Simulation Replicates

**Affects**: Statistical power of method comparison

| Replicates | Effect |
|------------|--------|
| 10 | Wide confidence intervals; differences may not be significant |
| 100 | Narrow CIs; robust statistical claims |

**Paper claim at risk if too few**: "DeSurv significantly outperforms NMF (p < 0.001)"

### Claims Requiring Specific Configurations

| Claim | Required Configuration |
|-------|------------------------|
| "DeSurv outperforms NMF" | Both `bo` AND `bo_alpha0` analysis methods |
| "Optimal k=3 for PDAC" | Sufficient bo_n_iter (≥50) |
| "Factors match known programs" | High ninit_full (≥100) for stable factors |
| "Generalizes to external cohorts" | External validation mode with multiple datasets |
| "Applicable to bladder cancer" | Bladder config with split validation |

---

## 8. Quick vs Full Mode Differences

### Parameter Comparison

| Parameter | Quick | Full | HPC Production |
|-----------|-------|------|----------------|
| `ninit` | 4 | 30 | 50 |
| `bo_n_init` | 4 | 20 | 20 |
| `bo_n_iter` | 4 | 50 | 60 |
| `bo_candidate_pool` | 200 | 4000 | 4000 |
| `bo_max_refinements` | 0 | 1 | 2 |
| `ninit_full` | 19 | 100 | 100 |
| `desurv_ncores_grid` | 1 | 4 | 19+ |

### What Quick Mode Affects

1. **Figure quality**: Results may differ from published figures
2. **Statistical claims**: Reduced replicates weaken significance
3. **Missing comparisons**: May skip `bo_alpha0` runs entirely
4. **Reproducibility**: Different seeds/inits → different factors

### Switching to Full Mode

```bash
# Copy full configs
cp local_slurm/full/targets_*.R local_slurm/

# Or for HPC production
cp .config_backup/targets_*.R .
```

---

## 9. Key File Locations

### Pipeline Files

| File | Purpose |
|------|---------|
| `_targets.R` | Main PDAC pipeline |
| `_targets_sims.R` | Simulation pipeline |
| `_targets_bladder.R` | Bladder cancer pipeline |
| `targets_common_pipeline.R` | Shared target definitions |
| `targets_setup.R` | Slurm/crew configuration |

### Configuration Files

| File | Purpose |
|------|---------|
| `targets_configs.R` | Combined configs (symlink) |
| `local_slurm/targets_bo_configs.R` | BO parameters |
| `local_slurm/targets_run_configs.R` | Run parameters |
| `local_slurm/targets_val_configs.R` | Validation parameters |
| `local_slurm/targets_figure_configs.R` | Figure settings |

### R Functions

| File | Purpose |
|------|---------|
| `R/targets_config.R` | Config validation/resolution |
| `R/figure_targets.R` | Figure generation |
| `R/bo_helpers.R` | BO utilities |
| `R/cluster_alignment.R` | Subtype matching |
| `sim_figs.R` | Simulation figures |

### Data Files

| Path | Contents |
|------|----------|
| `data/original/*.rds` | Raw expression matrices |
| `data/original/*.survival_data.rds` | Survival annotations |
| `data/derv/` | Processed data |
| `results/` | Model outputs |

### Manuscript Files

| File | Purpose |
|------|---------|
| `paper/paper.Rmd` | Main document |
| `paper/04_results.Rmd` | Results (loads targets) |
| `paper/_targets.yaml` | Store configuration |
| `figures/` | Generated figures |

---

## Summary: Ensuring Valid Paper Claims

To support all paper claims, ensure:

1. **Both DeSurv and NMF baselines run**: Analysis specs include `bo` AND `bo_alpha0`
2. **Sufficient BO iterations**: `bo_n_iter ≥ 50` for thorough search
3. **Adequate initializations**: `ninit_full ≥ 100` for stable factors
4. **All simulation scenarios**: R0_easy, R00_null, R0k6, R_mixed
5. **Full replicates**: 100 per scenario for statistical power
6. **External validation**: All PDAC cohorts + bladder cancer

**Quick mode limitations**:
- Missing NMF comparisons (only 1 box in boxplots)
- Reduced statistical power
- Potentially different k selection
- Figures may not match published results

---

*Generated: 2026-01-25*
*Repository: DeSurv-paper*
*Branch: naimedits0125*
