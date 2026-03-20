# DeSurv Clean Repository Restructuring Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create clean, reviewer-friendly repositories at `rashidlab/DeSurv` and `rashidlab/DeSurv-paper` for manuscript submission, replacing the `targets` pipeline with numbered R scripts.

**Architecture:** The paper repo ships pre-computed results (exported from the existing targets store) so reviewers can regenerate figures and compile the manuscript without HPC access. A separate numbered-script pipeline allows full re-computation. A quick mode uses DeSurv's built-in reduced-parameter settings to smoke-test the pipeline end-to-end on any laptop.

**Tech Stack:** R, DeSurv (Rcpp/RcppArmadillo), renv, GNU Make, rmarkdown, Zenodo (data hosting)

---

## Dependency Trace: What the Paper Actually Needs

### Traced from paper.Rmd → child .Rmd files → tar_load/tar_read/readRDS/image_read calls

**Child files included by paper.Rmd:**
- `paper/02_introduction_REVISED.Rmd` (text only)
- `paper/04_results_REVISED.Rmd` (main results — heavy dependencies)
- `paper/05_discussion_REVISED.Rmd` (text only)
- `paper/03_methods_REVISED.Rmd` (text only)

**Supplement: `paper/supplement.Rmd`** (additional dependencies)

### Targets store objects referenced by the paper

**From 04_results_REVISED.Rmd (tar_load):**
1. `sim_figs_by_scenario` — simulation figure objects (27 MB)
2. `fig_bo_heat_tcgacptac` — BO heatmap (34 MB)
3. `fig_gene_overlap_heatmap_desurv_tcgacptac` — DeSurv gene program heatmap
4. `fig_gene_overlap_heatmap_std_desurvk_tcgacptac` — NMF gene program heatmap
5. `fig_variation_explained_tcgacptac` — variance vs survival scatter (2.6 MB)
6. `fig_desurv_std_correlation_tcgacptac` — W-matrix correlation
7. `fig_hr_forest_tcgacptac` — forest plot (5.2 MB)
8. `fig_median_survival_desurv_tcgacptac` — DeSurv KM curves (7 MB)
9. `fig_median_survival_std_desurvk_tcgacptac` — NMF KM curves (5.9 MB)
10. `val_latent_desurv_tcgacptac` — validation latent scores (for inline stats)
11. `val_latent_std_desurvk_tcgacptac` — NMF validation latent scores

**From 04_results_REVISED.Rmd (tar_read, inline R):**
12. `tar_k_selection_tcgacptac` — k selection result
13. `desurv_bo_results_tcgacptac` — BO history (34 MB)
14. `tar_params_best_tcgacptac` — best hyperparameters
15. `tar_data_filtered_tcgacptac` — filtered training data (for sample counts)

**From supplement.Rmd (tar_load):**
16. `desurv_seed_fits_tcgacptac` — convergence trajectories (121 MB)
17. `fig_residuals_tcgacptac` — NMF residuals plot (4.7 MB)
18. `fig_cophenetic_tcgacptac` — cophenetic correlation (4.6 MB)
19. `fig_silhouette_tcgacptac` — silhouette width (4.6 MB)
20. `fit_std_tcgacptac` — standard NMF fit (2.2 MB)
21. `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` — NMF k=5 heatmap
22. `val_cindex_desurv_tcgacptac` — DeSurv C-index table
23. `val_cindex_std_desurvk_tcgacptac` — NMF k=3 C-index
24. `val_cindex_std_elbowk_tcgacptac` — NMF k=5 C-index
25. `val_cindex_desurv_alpha0_tcgacptac` — NMF k=7 (BO α=0) C-index
26. `val_latent_desurv_alpha0_tcgacptac` — NMF k=7 validation latent
27. `data_val_filtered_tcgacptac` — validation data with PurIST/DeCAF
28. `fit_std_elbowk_tcgacptac` — NMF k=5 fit
29. `tar_data_filtered_elbowk_tcgacptac` — k=5 filtered data

**From supplement.Rmd (direct readRDS from store):**
30. `fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac` — NMF k=7 heatmap (no metadata)

**Total: 30 distinct targets store objects needed**

### RDS files from results/cv_grid/ (produced by inst/ scripts)

Referenced by BOTH 04_results_REVISED.Rmd and supplement.Rmd:
1. `adj_p_270_matrix.rds`
2. `adj_p_all_matrix.rds`
3. `unadj_p_270_matrix.rds` (results only, not supplement)
4. `hcor_270_matrix.rds`
5. `hcor_all_matrix.rds`
6. `master_rows_270.rds`
7. `master_rows_all.rds`
8. `k3_k7_summary.rds`
9. `lam300_summary.rds`
10. `production_summary.rds`

### Static PDF figure files

Referenced by paper.Rmd or supplement.Rmd:
1. `figures/model_schematic_final.pdf` — hand-drawn schematic
2. `figures/cv_grid/cv_cindex_by_k_primary.pdf` — C-index sensitivity
3. `figures/cv_grid/k3_k7_factor_correlation_heatmap.pdf` — k=3 vs k=7
4. `figures/cutpoint_curve_logrank_tcgacptac.pdf` — cutpoint selection
5. `figures/km_dichot/km_val_pooled_logrank_tcgacptac.pdf` — pooled KM
6. `figures/km_dichot/km_val_Dijk_logrank_tcgacptac.pdf`
7. `figures/km_dichot/km_val_Moffitt_GEO_array_logrank_tcgacptac.pdf`
8. `figures/km_dichot/km_val_PACA_AU_logrank_tcgacptac.pdf`
9. `figures/km_dichot/km_val_Puleo_array_logrank_tcgacptac.pdf`
10. `figures/km_dichot/subtype_overlap_pooled_logrank_tcgacptac.pdf`

### Input data files needed

The `tcgacptac` config trains on TCGA_PAAD + CPTAC. Validation uses Dijk, Moffitt, PACA_AU (array+seq), Puleo.
Each dataset has 3 files: `{name}.rds`, `{name}.survival_data.rds`, `{name}_subtype.csv`.
Also needed: `{name}.caf_subtype.rds` (used in validation for PurIST/DeCAF).

**Training data (2 cohorts × 4 files = 8 files):**
- `TCGA_PAAD.rds`, `.survival_data.rds`, `_subtype.csv`, `.caf_subtype.rds`
- `CPTAC.rds`, `.survival_data.rds`, `_subtype.csv`, `.caf_subtype.rds`

**Validation data (5 cohorts × 4 files = 20 files):**
- `Dijk.*`, `Moffitt_GEO_array.*`, `PACA_AU_array.*`, `PACA_AU_seq.*`, `Puleo_array.*`

**Derived data needed by inst/ scripts:**
- `data/original/cmbSubtypes.RData` (for cv_grid analysis scripts)

**Total: 29 data files, ~369 MB**

### R packages required for paper rendering

targets, stringr, dplyr, ggplot2, survival, survminer, cowplot, enrichplot,
magick, Seurat, gt, VAM, grid, viridis, ggplotify, NMF, gtable, knitr,
kableExtra, ggrepel, rticles, RColorBrewer, tibble, purrr

### R packages required for full pipeline re-run

DeSurv (+ Rcpp, RcppArmadillo, RcppEigen, survival, parallel, cvwrapr, preprocessCore),
NMF, pheatmap, ComplexHeatmap, ggplot2, cowplot, survival, survminer, dplyr,
data.table, clusterProfiler, org.Hs.eg.db (optional, for enrichment)

### LaTeX dependencies for paper compilation

`algorithm.sty`, `algpseudocode.sty`, `algorithmicx.sty` (in `paper/`),
`pnas-new.cls`, `pnasresearcharticle.sty`, `pnas.csl`

---

## Phase 1: Export pre-computed results from targets store

**Objective:** Extract everything the paper needs from the targets store into portable RDS files, eliminating the dependency on `targets` for paper rendering.

### Task 1.1: Create export script

**Files:**
- Create: `export_precomputed.R`

**What it does:** Reads the 30 target objects from the active store and saves them as individual RDS files in `results/precomputed/`. Also copies the 10 `results/cv_grid/` RDS files and 10 static PDFs into a self-contained `results/` structure.

```r
#!/usr/bin/env Rscript
# export_precomputed.R — Extract paper-needed objects from targets store
#
# Run ONCE from the old DeSurv-paper repo to export all pre-computed
# results needed for the clean repo.

library(targets)
tar_config_set(store = "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main")

out_dir <- "results/precomputed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Targets objects needed by the paper ─────────────────────────────────────
paper_targets <- c(
  # Main text figures (04_results_REVISED.Rmd)
  "sim_figs_by_scenario",
  "fig_bo_heat_tcgacptac",
  "fig_gene_overlap_heatmap_desurv_tcgacptac",
  "fig_gene_overlap_heatmap_std_desurvk_tcgacptac",
  "fig_variation_explained_tcgacptac",
  "fig_desurv_std_correlation_tcgacptac",
  "fig_hr_forest_tcgacptac",
  "fig_median_survival_desurv_tcgacptac",
  "fig_median_survival_std_desurvk_tcgacptac",
  # Main text data (for inline statistics)
  "val_latent_desurv_tcgacptac",
  "val_latent_std_desurvk_tcgacptac",
  "tar_k_selection_tcgacptac",
  "desurv_bo_results_tcgacptac",
  "tar_params_best_tcgacptac",
  "tar_data_filtered_tcgacptac",
  # Supplement figures
  "desurv_seed_fits_tcgacptac",
  "fig_residuals_tcgacptac",
  "fig_cophenetic_tcgacptac",
  "fig_silhouette_tcgacptac",
  "fit_std_tcgacptac",
  "fig_gene_overlap_heatmap_std_elbowk_tcgacptac",
  # Supplement tables
  "val_cindex_desurv_tcgacptac",
  "val_cindex_std_desurvk_tcgacptac",
  "val_cindex_std_elbowk_tcgacptac",
  "val_cindex_desurv_alpha0_tcgacptac",
  # Supplement data (adjusted HR analysis)
  "val_latent_desurv_alpha0_tcgacptac",
  "data_val_filtered_tcgacptac",
  "fit_std_elbowk_tcgacptac",
  "tar_data_filtered_elbowk_tcgacptac"
)

for (tgt in paper_targets) {
  message("Exporting: ", tgt)
  obj <- tar_read_raw(tgt)
  saveRDS(obj, file.path(out_dir, paste0(tgt, ".rds")))
}

# Special case: object with no metadata (must read from store directly)
message("Exporting: fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac (direct)")
obj <- readRDS(file.path(
  tar_config_get("store"), "objects",
  "fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac"
))
saveRDS(obj, file.path(out_dir, "fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac.rds"))

message("Done. Exported ", length(paper_targets) + 1, " objects to ", out_dir)
```

**Step 1:** Run `Rscript export_precomputed.R` in the current repo.
**Step 2:** Verify all 31 RDS files exist: `ls results/precomputed/*.rds | wc -l` → 31
**Step 3:** Verify total size is reasonable: `du -sh results/precomputed/` → ~350 MB

### Task 1.2: Verify cv_grid and static figure files

**Step 1:** Confirm all 10 cv_grid RDS files exist in `results/cv_grid/`:
```bash
for f in adj_p_270_matrix adj_p_all_matrix unadj_p_270_matrix hcor_270_matrix \
         hcor_all_matrix master_rows_270 master_rows_all k3_k7_summary \
         lam300_summary production_summary; do
  test -f "results/cv_grid/${f}.rds" && echo "OK: ${f}.rds" || echo "MISSING: ${f}.rds"
done
```

**Step 2:** Confirm all 10 static PDF figures exist:
```bash
for f in figures/model_schematic_final.pdf \
         figures/cv_grid/cv_cindex_by_k_primary.pdf \
         figures/cv_grid/k3_k7_factor_correlation_heatmap.pdf \
         figures/cutpoint_curve_logrank_tcgacptac.pdf \
         figures/km_dichot/km_val_pooled_logrank_tcgacptac.pdf \
         figures/km_dichot/km_val_Dijk_logrank_tcgacptac.pdf \
         figures/km_dichot/km_val_Moffitt_GEO_array_logrank_tcgacptac.pdf \
         figures/km_dichot/km_val_PACA_AU_logrank_tcgacptac.pdf \
         figures/km_dichot/km_val_Puleo_array_logrank_tcgacptac.pdf \
         figures/km_dichot/subtype_overlap_pooled_logrank_tcgacptac.pdf; do
  test -f "$f" && echo "OK: $f" || echo "MISSING: $f"
done
```

---

## Phase 2: Modify paper .Rmd files to read from precomputed RDS instead of targets

**Objective:** Replace all `tar_load()` / `tar_read()` calls with `readRDS()` from `results/precomputed/`. This removes the `targets` dependency from paper rendering entirely.

### Task 2.1: Create a helper function that replaces tar_load/tar_read

**Files:**
- Create: `paper/load_precomputed.R`

```r
# paper/load_precomputed.R
# Drop-in replacement for tar_load/tar_read that reads from precomputed RDS files.
# Source this file at the top of each .Rmd instead of library(targets).

PRECOMPUTED_DIR <- file.path("..", "results", "precomputed")

load_result <- function(name, envir = parent.frame()) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", path,
         "\nRun the pipeline or download pre-computed results first.")
  }
  obj <- readRDS(path)
  assign(name, obj, envir = envir)
  invisible(obj)
}

read_result <- function(name) {
  path <- file.path(PRECOMPUTED_DIR, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Pre-computed result not found: ", path)
  }
  readRDS(path)
}
```

### Task 2.2: Update paper.Rmd

**Files:**
- Modify: `paper/paper.Rmd`

**Changes:**
1. Remove the `targets-store` setup chunk (lines 120-126)
2. Add `source("paper/load_precomputed.R")` in the setup chunk
3. Remove `library(targets)` from setup

### Task 2.3: Update 04_results_REVISED.Rmd

**Files:**
- Modify: `paper/04_results_REVISED.Rmd`

**Changes:**
1. Replace `library(targets)` → `source("paper/load_precomputed.R")`
2. Replace all `tar_load(X)` → `load_result("X")`
3. Replace all `tar_read(X)` → `read_result("X")`
4. ~15 tar_load calls + ~5 tar_read calls (mostly inline R)

### Task 2.4: Update supplement.Rmd

**Files:**
- Modify: `paper/supplement.Rmd`

**Changes:**
1. Replace `library(targets)` → `source("paper/load_precomputed.R")`
2. Replace all `tar_load(X)` → `load_result("X")`
3. Replace `tar_read(X)` inline calls → `read_result("X")`
4. Replace the direct `readRDS(file.path(tar_config_get("store"), ...))` call (line 485-487) with `load_result("fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac")`
5. ~18 tar_load calls + 1 tar_read + 1 direct readRDS

### Task 2.5: Verify paper renders with precomputed results

**Step 1:** `Rscript -e 'rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd())'`
**Step 2:** `Rscript -e 'rmarkdown::render("paper/supplement.Rmd", knit_root_dir = getwd())'`
**Step 3:** Visually compare output PDFs with existing versions to confirm no regressions.

---

## Phase 3: Design the numbered-script pipeline

**Objective:** Write standalone R scripts that reproduce all 30 pre-computed objects + 10 cv_grid RDS files + 10 static PDFs from raw data + DeSurv package.

### Task 3.1: Create the pipeline scripts

The pipeline has this dependency structure (traced from targets):

```
01_install.R (install DeSurv)
    ↓
02_load_data.R (load TCGA+CPTAC training data + validation cohorts)
    ↓
03_bayesian_optimization.R (BO for DeSurv + BO for alpha=0 NMF)
    ↓
04_fit_models.R (fit DeSurv k=3, std NMF k=3, std NMF k=5, std NMF k=7)
    ↓
05_external_validation.R (project to validation cohorts, compute C-index, HR)
    ↓
06_sensitivity_analysis.R (k-sensitivity grid, cv_grid analysis)
    ↓
07_simulations.R (3 scenarios × 100 replicates × 6 methods = 1800 runs)
    ↓
08_figures.R (generate all figure objects + static PDFs)
    ↓
09_render_paper.R (compile manuscript + supplement)
```

**Files to create:**
- `code/00_helpers.R` — shared utilities (cache_or_compute, config loading)
- `code/01_install.R`
- `code/02_load_data.R`
- `code/03_bayesian_optimization.R`
- `code/04_fit_models.R`
- `code/05_external_validation.R`
- `code/06_sensitivity_analysis.R`
- `code/07_simulations.R`
- `code/08_figures.R`
- `code/09_render_paper.R`

Each script follows this pattern:

```r
#!/usr/bin/env Rscript
# code/03_bayesian_optimization.R
# Runs Bayesian optimization for DeSurv and standard NMF (alpha=0)
# Inputs:  results/precomputed/tar_data_tcgacptac.rds (from step 02)
# Outputs: results/precomputed/desurv_bo_results_tcgacptac.rds
#          results/precomputed/tar_k_selection_tcgacptac.rds
#          results/precomputed/tar_params_best_tcgacptac.rds
#          (+ alpha0 variants)
# Runtime: ~4-8 hours on HPC (30 cores), ~5 min in quick mode

source("code/00_helpers.R")

# ── Quick mode reduces iterations ─────────────────────────────────────────
if (CONFIG$quick) {
  ninit        <- 2
  bo_n_init    <- 4
  bo_n_iter    <- 4
  ncores_grid  <- 1
  parallel     <- FALSE
} else {
  ninit        <- 30
  bo_n_init    <- 50
  bo_n_iter    <- 100
  ncores_grid  <- CONFIG$ncores
  parallel     <- CONFIG$ncores > 1
}

# ── Load training data ────────────────────────────────────────────────────
tar_data <- readRDS(precomputed_path("tar_data_tcgacptac"))

# ── Run DeSurv BO ─────────────────────────────────────────────────────────
desurv_bo_results <- cache_or_compute("desurv_bo_results_tcgacptac", {
  bounds <- list(
    k_grid     = list(lower = 2L, upper = 12L, type = "integer"),
    alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
    lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
    nu_grid    = list(lower = 0, upper = 1, type = "continuous")
  )
  DeSurv::desurv_cv_bayesopt_refine(
    X = tar_data$ex, y = tar_data$sampInfo$time, d = tar_data$sampInfo$event,
    dataset = tar_data$sampInfo$dataset, samp_keeps = tar_data$samp_keeps,
    preprocess = TRUE, method_trans_train = "rank", engine = "warmstart",
    nfolds = 5, tol = 1e-5, maxit = 4000,
    coarse_bounds = bounds,
    bo_fixed = list(n_starts = ninit, ngene = 3000L, lambdaW_grid = 0, lambdaH_grid = 0),
    max_refinements = 0, tol_gain = 0.002, plateau = 1, top_k = 10,
    shrink_base = 0.3, importance_gain = 0.1,
    verbose = TRUE,
    parallel_grid = parallel, ncores_grid = ncores_grid
  )
})

# ... (similar for alpha0 BO, k selection, param extraction)
```

### Task 3.2: Create the shared helpers file

**Files:**
- Create: `code/00_helpers.R`

```r
# code/00_helpers.R — Shared utilities for all pipeline scripts

# ── Parse command-line args / environment ──────────────────────────────────
CONFIG <- list(
  quick   = identical(Sys.getenv("DESURV_QUICK"), "TRUE"),
  ncores  = as.integer(Sys.getenv("DESURV_NCORES", "1")),
  recompute = identical(Sys.getenv("DESURV_RECOMPUTE"), "TRUE")
)

if (CONFIG$quick) {
  message("=== QUICK MODE: reduced data/iterations for smoke testing ===")
}

# ── Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR    <- "results/precomputed"
CV_GRID_DIR    <- "results/cv_grid"
FIGURE_DIR     <- "figures"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CV_GRID_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

# ── Cache-or-compute pattern ──────────────────────────────────────────────
cache_or_compute <- function(name, expr) {
  path <- precomputed_path(name)
  if (!CONFIG$recompute && file.exists(path)) {
    message("  Loading cached: ", name)
    return(readRDS(path))
  }
  message("  Computing: ", name, " ...")
  result <- force(expr)
  saveRDS(result, path)
  message("  Saved: ", path)
  result
}
```

### Task 3.3: Create the Makefile

**Files:**
- Create: `Makefile`

```makefile
# DeSurv-paper Makefile
# Usage:
#   make quick          # Smoke test (~10 min, any laptop)
#   make all            # Full pipeline (requires HPC or patience)
#   make figures        # Regenerate figures from precomputed results
#   make paper          # Compile manuscript from precomputed results
#   make from-precomputed  # Figures + paper only (no pipeline)

NCORES ?= 1

.PHONY: all quick figures paper from-precomputed clean

# ── Quick mode (smoke test) ───────────────────────────────────────────────
quick:
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/01_install.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/02_load_data.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/03_bayesian_optimization.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/04_fit_models.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/05_external_validation.R
	DESURV_QUICK=TRUE DESURV_RECOMPUTE=TRUE DESURV_NCORES=1 Rscript code/07_simulations.R
	@echo "=== Quick mode complete. Pipeline runs end-to-end. ==="

# ── Full pipeline ─────────────────────────────────────────────────────────
all:
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/01_install.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/02_load_data.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/03_bayesian_optimization.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/04_fit_models.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/05_external_validation.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/06_sensitivity_analysis.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/07_simulations.R
	DESURV_RECOMPUTE=TRUE DESURV_NCORES=$(NCORES) Rscript code/08_figures.R
	Rscript code/09_render_paper.R
	@echo "=== Full pipeline complete. ==="

# ── From precomputed (reviewer fast path) ─────────────────────────────────
from-precomputed: figures paper

figures:
	Rscript code/08_figures.R

paper:
	Rscript code/09_render_paper.R

clean:
	rm -rf results/precomputed/*.rds results/cv_grid/*.rds figures/*.pdf
```

### Task 3.4: Create the master R script (Make-free alternative)

**Files:**
- Create: `run_pipeline.R`

```r
#!/usr/bin/env Rscript
# run_pipeline.R — Run the DeSurv analysis pipeline
#
# Usage:
#   Rscript run_pipeline.R                    # From precomputed (figures + paper)
#   Rscript run_pipeline.R --quick            # Quick smoke test
#   Rscript run_pipeline.R --full --ncores 8  # Full re-computation
#   Rscript run_pipeline.R --step 8           # Run from step 8 onward
#   Rscript run_pipeline.R --step 8 --only    # Run only step 8

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(optparse)

option_list <- list(
  make_option("--quick", action = "store_true", default = FALSE,
              help = "Quick mode: reduced data/iterations for smoke testing"),
  make_option("--full", action = "store_true", default = FALSE,
              help = "Full re-computation from raw data"),
  make_option("--step", type = "integer", default = 1L,
              help = "Start from this step number [default: %default]"),
  make_option("--only", action = "store_true", default = FALSE,
              help = "Run only the specified step"),
  make_option("--ncores", type = "integer", default = 1L,
              help = "Number of cores for parallel steps [default: %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

# Set environment variables consumed by individual scripts
if (opts$quick) {
  Sys.setenv(DESURV_QUICK = "TRUE", DESURV_RECOMPUTE = "TRUE", DESURV_NCORES = "1")
} else if (opts$full) {
  Sys.setenv(DESURV_RECOMPUTE = "TRUE", DESURV_NCORES = as.character(opts$ncores))
}

steps <- c(
  "code/01_install.R",
  "code/02_load_data.R",
  "code/03_bayesian_optimization.R",
  "code/04_fit_models.R",
  "code/05_external_validation.R",
  "code/06_sensitivity_analysis.R",
  "code/07_simulations.R",
  "code/08_figures.R",
  "code/09_render_paper.R"
)

# Default: just figures + paper (from precomputed)
if (!opts$quick && !opts$full && opts$step == 1L) {
  opts$step <- 8L
  message("No --quick or --full specified. Running from step 8 (figures + paper).")
}

for (i in seq_along(steps)) {
  if (i < opts$step) next
  if (opts$only && i != opts$step) next
  cat(sprintf("\n=== Step %d/%d: %s ===\n", i, length(steps), steps[i]))
  t0 <- Sys.time()
  source(steps[i], local = new.env(parent = globalenv()))
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("    Completed in %.1f minutes\n", as.numeric(elapsed)))
}

cat("\n=== Pipeline complete ===\n")
```

---

## Phase 4: Build the clean repo structure

**Objective:** Assemble the clean `rashidlab/DeSurv-paper` repository with only necessary files.

### Task 4.1: Initialize clean repo locally

```bash
mkdir -p ~/Downloads/rashidlab-DeSurv-paper
cd ~/Downloads/rashidlab-DeSurv-paper
git init
```

### Task 4.2: Copy necessary files from current repo

**Pipeline code (new, from Phase 3):**
```
code/00_helpers.R
code/01_install.R
code/02_load_data.R
code/03_bayesian_optimization.R
code/04_fit_models.R
code/05_external_validation.R
code/06_sensitivity_analysis.R
code/07_simulations.R
code/08_figures.R
code/09_render_paper.R
```

**Helper R functions (subset of current R/):**
Only the functions actually called by the pipeline scripts and paper .Rmd files.
Exact list to be determined during Phase 3 implementation, but likely includes:
```
code/helpers/
  load_data.R                  # load_data(), load_data_internal()
  cluster_alignment.R          # align_clusters(), cluster-comparison logic
  figure_functions.R           # extracted from figure_targets.R — the build/plot functions
  enrichment_map.R             # ORA analysis functions
  plot_survival.R              # splot_cutpoint, KM plotting
  compare_models.R             # compute_hrs, variance-survival functions
  bo_helpers.R                 # select_bo_k_by_cv_se, standardize_bo_params
  get_top_genes.R              # gene ranking extraction
  cluster_validation_scores.R  # validation score computation
```

**Simulation code:**
```
code/helpers/simulation_functions/   # entire directory from R/simulation_functions/
sim_figs.R                           # simulation figure building
```

**Pre-computed results (from Phase 1):**
```
results/precomputed/    # 31 RDS files, ~350 MB → Zenodo
results/cv_grid/        # 10 RDS files, ~1 MB
```

**Static figures:**
```
figures/model_schematic_final.pdf
figures/cv_grid/cv_cindex_by_k_primary.pdf
figures/cv_grid/k3_k7_factor_correlation_heatmap.pdf
figures/cutpoint_curve_logrank_tcgacptac.pdf
figures/km_dichot/km_val_pooled_logrank_tcgacptac.pdf
figures/km_dichot/km_val_Dijk_logrank_tcgacptac.pdf
figures/km_dichot/km_val_Moffitt_GEO_array_logrank_tcgacptac.pdf
figures/km_dichot/km_val_PACA_AU_logrank_tcgacptac.pdf
figures/km_dichot/km_val_Puleo_array_logrank_tcgacptac.pdf
figures/km_dichot/subtype_overlap_pooled_logrank_tcgacptac.pdf
```

**Data (29 files, ~369 MB → Zenodo or download script):**
```
data/original/TCGA_PAAD.rds  (+ .survival_data.rds, _subtype.csv, .caf_subtype.rds)
data/original/CPTAC.*
data/original/Dijk.*
data/original/Moffitt_GEO_array.*
data/original/PACA_AU_array.*
data/original/PACA_AU_seq.*
data/original/Puleo_array.*
data/original/cmbSubtypes.RData
```

**Paper source:**
```
paper/paper.Rmd                   # modified (Phase 2)
paper/04_results_REVISED.Rmd      # modified (Phase 2)
paper/02_introduction_REVISED.Rmd
paper/03_methods_REVISED.Rmd
paper/05_discussion_REVISED.Rmd
paper/supplement.Rmd              # modified (Phase 2)
paper/supp_methods.Rmd
paper/load_precomputed.R          # new (Phase 2)
paper/references_30102025.bib
paper/pnas.csl
paper/pnas-new.cls
paper/pnasresearcharticle.sty
paper/algorithm.sty
paper/algpseudocode.sty
paper/algorithmicx.sty
paper/gene_lists_top270_k3.csv    # if referenced
```

**Top-level files:**
```
README.md                # new — reproduction instructions
LICENSE
Makefile                 # new (Phase 3)
run_pipeline.R           # new (Phase 3)
DESCRIPTION              # for renv dependency tracking
renv.lock                # generated by renv
.Rprofile                # renv bootstrap
.gitignore               # exclude results/precomputed/, data/original/, renv/library/
```

**HPC scripts:**
```
slurm/run_full_pipeline.sh       # single-job submission
slurm/run_step_by_step.sh        # chained job submissions
```

### Task 4.3: Files explicitly NOT included (and why)

| File/Directory | Reason for exclusion |
|---|---|
| `_targets.R`, `_targets_*.R` | Replaced by numbered scripts |
| `targets_*.R` (configs, setup, common_pipeline) | Pipeline internals, replaced |
| `targets_setup.R` | Library loading / crew controllers — not needed |
| `_targets/`, `store_*` (all 6 stores) | Replaced by `results/precomputed/` |
| `R/*.R` (64 files) | Only ~10-12 helper functions extracted to `code/helpers/` |
| `inst/*.R` (35 scripts) | Standalone scripts; cv_grid outputs already in `results/cv_grid/` |
| `old_sim_pipelines/` | Superseded |
| `logs/` | Runtime artifacts |
| `rebuild_figures.R`, `regenerate_figures.R` | Workarounds for ggplot2 version mismatch |
| `*.md` review docs (CODE_REVIEW, NARRATIVE_ARC, PNAS_REVIEW, etc.) | Internal review artifacts |
| `*.pptx`, `*.py`, `Rplots.pdf`, hex sticker | Non-essential |
| `docs/plans/` | Development planning docs |
| `review/` | Review artifacts |
| `data/derv/` (Elyada UMAPs, etc.) | scRNA-seq analysis commented out in supplement |
| `data/original/bladder*`, `Grunwald*`, `Hayashi*`, `Linehan*`, `Olive*` | Not in paper |
| `paper/old/`, `paper/*_files/` | Build artifacts |
| `paper/04_results.Rmd` (non-REVISED) | Superseded by REVISED version |
| `paper/paper.docx`, `paper.tex`, `paper.log` | Build outputs |
| `sim_figs.R`, `sim_figs_by_scenario.rds`, `sim_results_table.rds` | Top-level copies, superseded |
| `supp_bo.Rmd` | Not included in final paper |
| `CLAUDE.md`, `.claude/` | Development tooling |
| `local_slurm/` | Local development configs |
| `untitled_consensus_cluster/`, `NMF_*/` | Temp directories |
| `tests/` | Pipeline tests — rewrite for new pipeline |
| `data/original/IMmotion150*.rds`, `IMVigor210*.rds` | Bladder data, not in paper |

---

## Phase 5: Set up renv and Docker

### Task 5.1: Initialize renv

```r
# In the clean repo directory:
renv::init()
# Install all needed packages, then:
renv::snapshot()
```

### Task 5.2: Create Dockerfile (essential for reviewers)

**Files:**
- Create: `Dockerfile`

```dockerfile
FROM rocker/r-ver:4.3.1

# System dependencies for RcppArmadillo, preprocessCore, etc.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    liblapack-dev libblas-dev gfortran \
    libmagick++-dev libpoppler-cpp-dev \
    texlive-latex-base texlive-latex-extra texlive-fonts-recommended \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /desurv-paper

# Copy renv lock first (caching layer)
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/activate.R
RUN Rscript -e 'renv::restore()'

# Copy everything else
COPY . .

# Default: regenerate figures and compile paper from precomputed results
CMD ["make", "from-precomputed"]
```

### Task 5.3: Write README.md

**Files:**
- Create: `README.md`

Key sections:
1. **Quick start** (3 commands: clone, restore, make from-precomputed)
2. **Quick mode** (smoke test the full pipeline in ~10 min)
3. **Full reproduction** (HPC instructions)
4. **Docker** (single-command reproduction)
5. **Data availability** (Zenodo DOI, GEO accessions)
6. **Repository structure** (directory layout)
7. **Archive note** (link to Amber's repos for development history)

---

## Phase 6: Clean DeSurv package repo

### Task 6.1: Create rashidlab/DeSurv

**Steps:**
1. Copy package source from `~/Downloads/DeSurv/`
2. Merge `20260107bugfix` branch fixes into `main`
3. Update DESCRIPTION: version → 2.0.0, URL → rashidlab, BugReports → rashidlab
4. Update AUTHORS to include all co-authors
5. Verify: `R CMD check --as-cran`
6. Verify: `Rscript -e 'devtools::test()'` — all 14 test files pass
7. Build pkgdown site or write comprehensive README

**Files to include:**
```
R/                   # 21 source files
src/                 # functions.cpp, RcppExports.cpp
tests/testthat/      # 14 test files
man/                 # 25+ .Rd files
vignettes/           # desurv-intro.Rmd
DESCRIPTION
NAMESPACE
LICENSE
README.md
.Rbuildignore
```

**Files to exclude:**
```
example_cv_bayesopt.R  # not part of package
.claude/               # development tooling
DeSurv.Rproj           # optional, keep if desired
```

---

## Phase 7: Data hosting and Zenodo deposit

### Task 7.1: Create Zenodo deposit

Upload to Zenodo:
1. `data/original/` (29 files, ~369 MB) — or subset if licensing restricts
2. `results/precomputed/` (31 RDS files, ~350 MB)
3. `results/cv_grid/` (10 RDS files, ~1 MB)

Get DOI. Reference in README.md and paper.

### Task 7.2: Create data download script

**Files:**
- Create: `code/00_download_data.R`

```r
#!/usr/bin/env Rscript
# Downloads pre-computed results and data from Zenodo
# Only needed if you want to skip full pipeline re-computation

ZENODO_DOI <- "10.5281/zenodo.XXXXXXX"  # Fill after deposit
# ... download logic using httr or curl
```

### Task 7.3: Audit data provenance

**Files:**
- Create: `data/README.md`

Document for each dataset:
- Source (GEO, CPTAC, ICGC, etc.)
- Accession numbers
- Processing steps applied in `load_data_internal()`
- Any access restrictions

---

## Phase 8: Verification

### Task 8.1: Fresh-machine test

1. Clone `rashidlab/DeSurv-paper` to a clean machine (or Docker)
2. Run `make from-precomputed` — verify paper compiles
3. Run `make quick` — verify pipeline runs end-to-end
4. Compare rendered PDFs with submission versions

### Task 8.2: Lab member test

Have someone other than Amber:
1. Clone both repos
2. Follow README exactly
3. Report any failures or unclear instructions
4. Time the quick mode and from-precomputed paths

### Task 8.3: Tag releases

```bash
# In rashidlab/DeSurv:
git tag -a v2.0.0 -m "Submission version for DeSurv manuscript"
git push origin v2.0.0

# In rashidlab/DeSurv-paper:
git tag -a v1.0-submission -m "Manuscript submission version"
git push origin v1.0-submission
```

---

## Summary: What goes where

| Current location | Clean repo destination | Size |
|---|---|---|
| `store/objects/` (30 needed) | `results/precomputed/` → Zenodo | ~350 MB |
| `results/cv_grid/` (10 RDS) | `results/cv_grid/` | ~1 MB |
| `figures/` (10 static PDFs) | `figures/` | ~5 MB |
| `data/original/` (29 files) | `data/original/` → Zenodo | ~369 MB |
| `R/` (64 files) | `code/helpers/` (~12 files) | ~150 KB |
| `R/simulation_functions/` | `code/helpers/simulation_functions/` | ~20 KB |
| `R/figure_targets.R` (2590 lines) | `code/helpers/figure_functions.R` (extracted) | ~80 KB |
| `targets_common_pipeline.R` | Replaced by `code/02-08*.R` | N/A |
| `paper/*.Rmd` + templates | `paper/` (modified for precomputed) | ~2 MB |
| `_targets*.R`, `targets_*.R` | Not included | N/A |
| All 6 store directories | Not included | Saves ~37 GB |

**Total clean repo:** ~5 MB code + ~2 MB paper + ~5 MB static figures + references to Zenodo for data (~720 MB)
