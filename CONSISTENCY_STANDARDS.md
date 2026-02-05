# Consistency Standards

Guidelines to ensure alignment between code, configuration, documentation, and paper claims.

## Quick Reference Checklist

Before committing changes, verify:

- [ ] **Config changes** → Update CLAUDE.md "Current Branch Configuration" section
- [ ] **Function default changes** → Update roxygen docs AND DEFAULTS.md
- [ ] **Figure code changes** → Re-run targets, verify paper claims still hold
- [ ] **Paper text changes** → Ensure code produces what text claims
- [ ] **New targets** → Add to MANUSCRIPT_PIPELINE_ANALYSIS.md figure mapping

---

## 1. Configuration Consistency

### Single Source of Truth

| Config Type | Primary File | Must Sync With |
|-------------|--------------|----------------|
| BO hyperparameters | `targets_bo_configs.R` | CLAUDE.md, paper methods |
| Run parameters | `targets_run_configs.R` | CLAUDE.md |
| Validation datasets | `targets_val_configs.R` | Paper Figure 6 claims |
| Store path | `paper/_targets.yaml` | `paper/paper.Rmd` params |

### When Changing Configs

1. **Update the config file** (e.g., `targets_bo_configs.R`)
2. **Run validation**: `Rscript -e 'source("R/targets_config.R"); validate_desurv_configs(...)'`
3. **Invalidate affected targets**: `tar_invalidate(names = affected_targets)`
4. **Update documentation**:
   - `CLAUDE.md` → "Current Branch Configuration" section
   - `CHANGELOG.md` → Add entry with rationale
5. **If paper-relevant**: Check if paper claims still hold after re-running

### Config Validation Command

```r
# Run before committing config changes
source("R/targets_config.R")
bo <- targets_bo_configs()
run <- targets_run_configs()
val <- targets_val_configs()
validate_desurv_configs(bo, run, val)
```

---

## 2. Function Defaults Consistency

### Documentation Standard

Every function with user-facing defaults must have:

```r
#' @param param_name Description. Default: \code{actual_default_value}
```

Example:
```r
#' Run Bayesian Optimization for DeSurv
#'
#' @param bo_n_init Number of initial random samples. Default: \code{50}
#' @param bo_n_iter Number of BO iterations. Default: \code{100}
#' @param nfold Number of CV folds. Default: \code{5}
```

### Defaults Registry

Maintain `DEFAULTS.md` with all function defaults:

```markdown
## desurv_default_bo_config()
| Parameter | Default | Rationale |
|-----------|---------|-----------|
| bo_n_init | 50 | Balance exploration vs computation |
| bo_n_iter | 100 | Sufficient for convergence |
```

### When Changing Defaults

1. Update function code
2. Update roxygen `@param` documentation
3. Update `DEFAULTS.md`
4. Run `devtools::document()` if in DeSurv package
5. Add to `CHANGELOG.md`

---

## 3. Paper-Code Consistency

### Figure-Target Mapping

Every figure panel must have a documented target:

| Figure | Panel | Target | Claim |
|--------|-------|--------|-------|
| Fig 2 | A | `fig_residuals_tcgacptac` | Residuals decrease smoothly |
| Fig 2 | B | `fig_cophenetic_tcgacptac` | Cophenetic fluctuates |
| Fig 2 | C | `fig_silhouette_tcgacptac` | Silhouette highest at k=2-3 |
| Fig 2 | D | `fig_bo_heat_tcgacptac` | C-index surface shows bimodal structure |
| Fig 2 | E | `sim_figs_by_scenario` | DeSurv recovers true k=3 |
| Fig 3 | A-B | `sim_figs_by_scenario` | DeSurv > NMF for C-index and precision |
| Fig 4 | A-B | `fig_gene_overlap_heatmap_*` | Factor-program correlations |
| Fig 4 | C | `fig_variation_explained_*` | Variance vs survival contribution |
| Fig 5 | A | `fig_hr_forest_tcgacptac` | Forest plot of HRs |
| Fig 6 | A-C | `fig_extval_*` | Generalization to 5 datasets |

### Quantitative Claims Checklist

For each quantitative claim in the paper:

1. **Identify the claim** (e.g., "DeSurv recovers k=3 with concentrated distribution")
2. **Find the generating code** (e.g., `sim_figs.R:plot_sim_k_hist()`)
3. **Verify the data supports the claim** (mode of histogram = 3)
4. **Document the validation** in MANUSCRIPT_PIPELINE_ANALYSIS.md

### Before Submitting Paper

Run this audit:
```r
# Verify all figure targets are up-to-date
library(targets)
tar_config_set(store = "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full")

# Check target status
meta <- tar_meta()
figure_targets <- meta[grepl("^fig_", meta$name), ]
outdated <- figure_targets[figure_targets$time < file.info("R/figure_targets.R")$mtime, ]
if (nrow(outdated) > 0) {
  warning("Outdated figure targets: ", paste(outdated$name, collapse = ", "))
}
```

---

## 4. Store Consistency

### Current Store: `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`

All components must reference the same store:

| Component | File | Setting |
|-----------|------|---------|
| Main pipeline | `submit_targets.R` | Auto-generated from branch |
| Simulation pipeline | `_targets_sims_local.sh` | Explicit in script |
| Paper rendering | `paper/_targets.yaml` | `store:` key |
| Paper params | `paper/paper.Rmd` | `params$tar_store` |

### When Changing Stores

1. Update ALL four locations above
2. Verify with: `grep -r "store_PKG_VERSION" . --include="*.R" --include="*.Rmd" --include="*.yaml" --include="*.sh"`
3. Add to CHANGELOG.md

---

## 5. Validation Hooks

### Pre-Commit Checks (Manual)

Before committing, run:

```bash
# 1. R syntax check
Rscript -e 'lapply(list.files("R", pattern="\\.R$", full.names=TRUE), parse)'

# 2. Config validation
Rscript -e 'source("R/targets_config.R"); bo <- targets_bo_configs(); validate_config_labels(bo, "bo")'

# 3. Targets validation
Rscript -e 'targets::tar_validate()'
```

### Pre-Push Checks (Critical)

Before pushing, verify:

```bash
# 1. All tests pass
Rscript -e 'testthat::test_dir("tests/testthat")'

# 2. Paper renders without errors
Rscript -e 'rmarkdown::render("paper/paper.Rmd", run_pandoc=FALSE)'

# 3. No debug statements
grep -r "browser()" R/ --include="*.R" && echo "FAIL: Remove browser() statements"
```

---

## 6. Change Documentation

### CHANGELOG.md Format

```markdown
## [YYYY-MM-DD]

### Category (Config/Code/Paper/Figures)
- **What changed**: Brief description
  - File(s): `path/to/file.R`
  - Rationale: Why this change was made
  - Impact: What targets/figures are affected
```

### Commit Message Format

```
<type>: <short description>

<body with details>

Affects: <list of affected targets or figures>
```

Types: `config`, `fix`, `feat`, `docs`, `paper`, `sim`

Example:
```
config: update tcgacptac BO bounds to match student

- k_grid: 2-15 → 2-12
- alpha_grid: 0-0.95 → 0-1.0
- bo_n_init: 20 → 50, bo_n_iter: 50 → 100

Affects: desurv_bo_results_tcgacptac, fig_bo_heat_tcgacptac, Fig 2D
```

---

## 7. File Modification Triggers

When you modify these files, also update:

| Modified File | Also Update |
|---------------|-------------|
| `targets_bo_configs.R` | CLAUDE.md, CHANGELOG.md, invalidate BO targets |
| `R/figure_targets.R` | MANUSCRIPT_PIPELINE_ANALYSIS.md if new figures |
| `paper/04_results.Rmd` | Verify all `tar_load()` targets exist |
| `_targets_sims.R` | CHANGELOG.md, check simulation scenarios |
| `paper/_targets.yaml` | `paper/paper.Rmd` params |
| Any `.R` file | Run `tar_validate()` |

---

## 8. Consistency Verification Script

Save as `scripts/verify_consistency.R`:

```r
#!/usr/bin/env Rscript
# Run: Rscript scripts/verify_consistency.R

library(targets)

cat("=== DeSurv Consistency Verification ===\n\n")

# 1. Store consistency
cat("1. Checking store consistency...\n")
store_refs <- c(
  yaml = yaml::read_yaml("paper/_targets.yaml")$main$store,
  rmd = rmarkdown::yaml_front_matter("paper/paper.Rmd")$params$tar_store
)
if (length(unique(store_refs)) == 1) {
  cat("   ✓ All store references match:", unique(store_refs), "\n")
} else {
  cat("   ✗ MISMATCH:\n")
  print(store_refs)
}

# 2. Config validation
cat("\n2. Validating configs...\n")
tryCatch({
  source("R/targets_config.R")
  bo <- targets_bo_configs()
  cat("   ✓ BO configs valid\n")
}, error = function(e) cat("   ✗ BO config error:", e$message, "\n"))

# 3. Pipeline validation
cat("\n3. Validating targets pipeline...\n")
tryCatch({
  tar_validate(script = "_targets.R")
  cat("   ✓ Main pipeline valid\n")
}, error = function(e) cat("   ✗ Pipeline error:", e$message, "\n"))

# 4. Debug statement check
cat("\n4. Checking for debug statements...\n")
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
debug_found <- FALSE
for (f in r_files) {
  content <- readLines(f, warn = FALSE)
  if (any(grepl("browser\\(\\)", content))) {
    cat("   ✗ browser() found in:", f, "\n")
    debug_found <- TRUE
  }
}
if (!debug_found) cat("   ✓ No debug statements found\n")

# 5. Figure target status
cat("\n5. Checking figure targets...\n")
tryCatch({
  tar_config_set(store = store_refs[[1]])
  meta <- tar_meta()
  fig_targets <- meta[grepl("^fig_", meta$name), ]
  cat("   Figure targets:", nrow(fig_targets), "\n")
  errored <- fig_targets[fig_targets$error != "", ]
  if (nrow(errored) > 0) {
    cat("   ✗ Errored targets:", paste(errored$name, collapse = ", "), "\n")
  } else {
    cat("   ✓ No errored figure targets\n")
  }
}, error = function(e) cat("   Could not check targets:", e$message, "\n"))

cat("\n=== Verification Complete ===\n")
```

---

## Summary

**The core principle**: Every change should be traceable from code → config → documentation → paper.

**Enforcement**: Run `scripts/verify_consistency.R` before any commit touching configs, figures, or paper text.

**Documentation**: Always update CHANGELOG.md when making substantive changes.
