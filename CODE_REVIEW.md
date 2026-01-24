# Code Review Findings

Comprehensive code review of DeSurv-paper and DeSurv package repositories. Issues organized by severity and category.

---

## Critical Bugs

### 1. Variable Assignment Bug in compute_metrics.R
**File**: `R/compute_metrics.R:9`
**Severity**: Critical
**Category**: Copy/paste error

Line 9 incorrectly assigns `fit$meanZ` to `sdZ` instead of `fit$sdZ`:
```r
sdZ = fit$meanZ[,ns]  # WRONG - should be fit$sdZ
```

This causes incorrect metric calculations. The variable `sdZ` is used for standardization (line 12) and linear predictor calculation (line 13).

**Fix**: Change to `sdZ = fit$sdZ[,ns]`

---

### 2. Unused Filter Variable in select_best_init.R
**File**: `R/select_best_init.R:5`
**Severity**: Critical
**Category**: Logic error

Line 5 creates filtered dataframe `keep` but never uses it:
```r
keep <- df[!df$flag_nan & is.finite(df$loss), ]
# Lines 11-15 still use original `df` instead of `keep`
```

Invalid initializations (NaN/infinite loss) may be selected as "best".

**Fix**: Replace uses of `df` after line 5 with `keep`.

---

### 3. Beta Backtracking Logic Inverted (DeSurv Package)
**File**: `../DeSurv/src/functions.cpp:472-473`
**Severity**: Critical
**Category**: Logic error

The beta update backtracking checks if loss is **lower** (better), but `surv_loss` is a likelihood we want to **maximize**:
```cpp
if (!std::isfinite(surv_loss_candidate) ||
    surv_loss_candidate < surv_loss_prev) {
```

This triggers unnecessary backtracking when the fit improves, slowing convergence.

**Fix**: Change to `surv_loss_candidate > surv_loss_prev`

---

### 4. Debug Breakpoints in Production Code
**Files**: `targets_common_pipeline.R:1691,1803,2075,2097,2119,2141`
**Severity**: Critical
**Category**: Debug code left in

Multiple `browser()` statements will halt the pipeline and wait for user input.

**Fix**: Remove all `browser()` calls from production code.

---

## High Priority Issues

### 5. Missing Function Definitions in bo_helpers.R
**File**: `R/bo_helpers.R:62,66`
**Severity**: High
**Category**: Missing dependency

Functions `collect_bo_diagnostics()` and `compute_bo_eval_se()` are called but don't exist. The `tryCatch` silently masks this, causing fallback to simpler k-selection logic.

---

### 6. Undefined Global Variables
**Files**: `R/run_seed_fits.R:2,18-21`, `R/create_filepath_coldstart_runs.R:10`
**Severity**: High
**Category**: Scope error

References to `NINIT_FULL`, `TOL`, `MAXIT`, `NGENE` assume these globals exist. Functions fail when called outside targets pipeline.

**Fix**: Pass as function parameters.

---

### 7. Missing Gene Subset Check in CV (DeSurv Package)
**File**: `../DeSurv/R/cv_helpers.R:231`
**Severity**: High
**Category**: Array bounds

When genes are missing from validation matrix, R creates NA rows that are silently replaced with zeros, corrupting predictions.

---

### 8. Path Injection Risk
**File**: `R/load_data_internal.R:6-8`
**Severity**: High
**Category**: Security

`dataname` parameter directly concatenated into file paths:
```r
dat = readRDS(paste0("data/original/", dataname, ".rds"))
```

Path traversal possible (e.g., `dataname = "../../etc/passwd"`).

**Fix**: Validate `dataname` is alphanumeric only before use.

---

### 9. DEFAULT_NINIT Reference Before Definition
**File**: `targets_setup.R:110`
**Severity**: High
**Category**: Race condition

`DEFAULT_NINIT` referenced before defined in `_targets.R`, causing cryptic errors.

---

### 10. Resource Overcommitment
**File**: `targets_setup.R:89-131`
**Severity**: High
**Category**: Resource management

Controllers request 202 workers × 30-100 CPUs = 6,000-20,000 CPUs total. Slurm jobs will fail or queue indefinitely.

---

## Medium Priority Issues

### 11. Non-Finite Validation Inconsistency (DeSurv Package)
**File**: `../DeSurv/R/predict_methods.R:71-106`
**Severity**: Medium
**Category**: Data integrity

Training handles non-finite values via replacement, but prediction rejects them strictly. Breaks symmetry between workflows.

---

### 12. W Column Normalization Near-Zero Columns (DeSurv Package)
**File**: `../DeSurv/src/functions.cpp:232-236`
**Severity**: Medium
**Category**: Numerical stability

Near-zero W columns (slightly above ε=1e-12) can produce large H/beta values, destabilizing optimization.

**Fix**: Increase threshold to 1e-8 or re-initialize degenerate columns.

---

### 13. Unsafe Cox Gradient Scaling (DeSurv Package)
**File**: `../DeSurv/src/functions.cpp:201-205`
**Severity**: Medium
**Category**: Numerical stability

`cox_scale` can reach 1e6 when Cox gradient norm is small, massively amplifying noisy gradients.

**Fix**: Lower cap from 1e6 to 100-1000.

---

### 14. Zero-Padding Creates Information Leakage
**File**: `R/load_data.R:17-19`
**Severity**: Medium
**Category**: Data integrity

Missing genes zero-filled when combining datasets, creating artificial patterns NMF may learn.

---

### 15. Transform Drift in Validation
**File**: `R/targets_config.R:582-653`
**Severity**: Medium
**Category**: Data integrity

Rank transforms on validation data compute new ranks independently, creating distribution shift.

---

### 16. Error Suppression in Parallel Execution
**File**: `targets_common_pipeline.R:501-543`
**Severity**: Medium
**Category**: Error handling

Seed fits fail silently; only checked at aggregate level. Difficult to diagnose failures.

---

### 17. Division by Small SD in BO (DeSurv Package)
**File**: `../DeSurv/R/desurv_bayesopt.R:336-347`
**Severity**: Medium
**Category**: Numerical stability

EI acquisition uses threshold ε≈1.5e-8 which may be too conservative. Consider increasing to 1e-6.

---

### 18. Configuration Function Mismatch
**Files**: `targets_run_configs.R`, `local_slurm/targets_configs.R`
**Severity**: Medium
**Category**: Configuration

Main and local_slurm configs have incompatible structure. Cannot easily switch environments.

---

## Architectural Issues

### 19. Tight Coupling Through Bundles
**File**: `targets_common_pipeline.R:1-21,447-467`

16 fields bundled together; changing any field invalidates entire bundle, forcing expensive recomputation.

---

### 20. Implicit File Dependencies
**File**: `targets_common_pipeline.R:1123-1137`

File writes are side effects not tracked by targets. Manual cleanup required.

---

### 21. Triple Assignment Typo
**File**: `targets_common_pipeline.R:327`

```r
bo_fixed <- listbo_fixed <- listbo_fixed <- list(...)
```

Middle assignments are no-ops.

---

### 22. No Schema Validation for Configs
**File**: `R/targets_config.R:353-403`

Validation checks existence only, not types or ranges. Invalid configs propagate to runtime.

---

### 23. Error Continue Mode
**File**: `targets_setup.R:164`

`tar_option_set(error = "continue")` lets pipeline continue after critical failures.

---

## Test Coverage Gaps

Functions lacking adequate test coverage:
- `compute_metrics.R` - Critical metric calculations
- `select_best_init.R` - Model selection
- `load_data_internal.R` - Data validation
- `run_seed_fits.R` - Convergence failure handling
- `split_train_validation()` - Data splitting

---

## DeSurv Package: Critical Invariants

### Data Constraints
```
X: p×n, all finite, all >= 0, Xnorm > 0
y: n×1, all finite, all > 0
d: n×1, all in {0,1}, sum(d) > 0
k: scalar, 1 <= k <= min(p,n)
```

### Hyperparameter Constraints
```
alpha  : [0, 1)
nu     : [0, 1]
lambda, lambdaW, lambdaH: >= 0
```

### NMF Constraints
```
W, H: all non-negative
W columns L2-normalized after each update
H restored to original sample order before returning
```

---

## DeSurv Package: Key Code Paths

### Training Pipeline
```
desurv_fit()
  → .validate_desurv_hyperparams()
  → init() [multi-start, parallel optional]
    → optimize_loss_cpp() × ninit [short runs]
  → optimize_loss_cpp() [full optimization]
    → update_H_cpp()
    → update_W_damped_backtrack()
    → update_beta_cpp()
  → return desurv_fit object
```

### Cross-Validation Pipeline
```
desurv_cv()
  → .desurv_cv_grid_warmstart()
    → .desurv_make_folds_stratified()
    → desurv_alpha_warmstart() per job
  → select best (max or 1se rule)
  → desurv_fit() on full data
```

---

## Recommendations

### Immediate (Critical)
1. Fix `sdZ` assignment in `compute_metrics.R`
2. Use `keep` variable in `select_best_init.R`
3. Fix beta backtracking logic in C++
4. Remove all `browser()` statements

### High Priority
5. Fix `DEFAULT_NINIT` reference order
6. Add path validation to `load_data_internal()`
7. Recalculate resource controller CPUs
8. Pass globals as function parameters

### Medium Priority
9. Increase numerical stability thresholds in C++
10. Add schema validation for configs
11. Standardize naming conventions
12. Add integration tests for full pipeline

### Long-term
13. Split large bundles for granular invalidation
14. Add logging for debugging
15. Create config templates with JSON schema
16. Document error codes for k-selection
