# Bug 1 Impact Analysis: compute_metrics.R sdZ Assignment Error

## Bug Description
**File:** `R/compute_metrics.R`, Line 9
**Issue:** `sdZ = fit$meanZ[,ns]` should be `sdZ = fit$sdZ[,ns]`

## Investigation Summary

### What the Bug Does
When active, this bug would cause:
1. **Line 12:** `sweep(XtW, 2, sdZ, FUN="/")` divides by mean instead of standard deviation
2. **Line 13:** `lp <- XtW %*% (beta*sdZ)` multiplies beta by mean instead of standard deviation
3. **Result:** Incorrect linear predictor values → incorrect C-index calculations

### Code Path Analysis

#### 1. Simulation Pipeline (NOT AFFECTED)
The simulation pipeline in `_targets_sims.R` uses its own C-index calculation:
```r
# Lines 403-423 in _targets_sims.R
compute_dataset_cindex <- function(fit, dataset) {
  Z <- t(X_sub) %*% W_sub
  risk <- drop(Z %*% beta)
  survival::concordance(surv_obj ~ preds$risk_score)
}
```
This does **NOT** call `compute_metrics.R` - it computes risk directly as `Z %*% beta`.

#### 2. Real Data Validation Pipeline (NOT AFFECTED)
The validation pipeline uses `R/predict_validation_scores.R`:
```r
# Lines 145-158
if (!is.null(stats)) {
    mean_vec <- stats$mean[factor_track]
    sd_vec <- stats$sd[factor_track]  # Correctly uses sd_vec
    Z_centered <- sweep(Z, 2, mean_vec, FUN = "-")
    Z_scaled <- sweep(Z_centered, 2, sd_vec, FUN = "/")
}
```
This code correctly uses separate `mean_vec` and `sd_vec` variables.

#### 3. compute_metrics_bundle.R (UNUSED)
- `compute_metrics_bundle.R` calls `compute_metrics()` at lines 22 and 36
- **However**, no targets in the current pipeline call `compute_metrics_bundle()`
- The direct call in `_targets.R:399` is commented out

### Verification
```bash
# No active targets use compute_metrics
grep -rn "compute_metrics(" --include="*.R" | grep -v "function" | grep -v "_bundle" | grep -v "_CV" | grep -v "#"
# Result: Only the commented-out call in _targets.R
```

## Conclusion

### Impact: NONE on current pipeline results

The `compute_metrics.R` function with the bug is **dead code** in the current pipeline:
- Not called by simulation targets
- Not called by validation targets
- `compute_metrics_bundle.R` exists but is never invoked

### Recommendation

1. **No need to rerun results** - the bug does not affect any computed outputs
2. **Should still fix the bug** - prevent future issues if this code is reactivated
3. **Consider removing** `compute_metrics.R` and `compute_metrics_bundle.R` if confirmed unused

### Risk Assessment
| Pipeline Component | Uses compute_metrics.R? | Impact |
|-------------------|------------------------|--------|
| Simulation C-index | No | None |
| Validation C-index | No | None |
| Training metrics | No (commented out) | None |
| Paper figures | No | None |

---

## Resolution (2026-01-25)

Bug was fixed and dead code removed:
- Fixed `sdZ = fit$meanZ[,ns]` → `sdZ = fit$sdZ[,ns]` in `R/compute_metrics.R`
- Removed `R/compute_metrics.R` (unused)
- Removed `R/compute_metrics_bundle.R` (unused)
- Removed `R/compute_metrics_CV.R` (unused)

These files were confirmed to have no active callers in the pipeline.

---
*Analysis performed: 2026-01-25*
*Pipeline store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125*
