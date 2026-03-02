# Figure 4 Caption-Code Mismatch

**Date:** 2026-02-23
**Priority:** Low (code is correct; caption is slightly inaccurate)
**File:** `paper/04_results.Rmd`, line 122 (fig-bio caption)

## Issue

Caption says:
> "quantified by the change in partial log-likelihood from univariate Cox models"

Code (`R/figure_targets.R:2335-2359`, `compute_survival_explained()`) implements Type III leave-one-out:
```r
full_fit <- coxph(Surv(time, event) ~ XtW)  # all k factors
ll_full <- full_fit$loglik[2]
# For each factor j:
ll_full - reduced_fit$loglik[2]  # full minus (k-1)-factor model
```

## Suggested Fix

Replace "from univariate Cox models" with "when each factor is removed from the full multivariable Cox model (Type III)".

## Notes

- The code is correct and the Type III approach is stronger than univariate
- This strengthens the "concentration" claim (conditional, not marginal)
- No results change needed, only the caption text
