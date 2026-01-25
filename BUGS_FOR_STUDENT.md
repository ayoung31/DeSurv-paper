# Bug Fixes for Review

This document explains 4 bugs found during code review. Please fix each one and add a test to prevent regression.

---

## Bug 1: Variable Assignment Error in compute_metrics.R

**File:** `R/compute_metrics.R`
**Lines:** 8-9

### What's Wrong

Look at lines 8 and 9:
```r
meanZ = fit$meanZ[,ns]
sdZ = fit$meanZ[,ns]    # <-- BUG: This should be fit$sdZ
```

Line 9 assigns `fit$meanZ` to `sdZ`, but it should assign `fit$sdZ`. This is a classic copy-paste error where line 9 was copied from line 8 but the variable name wasn't updated.

### Why It Matters

The variable `sdZ` (standard deviation) is used for two critical operations:

1. **Line 12 - Standardization:** `XtW = sweep(XtW, 2, sdZ, FUN="/")` divides by standard deviation
2. **Line 13 - Linear predictor:** `lp <- XtW %*% (beta * sdZ)` multiplies by standard deviation

When you use `meanZ` instead of `sdZ`, you're:
- Dividing by the mean instead of the standard deviation
- Scaling beta by the mean instead of the standard deviation

This produces mathematically incorrect C-index values and all downstream metrics.

### How to Fix

Change line 9 from:
```r
sdZ = fit$meanZ[,ns]
```
to:
```r
sdZ = fit$sdZ[,ns]
```

### What to Learn

**Always review adjacent lines with similar patterns.** Copy-paste errors are extremely common when you have parallel operations like extracting `meanZ` and `sdZ`. A good practice is to double-check variable names on both sides of assignments when copying code.

---

## Bug 2: Unused Filter Variable in select_best_init.R

**File:** `R/select_best_init.R`
**Lines:** 5-11

### What's Wrong

Look at the code flow:
```r
select_best_init <- function(df, method_select="surv") {

  # Line 5: Creates a filtered dataframe
  keep <- df[!df$flag_nan & is.finite(df$loss), ]

  # Line 8: Modifies the ORIGINAL df, not keep
  df$method_select = method_select

  # Lines 10-11: Uses ORIGINAL df, not keep
  if(method_select=="surv"){
    best = df %>% dplyr::slice_max(order_by = sloss, n = 1, with_ties = FALSE)
  }
  # ...
}
```

The filtered dataframe `keep` is created but **never used**. The function continues to operate on the original `df` which may contain:
- Rows where `flag_nan == TRUE` (NaN errors occurred)
- Rows where `loss` is infinite

### Why It Matters

The purpose of line 5 is to exclude bad initializations before selecting the "best" one. But since `keep` is never used:

1. An initialization with `flag_nan = TRUE` could be selected as "best"
2. An initialization with `loss = Inf` could be selected as "best"
3. The model selection is unreliable

This could cause downstream failures or select a model that didn't converge properly.

### How to Fix

Replace all uses of `df` after line 5 with `keep`:

```r
select_best_init <- function(df, method_select="surv") {

  # only keep initializations with finite loss and no errors
  keep <- df[!df$flag_nan & is.finite(df$loss), ]

  # Check if any valid initializations remain
  if (nrow(keep) == 0) {
    stop("No valid initializations found (all had NaN or infinite loss)")
  }

  # save method for init selection in dataframe
  keep$method_select = method_select

  if(method_select=="surv"){
    best = keep %>% dplyr::slice_max(order_by = sloss, n = 1, with_ties = FALSE)
  }else if(method_select=="nmf"){
    best = keep %>% dplyr::slice_min(order_by = nloss, n = 1, with_ties = FALSE)
  }else if(method_select=="loss"){
    best = keep %>% dplyr::slice_min(order_by = loss, n = 1, with_ties = FALSE)
  }else{
    stop("This method for selecting the best initialization is not supported")
  }

  return(best)
}
```

### What to Learn

**Trace your variables through the function.** When you create a filtered/transformed variable, make sure you actually use it. This bug pattern is common when:
- You add defensive filtering as an afterthought
- You refactor code and forget to update all references

A good practice is to search for all uses of the original variable name after creating a filtered version.

---

## Bug 3: Debug Statements Left in Production Code

**File:** `targets_common_pipeline.R`
**Lines:** 1691, 1803, 2075, 2097, 2119, 2141

### What's Wrong

There are 6 `browser()` statements scattered through the pipeline code:

```r
# Line 1691
browser()

# Line 1803
browser()

# Line 2075
browser()

# ... etc
```

### Why It Matters

`browser()` is R's interactive debugger. When the pipeline hits one of these lines, it will:

1. **Stop execution completely**
2. **Wait for interactive input** from a user
3. **Never continue** in a batch/Slurm environment

Since this code runs on Slurm (non-interactive), the pipeline will hang indefinitely at these points, wasting compute resources and never completing.

### How to Fix

Remove all 6 `browser()` calls. You can find them with:
```bash
grep -n "browser()" targets_common_pipeline.R
```

Then delete each line or comment them out.

### What to Learn

**Never commit debug statements to production code.** Best practices:

1. Use `git diff` or `git status` before committing to review changes
2. Add a pre-commit hook that rejects files containing `browser()`
3. If you must keep debug code, use a conditional:
   ```r
   if (interactive() && getOption("debug_mode", FALSE)) {
     browser()
   }
   ```

---

## Bug 4: Beta Backtracking Logic (DeSurv Package) - NEEDS VERIFICATION

**File:** `../DeSurv/src/functions.cpp`
**Lines:** 472-495

### Context

This is in the C++ optimizer code. The `penalized_surv_loss` function returns:
```cpp
return cs.loglik * 2.0 / n_event - lambda * penalty_beta;
```

Since `loglik` is a Cox partial log-likelihood, **higher values are better** (we want to maximize this).

### The Code in Question

```cpp
double surv_loss_candidate = penalized_surv_loss(Z, y, delta, beta_full, ...);

if (!std::isfinite(surv_loss_candidate) ||
    surv_loss_candidate < surv_loss_prev) {
  // Enter backtracking...
  for (int bt = 0; bt < max_backtracks_beta; ++bt) {
    // ...
    if (std::isfinite(surv_loss_trial) &&
        surv_loss_trial >= surv_loss_prev) {
      beta_full = beta_trial;
      accepted = true;
      break;
    }
  }
}
```

### Your Task

Analyze this code carefully and determine:

1. What does `surv_loss_candidate < surv_loss_prev` mean in context? (candidate is worse/better?)
2. When should backtracking be triggered?
3. Is the condition correct for a quantity we want to maximize?

**Hint:** Think about what it means to "backtrack" - you backtrack when your step made things worse, not better.

Write up your analysis with:
- Your conclusion (bug or not?)
- If it's a bug, what the fix should be
- If it's not a bug, explain why the logic is correct

---

## Testing Your Fixes

After fixing each bug, add a test in `tests/testthat/` to prevent regression:

### For Bug 1 (compute_metrics.R):
```r
test_that("compute_metrics uses sdZ not meanZ", {
  # Create a mock fit object where meanZ != sdZ
  mock_fit <- list(
    W = matrix(1, nrow = 10, ncol = 2),
    beta = c(0.5, 0.5),
    meanZ = matrix(c(5, 5), nrow = 1),  # Different from sdZ
    sdZ = matrix(c(1, 1), nrow = 1),    # Standard deviation = 1
    convergence = TRUE,
    iter = 100,
    nan_flag = FALSE,
    loss = list(loss = 0.5, surv_loss = 0.3, nmf_loss = 0.2)
  )
  # ... test that sdZ is actually used
})
```

### For Bug 2 (select_best_init.R):
```r
test_that("select_best_init excludes NaN and infinite loss rows", {
  df <- data.frame(
    flag_nan = c(FALSE, TRUE, FALSE),
    loss = c(0.5, 0.3, Inf),
    sloss = c(0.8, 0.9, 0.7)  # The NaN row has highest sloss
  )
  result <- select_best_init(df, method_select = "surv")

  # Should select row 1, NOT row 2 (which has highest sloss but flag_nan=TRUE)
  expect_equal(nrow(result), 1)
  expect_false(result$flag_nan)
  expect_true(is.finite(result$loss))
})
```

---

## Summary Checklist

- [x] Fix `sdZ` assignment in `R/compute_metrics.R:9` - **RESOLVED**: File removed (dead code, see BUG1_IMPACT_ANALYSIS.md)
- [x] Use `keep` variable in `R/select_best_init.R` - **RESOLVED**: File removed (dead code, never called in pipeline)
- [x] Remove all `browser()` from `targets_common_pipeline.R` - **RESOLVED**: All 6 browser() statements removed
- [x] Analyze and document C++ backtracking logic (Bug 4) - **RESOLVED**: NOT A BUG - logic is correct (see BUG4_ANALYSIS.md)
- [ ] Add regression tests for each fix
- [ ] Run full test suite: `Rscript -e 'testthat::test_dir("tests/testthat")'`
