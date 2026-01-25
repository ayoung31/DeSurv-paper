# Bug 4 Analysis: Beta Backtracking Logic in DeSurv C++

## Status: NOT A BUG - Code is Correct

**File:** `../DeSurv/src/functions.cpp`
**Lines:** 472-495

## Code Under Review

```cpp
double surv_loss_candidate = penalized_surv_loss(Z, y, delta, beta_full,
                                                  n_event, lambda, nu, sdZ);

if (!std::isfinite(surv_loss_candidate) ||
    surv_loss_candidate < surv_loss_prev) {
  arma::vec direction = beta_full - beta0;
  double theta = theta_beta_init * rho_beta;
  bool accepted = false;

  for (int bt = 0; bt < max_backtracks_beta; ++bt) {
    arma::vec beta_trial = beta0 + theta * direction;
    double surv_loss_trial = penalized_surv_loss(
      Z, y, delta, beta_trial, n_event, lambda, nu, sdZ);
    if (std::isfinite(surv_loss_trial) &&
        surv_loss_trial >= surv_loss_prev) {
      beta_full = beta_trial;
      accepted = true;
      break;
    }
    theta *= rho_beta;
    if (theta < 1e-8) break;
  }

  if (!accepted) {
    beta_full = beta0;
  }
}
```

## Analysis

### What `penalized_surv_loss` Returns

```cpp
return cs.loglik * 2.0 / n_event - lambda * penalty_beta;
```

- `cs.loglik` is the Cox partial log-likelihood
- The Cox log-likelihood is a quantity we want to **MAXIMIZE** (higher = better fit)
- The penalty term reduces the score for larger coefficients (L1/L2 regularization)
- **Conclusion:** Higher `surv_loss` values are better

### Backtracking Trigger Condition (Line 472-473)

```cpp
if (!std::isfinite(surv_loss_candidate) ||
    surv_loss_candidate < surv_loss_prev)
```

- **Non-finite:** Always backtrack if we got NaN/Inf
- **`candidate < prev`:** Since higher is better, `candidate < prev` means the candidate is **WORSE**
- **This triggers backtracking when the full step made things worse**
- **CORRECT BEHAVIOR**

### Trial Acceptance Condition (Line 482-483)

```cpp
if (std::isfinite(surv_loss_trial) &&
    surv_loss_trial >= surv_loss_prev)
```

- **Finite check:** Must have a valid number
- **`trial >= prev`:** Accept if trial is **at least as good** as the previous value
- **This is the standard Armijo sufficient decrease condition for MAXIMIZATION**
- **CORRECT BEHAVIOR**

### Algorithm Flow

1. Take a full Newton/gradient step to get `beta_full`
2. Evaluate the objective at `beta_full`
3. If objective got worse (or is non-finite):
   - Start backtracking: try smaller steps in the same direction
   - `theta` starts at 0.9 and shrinks by factor of 0.9 each iteration
   - Accept the first step that achieves `trial >= prev`
4. If no acceptable step found after 30 iterations (or theta < 1e-8):
   - Revert to `beta0` (the starting point)
5. Return the best `beta_full`

This is a standard **Armijo backtracking line search** adapted for maximization.

## Why the Initial Assessment Was Wrong

The bug report suggested the condition might be "inverted." This confusion likely arose because:

1. Most textbooks present backtracking for **minimization** where `candidate > prev` triggers backtracking
2. Cox models **maximize** log-likelihood, so the inequality flips
3. The code correctly implements backtracking for a **maximization** problem

## Verification

The condition is consistent with:
- Standard Armijo line search for maximization
- The Wolfe conditions adapted for maximization
- Common implementations in optimization libraries (e.g., scipy's line search for maximization)

## Conclusion

**No code change needed.** The backtracking logic is mathematically correct for maximizing the penalized Cox partial log-likelihood.

| Condition | Meaning | Correct? |
|-----------|---------|----------|
| `candidate < prev` triggers backtracking | Full step made things worse | Yes |
| `trial >= prev` accepts step | Trial is at least as good | Yes |

---

*Analysis performed: 2026-01-25*
