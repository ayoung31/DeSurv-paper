# Plan: CPTAC-Integrated Bayesian Optimization

## Motivation

The current pipeline selects hyperparameters (k, alpha, lambda, nu) by maximizing
5-fold CV C-index on the training data (TCGA or TCGA+CPTAC). Even with the 1SE rule
applied to k, alpha, and lambda, this creates selection bias: the reported CV C-index
is optimistically biased because the parameters were chosen to maximize it on the same
data used to evaluate it.

The goal is to use CPTAC C-index as the **BO objective itself**, so parameter selection
is driven entirely by held-out performance. The remaining datasets (Dijk, Moffitt,
PACA-AU, Puleo) serve as a completely untouched external test set.

---

## Evaluation Chain

| Stage | Data |
|---|---|
| BO training (model fitting) | TCGA |
| BO objective (parameter selection) | CPTAC C-index |
| External test (completely held out) | Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, Puleo_array |

---

## Core Design

For each BO candidate (k, alpha, lambda, nu):
1. Fit `desurv_fit` on all TCGA with `n_starts` random initializations, keep best
2. Predict on CPTAC
3. Compute CPTAC C-index → this is the value the GP models

No internal CV on TCGA. The GP surrogate smooths over the noise in CPTAC C-index
estimates across candidate evaluations. The 1SE rule for model selection uses GP
posterior SD (`preds$sd` from `extract_gp_curve`) as the uncertainty estimate,
which already happens in the existing selection logic.

---

## Changes Required

### 1. DeSurv package — external validation objective

Modify `desurv_cv_bayesopt` (and by extension `desurv_cv_bayesopt_refine`) to accept:

```r
X_val = NULL, y_val = NULL, d_val = NULL
```

When provided, each BO evaluation becomes:
- Fit `desurv_fit` on `(X, y, d)` with `n_starts` random inits, keep best by training likelihood
- Predict on `X_val`, compute C-index against `(y_val, d_val)`
- Return that C-index as the BO objective (replacing CV mean)

The internal CV path (`desurv_cv`) is bypassed entirely when `X_val` is provided.

Relevant files in `../DeSurv/R/`:
- `desurv_cv_bayesopt.R` — main BO loop, add val data path in `eval_point()`
- `desurv_bo_refine.R` — pass `X_val/y_val/d_val` through `common_args`

### 2. DeSurv package — SE for 1SE rule

Currently `c_se` comes from fold-to-fold variance in CV. With CPTAC as the objective
there are no folds. Use GP posterior SD (`preds$sd`) as the uncertainty, which is
already the basis for `lower`/`upper` in `extract_gp_curve`. No code change needed
in the selection logic — the 1SE rule via `select_1se_from_gp_curve` already uses
`best$lower = best$mean - z_value * best$sd`.

### 3. Pipeline — new BO config

Add a new entry `tcga_cptacval` to `targets_bo_configs.R`:

```r
tcga_cptacval = list(
  data_mode = "external",
  data_loader = "load_data",
  train_datasets = c("TCGA_PAAD"),
  val_dataset = "CPTAC",           # new field — signals CPTAC-as-BO-objective
  method_trans_train = "rank",
  desurv_bo_bounds = list(
    k_grid = list(lower = 2L, upper = 12L, type = "integer"),
    alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
    lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
    nu_grid = list(lower = 0, upper = 1, type = "continuous")
  ),
  ngene_config = c(3000),
  ntop_config = c(3000),           # all genes (fixed, not tuned)
  lambdaw_config = c(0),
  lambdah_config = c(0),
  ninit = ???,                     # see open questions
  bo_n_init = 50,
  bo_n_iter = ???,                 # see open questions
  ...
)
```

New targets needed in `targets_common_pipeline.R`:
- `tar_cptac_val_data` — load and preprocess CPTAC, gene-space aligned to TCGA
- Conditional branch in `desurv_bo_results` target: when `bo_config$val_dataset`
  is set, pass `X_val/y_val/d_val` to `desurv_cv_bayesopt_refine`

### 4. Pipeline — validation config

For `tcga_cptacval`, the val config must exclude CPTAC (it was consumed by BO):

```r
# in targets_val_configs.R, add a branch for tcga_cptacval:
val_datasets = c(
  "Dijk",
  "Moffitt_GEO_array",
  "PACA_AU_array",
  "PACA_AU_seq",
  "Puleo_array"
)
```

---

## Key Differences from Two-Phase Approach

The previously considered two-phase approach (BO on TCGA CV, then post-hoc CPTAC
evaluation for parameter selection) is replaced here. The integrated approach is
preferred because:
- CPTAC signal directly shapes the BO search trajectory, not just final selection
- No risk of the GP curve misrepresenting CPTAC performance (it IS the GP objective)

| | Two-phase | CPTAC-integrated BO |
|---|---|---|
| BO objective | TCGA 5-fold CV C-index | CPTAC C-index |
| SE for 1SE rule | Fold variance | GP posterior SD |
| Noise level | Low (averaged over 5 folds) | Higher (~140 CPTAC samples) |
| BO iterations needed | ~150 (current) | Likely more (noisier objective) |
| Cost per BO iteration | 5 × `desurv_cv` folds | 1 × `desurv_fit` (cheaper per iter) |
| Risk | Selection bias from TCGA CV | Noisier GP convergence |

---

## Open Questions

### Q1: Is CPTAC n large enough for a stable BO objective?

CPTAC has ~140 samples. C-index SE at that n is roughly 0.04–0.06. This means
adjacent BO candidates with similar true performance will be hard to distinguish.
The GP can model this noise, but convergence may be slow or the selected parameters
may vary substantially across runs.

**Decision needed:** Is the noise acceptable, or should we use a combined objective
(see Q3)?

### Q2: n_starts (ninit) per BO evaluation

Currently `ninit = 30` random starts per CV evaluation. With CPTAC as the objective
(fitting on all TCGA, predicting on CPTAC), each BO eval is a full `desurv_fit`
with 30 inits. This is computationally similar to the current cost but without the
5-fold multiplier — so overall BO should be cheaper per iteration.

**Decision needed:** Keep `ninit = 30`, or reduce to 10–15 since we're not averaging
over folds and the noise floor is already set by CPTAC sample size?

### Q3: Pure CPTAC C-index vs. combined objective

Two options for the BO objective:
- **Pure CPTAC:** `objective = CPTAC_cindex` — cleanest, fully external
- **Combined:** `objective = 0.5 * TCGA_CV_cindex + 0.5 * CPTAC_cindex` — more
  stable signal, but reintroduces partial selection bias on TCGA

**Decision needed:** Which do you prefer? Pure CPTAC is more defensible scientifically.
Combined is more statistically stable.

### Q4: How to handle gene-space alignment for CPTAC val data

TCGA and CPTAC may not share all genes. The BO fitting happens on TCGA gene space
(after filtering to `ngene = 3000`). The CPTAC prediction target must be restricted
to those same genes.

Two options:
- Compute gene filter on TCGA, then intersect with CPTAC before BO starts
  (cleaner — done once in a pre-processing target)
- Let `predict.desurv_fit` handle the intersection at prediction time

**Decision needed:** Pre-processing target approach is preferred to avoid repeated
intersection overhead inside the BO loop.

---

## Files to Modify

| File | Change |
|---|---|
| `../DeSurv/R/desurv_cv_bayesopt.R` | Add `X_val/y_val/d_val` params, val C-index path in `eval_point()` |
| `../DeSurv/R/desurv_bo_refine.R` | Pass val data through `common_args` |
| `targets_bo_configs.R` | Add `tcga_cptacval` config with `val_dataset` field |
| `targets_val_configs.R` | Exclude CPTAC from val datasets for `tcga_cptacval` |
| `targets_common_pipeline.R` | Add `tar_cptac_val_data` target; conditional val-data branch in BO target |
| `R/targets_config.R` | Add `val_dataset` to config schema and validation |
