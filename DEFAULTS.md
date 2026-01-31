# Parameter Defaults Registry

Single source of truth for all configurable parameters. When changing defaults, update both the code AND this file.

**Last updated**: 2026-01-31
**Branch**: `naimedits0125`

---

## Bayesian Optimization (BO) Config

Source: `targets_bo_configs.R`

### tcgacptac (Main PDAC Analysis)

| Parameter | Current Value | Student Original | Notes |
|-----------|---------------|------------------|-------|
| `k_grid` | 2-12 | 2-12 | Number of latent factors |
| `alpha_grid` | 0-1.0 | 0-1.0 | Supervision strength |
| `lambda_grid` | 1e-3 to 1e3 (log10) | same | Cox regularization |
| `nu_grid` | 0-1 | 0-1 | NMF regularization |
| `ngene_config` | 3000 | 3000 | Top variable genes |
| `ntop_config` | c(50, 300) | c(50, 300) | Top genes per factor |
| `bo_n_init` | 50 | 50 | Initial random BO samples |
| `bo_n_iter` | 100 | 100 | BO optimization iterations |
| `bo_candidate_pool` | 4000 | 4000 | Candidates per iteration |
| `ninit` | 30 | 30 | NMF random initializations |
| `nfold` | 5 | 5 | CV folds |
| `desurv_ncores_grid` | 5 | varies | Cores per CV fold |

### bladder (Bladder Cancer Analysis)

| Parameter | Current Value | Notes |
|-----------|---------------|-------|
| `k_grid` | 2-10 | Smaller upper bound |
| `alpha_grid` | 0-0.95 | Slightly constrained |
| `ngene_config` | 2000 | Fewer genes (smaller dataset) |
| `ntop_config` | c(100, 200) | Different range |
| `bo_n_init` | 20 | Fewer iterations |
| `bo_n_iter` | 50 | Fewer iterations |
| `ninit` | 30 | Same as tcgacptac |

---

## Simulation Config

Source: `_targets_sims.R`

| Parameter | Value | Notes |
|-----------|-------|-------|
| `SIM_DATASETS_PER_SCENARIO` | 100 | Replicates per scenario |
| `SIM_CV_NSTARTS` | 30 | NMF initializations in sims |
| `SIM_CV_NFOLDS` | 5 | CV folds |
| `SIM_BO_N_INIT` | 20 | BO initial samples |
| `SIM_BO_N_ITER` | 40 | BO iterations |
| `SIM_BO_CANDIDATE_POOL` | 1500 | Candidates |
| `SIM_DEFAULT_K` | 3 | True k in simulations |
| `SIM_DEFAULT_ALPHA` | 0.6 | Default supervision |
| `SIM_TRAIN_FRACTION` | 0.7 | Train/test split |

### Simulation Scenarios

| Scenario ID | Description | True k | Beta |
|-------------|-------------|--------|------|
| `R0_easy` | Default easy | 3 | default |
| `R0k6` | Higher complexity | 6 | c(2.0,-1.5,1.0,0,0,0) |
| `R00_null` | No survival signal | 3 | all zeros |
| `R_mixed` | Mixed markers | 3 | default |

### Analysis Specs

| Analysis ID | Mode | Alpha | Purpose |
|-------------|------|-------|---------|
| `fixed` | fixed | 0.6 | DeSurv with fixed params |
| `fixed_alpha0` | fixed | 0 | NMF baseline |
| `bo` | bayesopt | tuned | DeSurv with BO |
| `bo_alpha0` | bayesopt | 0 | NMF with BO (k selection) |
| `bo_tune_ntop` | bayesopt | tuned | BO including ntop |
| `bo_tune_ntop_alpha0` | bayesopt | 0 | NMF BO with ntop |

---

## Run Config

Source: `targets_run_configs.R`

| Parameter | Value | Notes |
|-----------|-------|-------|
| `ninit_full` | 19 | Full model initializations |
| References | `bo_key` | Links to BO config |

---

## Validation Config

Source: `targets_val_configs.R`

### tcgacptac Validation Datasets

| Dataset | Mode | Notes |
|---------|------|-------|
| CPTAC | external | Held-out during training |
| Dijk | external | External cohort |
| Moffitt_GEO_array | external | External cohort |
| PACA_AU_array | external | External cohort |
| PACA_AU_seq | external | External cohort |
| Puleo_array | external | External cohort |

### bladder Validation

| Dataset | Mode | Notes |
|---------|------|-------|
| train_split | train_split | 30% held out |

---

## Controller Config (Local Desktop)

Source: `targets_setup.R`, `_targets_sims.R`

| Controller | Workers | Purpose |
|------------|---------|---------|
| `crew_controller_local` | 2 | Main pipeline |
| `sim_analysis_controller` | 2 | Simulations |

---

## Updating This File

When changing any default:

1. Update the code file
2. Update this file with new value
3. Add entry to CHANGELOG.md
4. Run `scripts/verify_consistency.R`
