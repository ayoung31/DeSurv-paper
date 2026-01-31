# Changelog

All notable changes to the DeSurv paper pipeline on the `naimedits0125` branch.

## [2026-01-31] - Today

### Pipeline Configuration
- **Updated tcgacptac BO config to match student's original settings:**
  - `k_grid`: 2-15 → 2-12 (matching student's upper bound)
  - `alpha_grid`: 0-0.95 → 0-1.0 (allowing full supervision)
  - `ngene_config`: 2000 → 3000
  - `ntop_config`: c(100, 200) → c(50, 300)
  - `bo_n_init`: 20 → 50
  - `bo_n_iter`: 50 → 100
- **Enabled parallel processing for bladder config:**
  - `desurv_parallel_grid`: FALSE → TRUE
  - `desurv_ncores_grid`: 5 (one core per CV fold)

### Bug Fixes
- Fixed store path mismatch between `paper/paper.Rmd` and `paper/_targets.yaml`
- Added missing packages to `TARGET_PACKAGES` in `targets_setup.R`:
  - `ggrepel` (for `fig_variation_explained_*` targets)
  - `survminer` (for `fig_median_survival_*` targets)

### Simulation Pipeline
- **Diagnosed missing NMF comparison in Figure 3:**
  - Root cause: `bo_alpha0` analysis spec was never executed
  - Stored results only had `fixed` and `bo` (2 of 6 specs)
  - Only 2 scenarios with 20 replicates (should be 4 scenarios × 100 replicates)
- **Submitted full simulation reproduction (Job 267):**
  - 2,400 runs: 6 analysis specs × 4 scenarios × 100 replicates
  - Invalidated stale simulation targets in `_full` store
  - Created `_targets_sims_local.sh` for local desktop execution

### Store Alignment
- Aligned all configurations to use `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full`:
  - `paper/paper.Rmd` params
  - `paper/_targets.yaml`
  - `_targets_sims_local.sh`

### Jobs Submitted
- Job 259: Main BO with student's config (tcgacptac + bladder)
- Job 267: Full simulation reproduction (queued, ~60-80 hours estimated)

---

## [2026-01-30] - Yesterday

### Pipeline Updates
- Multiple configuration updates for main analysis pipeline
- Figure regeneration with updated parameters

---

## [2026-01-29]

### Pipeline Maintenance
- Continued pipeline execution and monitoring
- Configuration adjustments

---

## [2026-01-28]

### Pipeline Updates
- Multiple incremental updates to pipeline configuration
- Bug fixes and parameter tuning

---

## [2026-01-27]

### Updates
- Pipeline configuration refinements

---

## [2026-01-26]

### Updates
- Pipeline execution and result collection
- Configuration updates

---

## [2026-01-25]

### Major Setup
- **Local Desktop Pipeline Configuration:**
  - Configured local Slurm cluster for job scheduling
  - Set up `crew_controller_local` for parallel execution
  - Created quick-mode test configuration

- **Bug Fixes (from CODE_REVIEW.md):**
  - Bug 1: Fixed `compute_metrics.R:9` - removed dead code file
  - Bug 2: Fixed `select_best_init.R:5` - removed dead code file
  - Bug 3: Removed all 6 `browser()` debug statements from `targets_common_pipeline.R`
  - Bug 4: Analyzed `src/functions.cpp:472` - confirmed NOT a bug (correct maximization logic)

- **Documentation:**
  - Added comprehensive `MANUSCRIPT_PIPELINE_ANALYSIS.md`
  - Updated `CODE_REVIEW.md` with resolution status

- **Configuration:**
  - Updated config files to reference `naimedits0125` branch store
  - Set `DEFAULT_NINIT=19`, `DEFAULT_NINIT_FULL=19` for local desktop

---

## [2026-01-24]

### Initial Branch Setup
- Created `naimedits0125` branch from student's work
- Initial configuration adjustments for local development

---

## Pending Work

### Figures Requiring Updates
| Figure | Status | Blocker |
|--------|--------|---------|
| Fig 2D | ⏳ Pending | Main BO (Job 259) completion |
| Fig 2E | ⏳ Pending | Simulation results (Job 267) |
| Fig 3A-B | ⏳ Pending | Simulation results (Job 267) |

### Expected Timeline
- Main BO completion: ~hours (Job 259 running)
- Simulation completion: ~60-80 hours (Job 267 queued)
- Paper render: After all jobs complete
