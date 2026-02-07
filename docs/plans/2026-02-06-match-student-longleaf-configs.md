# Match Student Longleaf Configs Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make `local_slurm/longleaf/` configs produce identical BO parameters, controllers, and validation settings to the student's `master` branch, so that running on Longleaf reproduces the student's results exactly.

**Architecture:** The student's `master` uses a monolithic `targets_setup.R` + `_targets.R` where all parameters are globals. Our branch uses modular config functions (`targets_bo_configs()`, etc.) that feed into `resolve_desurv_bo_config()` in `R/targets_config.R`, which merges user values onto defaults and auto-generates `bo_coarse_control`/`bo_refine_control` when NULL. We only need to change the 4 config files in `local_slurm/longleaf/` plus `targets_setup.R` — no changes to pipeline code or `R/` helpers.

**Tech Stack:** R, targets, crew.cluster, Slurm

---

## Summary of All Differences

Student's `master` _targets.R values (the source of truth) are shown in the "Student" column. These are the values set in `_targets.R` (which override `targets_setup.R` defaults via assignment before `source("targets_setup.R")`).

### BO Parameters (TCGA+CPTAC only — student has no bladder config)

| Parameter | Student (`master`) | Current Longleaf | Action |
|---|---|---|---|
| k_grid upper | **10L** | 12L | Change to 10L |
| ngene_config | **c(500, 5000)** | c(3000) | Change to c(500, 5000) |
| ntop_config | **c(50, 250)** | c(50, 300) | Change to c(50, 250) |
| lambdah_config | **c(1e-7, 1e2)** | c(0) | Change to c(1e-7, 1e2) |
| bo_n_init | **20** | 50 | Change to 20 |
| bo_n_iter | **100** | 100 | OK |
| bo_candidate_pool | **4000** | 4000 | OK |
| bo_max_refinements | **2** | 0 | Change to 2 |
| bo_shrink_base | **0.5** | 0.3 | Change to 0.5 |
| bo_importance_gain | **0.3** | 0.1 | Change to 0.3 |
| ninit | **50** | 50 | OK |
| desurv_ncores_grid | **50** (= NINIT) | 50 | OK |
| desurv_parallel_grid | **TRUE** | TRUE | OK |
| nfold | **5** | 5 | OK |
| bo_tol (TOL) | **1e-5** | 1e-5 | OK |
| bo_maxit (MAXIT) | **4000** | 4000 | OK |
| bo_tol_gain | **0.002** | 0.002 | OK |
| bo_plateau | **1** | 1 | OK |
| bo_top_k | **10** | 10 | OK |

### Bladder Config

| Parameter | Student | Current Longleaf | Action |
|---|---|---|---|
| (entire config) | **does not exist** | exists | Remove bladder config |

### Validation Datasets

| Parameter | Student | Current Longleaf | Action |
|---|---|---|---|
| val_datasets | **Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, Puleo_array** (5) | adds CPTAC (6) | Remove CPTAC from val list |

### Controllers (targets_setup.R)

| Controller | Student | Current Longleaf | Action |
|---|---|---|---|
| default | **sequential** | 50 Slurm workers | Change to sequential |
| low_mem | **202 workers, 1 GB/cpu, 120 min** | absent | Add |
| cv | **202 workers, 2 GB/cpu, NINIT CPUs, 720 min** | same | OK |
| full | **202 workers, 32 GB required, NINIT_FULL CPUs, 600 min** | absent | Add |
| nmf | absent | 10 workers, 10 CPUs | Remove |
| med_mem | **120 workers, 8 GB, 200 min** | 120 workers, 8 GB, 240 min | Change to 200 min |
| controller group | default, low_mem, cv, full, med_mem | default, cv, nmf, med_mem | Match student's group |

### TARGET_PACKAGES

| Student | Current Longleaf | Action |
|---|---|---|
| missing `digest`, `caret`, `ggrepel`, `survminer` | has all four | Remove the 4 extra packages |

### Other

| Parameter | Student | Current Longleaf | Action |
|---|---|---|---|
| `LOCAL_RENDER` env var | **present** (guards crew.cluster import) | absent | Add LOCAL_RENDER guard |

---

## Tasks

### Task 1: Update `targets_bo_configs.R` — match student's BO parameters

**Files:**
- Modify: `local_slurm/longleaf/targets_bo_configs.R`

**Step 1: Replace file contents**

The student only has TCGA+CPTAC (no bladder). Match every parameter to `master:_targets.R`:

```r
# targets_bo_configs.R - UNC Longleaf HPC mode
# Matched to student's master branch for exact reproduction.

targets_bo_configs <- function() {
  list(
    tcgacptac = list(
      data_mode = "external",
      data_loader = "load_data",
      train_datasets = c("TCGA_PAAD", "CPTAC"),
      method_trans_train = "rank",
      desurv_bo_bounds = list(
        k_grid = list(lower = 2L, upper = 10L, type = "integer"),
        alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
        lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
        nu_grid = list(lower = 0, upper = 1, type = "continuous")
      ),
      ngene_config = c(500, 5000),
      ntop_config = c(50, 250),
      lambdaw_config = c(0),
      lambdah_config = c(1e-7, 1e2),
      ninit = 50,
      bo_n_init = 20,
      bo_n_iter = 100,
      bo_candidate_pool = 4000,
      bo_max_refinements = 2,
      bo_tol_gain = 0.002,
      bo_plateau = 1,
      bo_top_k = 10,
      bo_shrink_base = 0.5,
      bo_importance_gain = 0.3,
      bo_coarse_control = NULL,
      bo_refine_control = NULL,
      bo_tol = 1e-5,
      bo_maxit = 4000,
      nfold = 5,
      desurv_parallel_grid = TRUE,
      desurv_ncores_grid = 50
    )
  )
}
```

**Step 2: Verify auto-generated coarse/refine controls match student**

With `bo_coarse_control = NULL`, `resolve_desurv_bo_config()` in `R/targets_config.R:140-148` will generate:
```r
list(n_init = 20, n_iter = 100, candidate_pool = 4000,
     exploration_weight = 0.01, seed = 123, cv_verbose = FALSE)
```
This exactly matches student's `BO_COARSE_CONTROL`. Similarly for refine (seed = 456). No changes needed.

**Step 3: Commit**

```bash
git add local_slurm/longleaf/targets_bo_configs.R
git commit -m "match student BO params: ngene, lambdah, ntop, k_grid, refinement"
```

---

### Task 2: Update `targets_val_configs.R` — remove CPTAC from validation

**Files:**
- Modify: `local_slurm/longleaf/targets_val_configs.R`

**Step 1: Replace file contents**

Student's validation list does NOT include CPTAC (it's a training dataset):

```r
# targets_val_configs.R - UNC Longleaf HPC mode
# Matched to student's master branch for exact reproduction.

targets_val_config <- function(label) {
  list(
    mode = "external",
    val_datasets = c(
      "Dijk",
      "Moffitt_GEO_array",
      "PACA_AU_array",
      "PACA_AU_seq",
      "Puleo_array"
    )
  )
}

targets_val_configs <- function(bo_labels = names(targets_bo_configs())) {
  if (is.null(bo_labels)) {
    bo_labels <- character(0)
  }
  configs <- lapply(bo_labels, targets_val_config)
  names(configs) <- bo_labels
  configs
}
```

**Step 2: Commit**

```bash
git add local_slurm/longleaf/targets_val_configs.R
git commit -m "remove CPTAC from validation to match student"
```

---

### Task 3: Update `targets_setup.R` — match student's controllers and packages

**Files:**
- Modify: `local_slurm/longleaf/targets_setup.R`

**Step 1: Replace file contents**

Key changes:
- `default_controller` → sequential (not Slurm)
- Add `low_mem_controller` (202 workers, 1 GB/cpu, 120 min)
- Add `full_run_controller` (202 workers, 32 GB, NINIT_FULL CPUs, 600 min)
- Remove `nmf_controller`
- `med_mem_controller` → 200 min (was 240)
- Add `LOCAL_RENDER` guard around `crew.cluster`
- Remove `digest`, `caret`, `ggrepel`, `survminer` from TARGET_PACKAGES
- Controller group: default, low_mem, cv, full, med_mem

```r
# targets_setup.R - UNC Longleaf HPC mode
# Matched to student's master branch for exact reproduction.
#
# Controller assignments:
#   cv       - BO targets: mclapply with ncores_grid=50
#   full     - Full-model seed fits: ninit_full=100 CPUs, 32 GB
#   med_mem  - Medium memory targets
#   low_mem  - Light tasks
#   default  - Sequential fallback

library(targets)
library(tarchetypes)
LOCAL_RENDER <- identical(Sys.getenv("DESURV_LOCAL_RENDER"), "1")
library(crew)
if (!LOCAL_RENDER) {
  library(crew.cluster)
}
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION <- "HEAD"

GIT_BRANCH <- tryCatch(
  gert::git_branch(),
  error = function(e) {
    branch <- tryCatch(
      trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE, ignore.stderr = TRUE)),
      error = function(e2) "unknown"
    )
    if (length(branch) == 0 || !nzchar(branch)) branch <- "unknown"
    branch
  }
)

DEFAULT_NINIT <- 50
DEFAULT_NINIT_FULL <- 100

# ------ Slurm controllers (Longleaf HPC) ------
SLURM_SCRIPT_LINES <- c("module load r/4.4.0")

default_controller <- crew_controller_sequential()

if (LOCAL_RENDER) {
  active_controller <- default_controller
} else {
  low_mem_controller <- crew_controller_slurm(
    name = "low_mem",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_per_cpu = 1,
      time_minutes = 120,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  cv_comp_controller <- crew_controller_slurm(
    name = "cv",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_per_cpu = 2,
      cpus_per_task = DEFAULT_NINIT,
      time_minutes = 720,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  full_run_controller <- crew_controller_slurm(
    name = "full",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = 32,
      cpus_per_task = DEFAULT_NINIT_FULL,
      time_minutes = 600,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  med_mem_controller <- crew_controller_slurm(
    name = "med_mem",
    workers = 120,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = 8,
      cpus_per_task = 1,
      time_minutes = 200,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  active_controller <- crew_controller_group(
    default_controller,
    low_mem_controller,
    cv_comp_controller,
    full_run_controller,
    med_mem_controller
  )
}

# ---- Global options ----
TARGET_PACKAGES <- c(
  "clusterProfiler","org.Hs.eg.db","DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr",
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2"
)

tar_option_set(
  packages = TARGET_PACKAGES,
  format = "rds",
  controller = active_controller,
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# Load subtype data if available
CMB_SUBTYPES_PATH <- "data/derv/cmbSubtypes_formatted.RData"
if (file.exists(CMB_SUBTYPES_PATH)) {
  load(CMB_SUBTYPES_PATH)
}
```

**Step 2: Commit**

```bash
git add local_slurm/longleaf/targets_setup.R
git commit -m "match student controllers: add low_mem/full, remove nmf, sequential default"
```

---

### Task 4: Update `targets_run_configs.R` — no changes needed (verify)

**Files:**
- Verify: `local_slurm/longleaf/targets_run_configs.R`

The student's values: `ninit_full = 100`, `std_nmf_k_grid = 2:12`, `coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10)`, `coxnet_alpha_grid = seq(0, 1, by = 0.1)`. Current longleaf config matches exactly. **No changes needed.**

---

### Task 5: Verify `targets_figure_configs.R` — no changes needed

**Files:**
- Verify: `local_slurm/longleaf/targets_figure_configs.R`

Figure configs are identical between student and longleaf. **No changes needed.**

---

### Task 6: Handle `nmf` controller reference in `targets_common_pipeline.R`

**Files:**
- Verify: `targets_common_pipeline.R:786`

**Critical issue:** `targets_common_pipeline.R` assigns `fit_std` to controller `"nmf"`, but the student's `targets_setup.R` has no `nmf` controller. Two possibilities:

1. The student's `master` does not use `targets_common_pipeline.R` at all — their `_targets.R` defines targets inline (confirmed: `master:_targets.R` has inline target definitions, not sourcing `targets_common_pipeline.R`)
2. Our branch DOES use `targets_common_pipeline.R`

**Resolution:** Since our pipeline goes through `targets_common_pipeline.R` but the student's doesn't, the `"nmf"` controller reference comes from our refactored code. The student's `fit_std` target has no explicit controller assignment (it falls to the default group). We need to either:
- (a) Change the `fit_std` controller in `targets_common_pipeline.R` from `"nmf"` to nothing/default, OR
- (b) Keep the `nmf` controller in `targets_setup.R` so that `targets_common_pipeline.R` doesn't error

**Recommended: option (b)** — keep the `nmf` controller definition but the student's `fit_std` would run on whatever controller the crew group falls back to. Actually, since crew dispatches to the FIRST controller with available capacity when an explicitly-named controller is not in the group, removing "nmf" from the group will cause a runtime error.

**Actual fix: change `targets_common_pipeline.R` line ~786** to remove the explicit `"nmf"` controller assignment so `fit_std` dispatches to the default (sequential), matching the student's behavior. OR keep `nmf` in the controller group since the student didn't have this refactored pipeline file.

**Decision needed from PI:** This is a pipeline code change, not just a config change. See Task 7 for the two options.

---

### Task 7: Decision — `fit_std` controller assignment

**Option A (minimal config-only change):** Keep the `nmf` controller in `targets_setup.R` despite the student not having it. The NMF target runs on a dedicated 10-CPU Slurm allocation. This is an infrastructure improvement that doesn't change results — NMF output is deterministic regardless of which controller runs it.

**Option B (exact reproduction):** Remove `nmf` from `targets_setup.R` AND change `targets_common_pipeline.R:~786` to remove the `controller = "nmf"` resource assignment. The `fit_std` target would then run on the sequential default controller, exactly matching the student.

**Recommendation:** Option A — the `nmf` controller is purely an execution resource. It doesn't affect the target hash or results. The student's NMF runs on whatever Slurm allocation the default group provides; ours runs on a dedicated one. Same output either way.

---

### Task 8: Verify resolved config hashes match

**Step 1: Print resolved config from student defaults**

After making changes in Tasks 1-3, run this verification:

```r
source("targets_bo_configs.R")
source("targets_val_configs.R")
source("targets_run_configs.R")
source("R/targets_config.R")

bo <- resolve_desurv_bo_configs(targets_bo_configs())
cat("BO config hash:", bo$tcgacptac$config_id, "\n")

# Print key resolved values to verify
cfg <- bo$tcgacptac
cat("k_grid upper:", cfg$desurv_bo_bounds$k_grid$upper, "\n")
cat("ngene_config:", cfg$ngene_config, "\n")
cat("tune_ngene:", cfg$tune_ngene, "\n")
cat("lambdah_config:", cfg$lambdah_config, "\n")
cat("tune_lambdah:", cfg$tune_lambdah, "\n")
cat("ntop_config:", cfg$ntop_config, "\n")
cat("bo_n_init:", cfg$bo_n_init, "\n")
cat("bo_max_refinements:", cfg$bo_max_refinements, "\n")
cat("bo_shrink_base:", cfg$bo_shrink_base, "\n")
cat("bo_importance_gain:", cfg$bo_importance_gain, "\n")
cat("coarse seed:", cfg$bo_coarse_control$seed, "\n")
cat("refine seed:", cfg$bo_refine_control$seed, "\n")
```

Expected output:
```
k_grid upper: 10
ngene_config: 500 5000
tune_ngene: TRUE
lambdah_config: 1e-07 100
tune_lambdah: TRUE
ntop_config: 50 250
bo_n_init: 20
bo_max_refinements: 2
bo_shrink_base: 0.5
bo_importance_gain: 0.3
coarse seed: 123
refine seed: 456
```

**Step 2: Commit all changes**

```bash
git add local_slurm/longleaf/
git commit -m "match longleaf configs to student master for exact reproduction"
```
