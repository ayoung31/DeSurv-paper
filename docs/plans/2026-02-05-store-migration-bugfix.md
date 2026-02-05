# Store Migration: Rerun Non-Simulation Targets with Updated DeSurv

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a new targets store that uses the updated DeSurv package (branch `20260107bugfix`) for all non-simulation targets, while preserving simulation results from the student's HPC runs. Switch from `crew_controller_local` to `crew_controller_slurm` with resource-aware scheduling so workers only start when CPUs/RAM are available.

**Architecture:** Reconfigure crew controllers for Slurm with profiled resource limits, copy simulation objects into a new store, update all store path references, rebuild sim figures from raw data, then run `tar_make()` to recompute everything non-simulation with the bugfix DeSurv.

**Tech Stack:** R targets, DeSurv (branch 20260107bugfix, v1.0.1), crew + crew.cluster (Slurm controllers)

---

## Context

- **Old store:** `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full` (935 objects, ~13 GB)
- **New store:** `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`
- **Sim objects to copy:** 410 objects (~3.5 GB), prefixed `sim_*`
- **Non-sim objects to recompute:** 525 objects (~9.5 GB) — BO, fits, validation, figures, paper
- **Disk space available:** ~392 GB (plenty)
- **DeSurv installed version:** 1.0.1, branch `20260107bugfix`, commit `370c88a`
- **System:** 20 CPUs, 31 GB RAM, local Slurm with cgroup enforcement
- **Note:** BATON workers may be running concurrently (4 workers × 4 CPUs = 16 CPUs, 8 GB). Slurm will queue DeSurv jobs until resources free up.

## Resource Profiling (from old store metadata)

| Controller | Target Pattern | Count | Runtime | CPUs Needed | Memory Needed |
|-----------|---------------|-------|---------|-------------|---------------|
| **cv** (BO) | `desurv_bo_results_*` | 6 | 2.7-4.2 hrs | 6 (1 parent + 5 mclapply forks from `ncores_grid=5`) | ~5 GB |
| **nmf** (new) | `fit_std_*` | 6 | 3 sec - 20 min | 10 (NMF `.options="p30"` forks 30, but capped to allocation) | ~4 GB |
| **med_mem** | `desurv_seed_fits_*` | 6 | 6-29 min | 1 (sequential for-loop, `parallel_init=FALSE`) | ~4 GB (350 MB output) |
| **med_mem** | `tar_fit_desurv_*` | 6 | 1-14 sec | 1 | ~2 GB |
| **default** | clusters, ORA, figs, validation, data | ~400+ | <2 min each | 1 | ~2 GB |

**Key insight:** BO targets and fit_std currently share the "cv" controller but need very different CPUs (6 vs up to 30). Splitting into separate controllers prevents waste. Changing `resources = tar_resources(crew = ...)` does NOT invalidate target hashes — it's deployment metadata only.

## Files That Need Store Path Updates

### Critical (pipeline execution & paper rendering):

| # | File | Line | Type |
|---|------|------|------|
| 1 | `_targets.yaml` | 2 | Root pipeline store config |
| 2 | `paper/_targets.yaml` | 2 | Paper rendering store config |
| 3 | `paper/paper.Rmd` | 51 | Paper params `tar_store` |
| 4 | `_targets_sims_local.sh` | 28 | Sim pipeline store config |
| 5 | `submit_targets_full.R` | 2, 5 | Full submission script |
| 6 | `submit_targets_sims_full.R` | 2, 5 | Full sim submission script |
| 7 | `_targets_full.sh` | 9 | Full pipeline shell echo |
| 8 | `_targets_sims_full.sh` | 9 | Full sim pipeline shell echo |

### Documentation (update for consistency):

| # | File | Lines | Type |
|---|------|-------|------|
| 9 | `CLAUDE.md` | 33, 186, 211, 225, 295 | Project instructions |
| 10 | `CONSISTENCY_STANDARDS.md` | 126, 141 | Standards doc |

### Skip (old/obsolete/dynamic — no update needed):

- `paper/old/*.Rmd` — archived, not active
- `submit_targets.R`, `submit_targets_sims.R` — dynamic (use `PKG_VERSION`/`GIT_BRANCH` vars)
- `CHANGELOG.md`, `BUG1_IMPACT_ANALYSIS.md`, `MANUSCRIPT_PIPELINE_ANALYSIS.md` — historical records
- `local_slurm/quick/*`, `local_slurm/submit_targets_quick.R` — separate quick configs
- `paper/algorithms/_targets.yaml` — separate sub-pipeline
- `inst/render.R` — old render script
- `scripts/watchdog.sh` — uses glob pattern, no hardcoded store name

---

### Task 1: Switch targets_setup.R to Slurm Controllers with Profiled Resources

**Files:**
- Modify: `targets_setup.R` (lines 62-125)
- Modify: `targets_common_pipeline.R` (line 786: change fit_std controller from "cv" to "nmf")

**Rationale:** The current `targets_setup.R` uses `crew_controller_local` which bypasses Slurm. Switching to `crew_controller_slurm` ensures workers only start when enough CPUs/RAM are available. The resource limits per controller are derived from actual runtime profiling of the old store.

**Step 1: Replace crew controllers in `targets_setup.R`**

Replace lines 62-125 (the local controller block) with:

```r
# ------ Slurm controllers (resource-profiled for 20 CPU / 31 GB desktop) ------
# Slurm + cgroups enforces CPU/memory limits per worker.
# Workers queue until resources are available, preventing OOM crashes.
#
# Resource budget (from tar_meta() profiling of previous full run):
#   cv (BO):      6 CPUs × 5 GB = fits 3 concurrent (18 CPUs, 15 GB)
#   nmf:         10 CPUs × 4 GB = fits 2 concurrent (20 CPUs, 8 GB)
#   med_mem:      1 CPU  × 4 GB = fits 5 concurrent (5 CPUs, 20 GB)
#   default:      1 CPU  × 2 GB = fits 10 concurrent (10 CPUs, 20 GB)
# Slurm schedules dynamically — these are per-worker requests, not fixed reservations.

default_controller <- crew_controller_slurm(
  name = "default",
  workers = 10,
  seconds_idle = 120,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_per_cpu = 2,
    time_minutes = 30,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# BO targets: desurv_cv_bayesopt uses mclapply(ncores_grid=5) internally
# Each worker needs 6 CPUs (1 parent + 5 forks). Longest BO run: 4.2 hrs.
cv_comp_controller <- crew_controller_slurm(
  name = "cv",
  workers = 3,
  seconds_idle = 300,
  options_cluster = crew_options_slurm(
    cpus_per_task = 6,
    memory_gigabytes_per_cpu = 1,
    time_minutes = 300,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# NMF k-selection: NMF::nmf(.options="p30") forks up to 30 processes.
# With 10 CPUs allocated, NMF parallelism is capped by cgroups.
# Only 2 of 6 fit_std targets are slow (20 min); rest are <6 sec.
nmf_controller <- crew_controller_slurm(
  name = "nmf",
  workers = 2,
  seconds_idle = 120,
  options_cluster = crew_options_slurm(
    cpus_per_task = 10,
    memory_gigabytes_per_cpu = 1,
    time_minutes = 60,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# Seed fits: sequential for-loop (ninit_full=100, parallel_init=FALSE).
# Single-threaded but memory-heavy (up to 350 MB output). Longest: 29 min.
med_mem_controller <- crew_controller_slurm(
  name = "med_mem",
  workers = 4,
  seconds_idle = 120,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_required = 4,
    time_minutes = 120,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

active_controller <- crew_controller_group(
  default_controller,
  cv_comp_controller,
  nmf_controller,
  med_mem_controller
)
```

**Step 2: Reassign fit_std from "cv" to "nmf" controller**

In `targets_common_pipeline.R`, change line 786:

```r
# Old:
      crew = tar_resources_crew(controller = "cv")
# New:
      crew = tar_resources_crew(controller = "nmf")
```

This is the `fit_std` target's resource annotation. Changing controller assignment does NOT invalidate the target hash (only `resources` metadata changes, not the `command`).

**Step 3: Ensure logs directory exists**

```bash
mkdir -p ~/Downloads/DeSurv-paper/logs
```

**Step 4: Verify Slurm is working**

```bash
sinfo
# Expected: debug partition UP with 1 node
```

---

### Task 2: Create New Store Directory and Copy Sim Objects

**Files:**
- Create: `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/` (directory)
- Read from: `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full/objects/sim_*`

**Step 1: Create the new store directory structure**

```bash
cd ~/Downloads/DeSurv-paper
mkdir -p "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/objects"
mkdir -p "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/meta"
mkdir -p "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/user"
mkdir -p "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/workspaces"
```

**Step 2: Copy all sim objects from old store to new store**

```bash
cd ~/Downloads/DeSurv-paper
cp store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full/objects/sim_* \
   "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/objects/"
```

**Step 3: Verify the copy**

```bash
# Count should be 410
ls "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/objects/" | wc -l
# Spot check a key object exists
ls -la "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full/objects/sim_results_table"
```

Expected: 410 objects, `sim_results_table` present (~2 MB).

---

### Task 3: Update Critical Store Path References

**Step 1: Update root `_targets.yaml`**

- Modify: `_targets.yaml:2`

```yaml
main:
  store: store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full
  script: _targets.R
```

**Step 2: Update `paper/_targets.yaml`**

- Modify: `paper/_targets.yaml:2`

```yaml
main:
  store: store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full
```

**Step 3: Update `paper/paper.Rmd` params**

- Modify: `paper/paper.Rmd:51`

Change:
```yaml
  tar_store: "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"
```
To:
```yaml
  tar_store: "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full"
```

**Step 4: Update `_targets_sims_local.sh`**

- Modify: `_targets_sims_local.sh:28`

Change:
```r
tar_config_set(store = "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full")
```
To:
```r
tar_config_set(store = "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full")
```

**Step 5: Update `submit_targets_full.R`**

- Modify: `submit_targets_full.R:2,5`

Change store name in comment (line 2) and `STORE_NAME` variable (line 5) to `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`.

**Step 6: Update `submit_targets_sims_full.R`**

- Modify: `submit_targets_sims_full.R:2,5`

Same change as Step 5.

**Step 7: Update `_targets_full.sh`**

- Modify: `_targets_full.sh:9`

Change echo line to reference new store name.

**Step 8: Update `_targets_sims_full.sh`**

- Modify: `_targets_sims_full.sh:9`

Change echo line to reference new store name.

**Step 9: Verify all critical references are updated**

```bash
cd ~/Downloads/DeSurv-paper
grep -rn "naimedits0125_full" \
  _targets.yaml paper/_targets.yaml paper/paper.Rmd \
  _targets_sims_local.sh submit_targets_full.R \
  submit_targets_sims_full.R _targets_full.sh _targets_sims_full.sh
```

Expected: **zero matches** (all updated to `20260107bugfix_full`).

---

### Task 4: Update Documentation Store References

**Step 1: Update `CLAUDE.md`**

- Modify: `CLAUDE.md` — 5 occurrences (lines 33, 186, 211, 225, 295)

Replace all instances of `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full` with `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`.

Also update line 294 branch reference from `naimedits0125` to `20260107bugfix`.

**Step 2: Update `CONSISTENCY_STANDARDS.md`**

- Modify: `CONSISTENCY_STANDARDS.md` — 2 occurrences (lines 126, 141)

Replace store name similarly.

**Step 3: Verify documentation updates**

```bash
cd ~/Downloads/DeSurv-paper
grep -rn "naimedits0125_full" CLAUDE.md CONSISTENCY_STANDARDS.md
```

Expected: **zero matches** in these files.

---

### Task 5: Rebuild Simulation Figures in New Store

Per CLAUDE.md instructions, HPC-built ggplot objects are not renderable locally due to version mismatch. Must rebuild from raw `sim_results_table`.

**Step 1: Rebuild `sim_figs_by_scenario` and save to new store**

```bash
cd ~/Downloads/DeSurv-paper
Rscript -e '
source("sim_figs.R")
store <- "store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full"
sim_results_table <- readRDS(file.path(store, "objects", "sim_results_table"))
sim_figs_by_scenario <- build_sim_figs_by_scenario(sim_results_table)

# Save rebuilt object back to store (REQUIRED for paper rendering)
saveRDS(sim_figs_by_scenario, file.path(store, "objects", "sim_figs_by_scenario"))

# Also save standalone PDFs
save_sim_figs_by_scenario(sim_figs_by_scenario,
                          sim_dir = "figures/sim",
                          figure_configs = list())
cat("Done. sim_figs_by_scenario rebuilt and saved.\n")
'
```

Expected: Script completes without error, `sim_figs_by_scenario` updated in new store, PDFs written to `figures/sim/`.

---

### Task 6: Validate Setup Before Running Pipeline

**Step 1: Confirm DeSurv package version is from bugfix branch**

```bash
Rscript -e '
cat("DeSurv version:", packageDescription("DeSurv", fields = "Version"), "\n")
cat("Built:", packageDescription("DeSurv", fields = "Built"), "\n")
'
```

Expected: Version `1.0.1`, built today (2026-02-05).

**Step 2: Confirm targets sees the new store and Slurm controllers**

```bash
cd ~/Downloads/DeSurv-paper
Rscript -e '
library(targets)
cat("Active store:", tar_config_get("store"), "\n")
cat("Store exists:", dir.exists(tar_config_get("store")), "\n")
cat("Sim objects:", length(list.files(file.path(tar_config_get("store"), "objects"), pattern = "^sim_")), "\n")
'
```

Expected: Store is `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full`, exists, 410 sim objects.

**Step 3: Dry-run to verify Slurm controller configuration loads**

```bash
cd ~/Downloads/DeSurv-paper
Rscript -e '
source("targets_setup.R")
cat("Controllers in group:", paste(names(active_controller$controllers), collapse=", "), "\n")
cat("Setup OK\n")
'
```

Expected: Controllers listed as `default, cv, nmf, med_mem`.

---

### Task 7: Run the Pipeline (user executes)

**Step 1: Run `tar_make()` for the main pipeline**

```bash
cd ~/Downloads/DeSurv-paper
Rscript -e '
library(targets)
tar_make(script = "_targets.R")
'
```

This will recompute all non-simulation targets (BO, model fits, validation, figures, paper) using the updated DeSurv package. Slurm will schedule workers based on available resources.

**Note:** This is a long-running step. Monitor with:
```bash
# Target-level progress
Rscript -e 'p <- targets::tar_progress(); print(p[p$progress != "skipped",])'

# Slurm job queue
squeue -u $USER -o "%.8i %.9P %.12j %.2t %.10M %.6C %.10m"
```

**Expected total runtime:** ~25 hrs (dominated by 6 BO targets at 3-4 hrs each, running up to 3 concurrently → ~8-12 hrs for BO phase, then seed_fits, then downstream).

**Step 2: Verify pipeline completed successfully**

```bash
cd ~/Downloads/DeSurv-paper
Rscript -e '
library(targets)
p <- tar_progress()
cat("Completed:", sum(p$progress == "completed"), "\n")
cat("Errored:", sum(p$progress == "errored"), "\n")
cat("Skipped:", sum(p$progress == "skipped"), "\n")
if (any(p$progress == "errored")) {
  cat("\nErrored targets:\n")
  print(p$name[p$progress == "errored"])
}
'
```

Expected: All targets completed, zero errored.

---

## Slurm Controller Summary

| Controller | Targets | Workers | CPUs/Worker | RAM/Worker | Time Limit | Rationale |
|-----------|---------|---------|-------------|------------|------------|-----------|
| **default** | clustering, ORA, figs, validation, data (~400) | 10 | 1 | 2 GB | 30 min | Single-threaded, light tasks. 10 workers keeps queue moving. |
| **cv** | `desurv_bo_results_*` (6) | 3 | 6 | 6 GB | 5 hrs | mclapply forks 5 children (ncores_grid=5). 3×6=18 CPUs max. |
| **nmf** | `fit_std_*` (6) | 2 | 10 | 10 GB | 1 hr | NMF forks up to 30 procs, capped by cgroup to 10 CPUs. |
| **med_mem** | `desurv_seed_fits_*`, `tar_fit_desurv_*` (12) | 4 | 1 | 4 GB | 2 hrs | Sequential fitting loop, single-threaded but memory-heavy. |

## Summary of Changes

| What | Old Value | New Value |
|------|-----------|-----------|
| Store name | `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full` | `store_PKG_VERSION=NA_GIT_BRANCH=20260107bugfix_full` |
| DeSurv branch | `naimedits0125` (or older) | `20260107bugfix` |
| DeSurv version | older | 1.0.1 (commit 370c88a) |
| Crew controllers | `crew_controller_local` (no resource limits) | `crew_controller_slurm` (CPU/RAM/time enforced) |
| Controller count | 5 (default, low_mem, cv, full, med_mem) | 4 (default, cv, nmf, med_mem) |
| fit_std controller | cv (shared with BO) | nmf (dedicated, 10 CPUs) |
| Files modified | 12 files (8 store refs + 2 docs + 2 controller config) | |
| Sim objects | Copied as-is (410 objects, ~3.5 GB) | |
| Non-sim objects | Recomputed by `tar_make()` | |

## Rollback

If issues arise, the old store is untouched:
```bash
# Revert all file changes
cd ~/Downloads/DeSurv-paper && git checkout -- .

# Old store is still at:
# store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full/
```
