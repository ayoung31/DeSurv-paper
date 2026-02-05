# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Reproducible research repository for **DeSurv**: a survival-driven deconvolution tool for discovering prognostic cancer subtypes. Contains targets pipelines to reproduce all analyses, simulations, and manuscript figures from raw data.

**See [CHANGELOG.md](CHANGELOG.md) for recent changes and development history.**

**See [CONSISTENCY_STANDARDS.md](CONSISTENCY_STANDARDS.md) for code-documentation-paper alignment guidelines.**

## Essential Commands

```bash
# Install DeSurv package (required before running pipelines)
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'

# Run tests
Rscript -e 'testthat::test_dir("tests/testthat")'

# Pre-submission check (REQUIRED before sbatch)
./scripts/preflight_check.sh

# Submit pipelines to Slurm
sbatch _targets_wait.sh         # Main analysis (tcgacptac + bladder via config)
sbatch _targets_sims_local.sh   # Simulation studies

# Run specific pipeline directly
Rscript -e 'targets::tar_make(script = "_targets.R")'

# Rebuild simulation figures from pre-computed results (see below)
Rscript -e 'source("sim_figs.R"); sim_results_table <- readRDS("store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full/objects/sim_results_table"); save_sim_figs_by_scenario(build_sim_figs_by_scenario(sim_results_table), sim_dir="figures/sim", figure_configs=list())'

# Render manuscript
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'
```

**Note:** Bladder analysis runs via the same `_targets.R` pipeline with bladder configs in `targets_bo_configs.R`. There is no separate `_targets_bladder.R` file.

## Architecture

### Pipeline System

Uses `targets` R package for declarative workflow management with Slurm distribution via `crew.cluster`.

**Pipeline files:**
- `_targets.R` - Main TCGA/CPTAC pancreatic cancer analysis
- `_targets_bladder.R` - Bladder cancer analysis
- `_targets_sims.R` - Method validation simulations

**Configuration system:** All hyperparameters flow through separate config files (symlinks to `local_slurm/`):
- `targets_bo_configs.R` → `targets_bo_configs()` - Bayesian optimization bounds and settings
- `targets_run_configs.R` → `targets_run_configs()` - Full-model run parameters (references BO via `bo_key`)
- `targets_val_configs.R` → `targets_val_configs()` - External validation datasets and modes
- `targets_figure_configs.R` - Figure generation settings

Each config gets a hash-based `config_id` and readable `path_tag`. Validation runs at pipeline startup via `validate_desurv_configs()`.

**IMPORTANT:** There is NO `targets_configs.R` (monolithic file). The config files are split into separate modules. Old scripts referencing `targets_configs.R` are obsolete.

### Three-Phase Optimization Flow

1. **Bayesian Optimization** - Hyperparameter search (k, alpha, lambda, nu, ngene, ntop, penalties)
2. **Full-Model Run** - Best params with larger initialization count
3. **External Validation** - Evaluation on held-out datasets

### Key Modules

| File | Purpose |
|------|---------|
| `R/targets_config.R` | Config hashing, validation, diff utilities |
| `R/bo_helpers.R` | BO k-selection logic, bound management |
| `R/figure_targets.R` | All manuscript figures |
| `R/cluster_alignment.R` | Align discovered subtypes to references |
| `R/enrichment_map.R` | Pathway enrichment analysis |
| `targets_setup.R` | Library loading, Slurm controller definitions |
| `targets_common_pipeline.R` | Shared target definitions |

### Crew Resource Controllers (Local Desktop)

Defined in `targets_setup.R`. Tuned for 20 CPUs / 30GB RAM:

| Controller | Workers | Peak Memory | CPUs | Used By |
|-----------|---------|-------------|------|---------|
| `cv` | 1 | ~1.5 GB (+ 5 mclapply forks) | 6 | BO, NMF (`fit_std`), elbowk BO |
| `default` | 2 | ~1 GB | 2 | Figures, clustering, validation, ORA |
| `med_mem` | 2 | ~2 GB | 2 | Seed fits, consensus init |
| `full` | 2 | ~1 GB | 2 | Full model runs |
| `low_mem` | 2 | ~0.6 GB | spare | Light tasks |

**Memory budget:** ~12 GB peak (24 GB with 2× safety), leaving headroom on 30 GB system.

**CPU constraint:** Each cv worker forks `ncores_grid` sub-processes via `parallel::mclapply` inside DeSurv (one per CV fold). With `ncores_grid=5` and 2 cv workers, BO alone uses 12 CPUs. Do NOT increase cv workers without reducing `ncores_grid` in `targets_bo_configs.R`.

**IMPORTANT — avoid OOM crashes:**
- Changing `desurv_parallel_grid` or `desurv_ncores_grid` in `targets_bo_configs.R` **invalidates all BO targets** (they are inside `bo_config`, which is hashed). Only change these if you're willing to rerun BO from scratch.
- Worker counts and `seconds_idle` in `targets_setup.R` do NOT invalidate targets — safe to adjust anytime.
- All controllers use `seconds_idle = Inf` — workers stay alive for the entire pipeline run to prevent starvation during long BO tasks. Idle workers use no CPU and minimal RSS; `tar_make()` kills all on exit.
- If the pipeline crashes with crew worker crashes (6 consecutive), reduce non-cv workers first.

**HPC mode (Longleaf):** Uses `crew_controller_slurm` with 202 workers and higher memory. See Slurm section below.

### Data Modes

- `data_mode = "external"` - Separate train/validation files
- `data_mode = "split"` - Train/test split from single dataset
- `mode = "external"` / `mode = "train_split"` - Validation source

## DeSurv Package Dependency

The pipeline requires the DeSurv package to be **installed** (not just checked out). The pipeline uses `library(DeSurv)`, not `pkgload::load_all()`.

```bash
# Install before running pipelines
Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'
```

**Note:** If you modify the DeSurv package, you must reinstall it for changes to take effect in the pipeline.

**Core DeSurv functions used:**
- `desurv_data()` - Validate/encapsulate expression matrix + survival
- `desurv_cv_bayesopt_refine()` - BO with iterative bound refinement
- `desurv_fit()` - Main model fitting
- `predict.desurv_fit()` - Risk score predictions
- `desurv_get_top_genes()` - Factor interpretation

## Configuration Patterns

**Fixed vs tuned parameters:**
```r
# Single value = fixed
ngene_config = c(2000)

# Range = BO will tune
ngene_config = c(2000, 5000)
```

The helper `maybe_add_numeric_bound()` handles this logic.

**Compare configs:**
```r
desurv_config_diff(old_config, new_config, type = "bo")
```

## Paper & Manuscript

### File Structure

```
paper/
├── paper.Rmd                 # Main document (combines child sections)
├── 02_introduction_30102025.Rmd  # Introduction
├── 03_methods.Rmd            # Methods section
├── 04_results.Rmd            # Results section (loads most figures)
├── 05_discussion.Rmd         # Discussion section
├── supplement.Rmd            # Supplementary document
├── supp_methods.Rmd          # Supplementary methods (detailed algorithms)
├── references_30102025.bib   # Bibliography
├── pnas.csl                  # Citation style
├── pnas-new.cls              # PNAS LaTeX class
├── paper.pdf                 # Rendered output
└── supp_methods.pdf          # Rendered supplement
```

### Dependency Chain

```
Pipeline Targets (targets store)
       ↓
   tar_load() calls in .Rmd files
       ↓
paper/04_results.Rmd ──────────────────────────────┐
paper/supplement.Rmd                               │
       ↓                                           │
figures/panels/*.pdf (555 panel files)             │
figures/*.pdf (combined figures)                   │
       ↓                                           │
paper/paper.Rmd ← child includes ─────────────────┤
  ├── 02_introduction_30102025.Rmd                 │
  ├── 04_results.Rmd                               │
  ├── 05_discussion.Rmd                            │
  └── 03_methods.Rmd                               │
       ↓                                           │
paper/paper.pdf                                    │
paper/supp_methods.pdf                             │
```

### Key Targets Loaded in Results

The results section (`04_results.Rmd`) loads these targets. Most use `tar_load()`, but `sim_figs_by_scenario` uses `readRDS()` directly (see "Simulation Figures" below).

| Target | Figure | Content |
|--------|--------|---------|
| `sim_figs_by_scenario` | Fig 2E, 3 | Simulation results (loaded via `readRDS`, not `tar_load`) |
| `fig_bo_heat_tcgacptac` | Fig 2D | BO heatmap |
| `fig_residuals_tcgacptac` | Fig 2A | NMF residuals |
| `fig_cophenetic_tcgacptac` | Fig 2B | Cophenetic correlation |
| `fig_silhouette_tcgacptac` | Fig 2C | Silhouette width |
| `fig_gene_overlap_heatmap_*` | Fig 4A-B | Gene program correlations |
| `fig_variation_explained_*` | Fig 4C, 6A | Variance vs survival |
| `fig_hr_forest_tcgacptac` | Fig 5A | Forest plot |
| `fig_median_survival_*` | Fig 5B-C | KM curves |
| `tar_fit_desurv_tcgacptac` | Fig 6B | Fitted model for projection |

### Figure Generation

Figure targets are defined in `R/figure_targets.R` with this pattern:
- `build_fig_*_panels()` - Creates individual panel plots
- `combine_fig_*_panels()` - Arranges panels into composite figure
- `save_fig_*()` - Saves to `figures/` and `figures/panels/`

### Simulation Figures (HPC Store Workaround)

Simulation analyses (`sim_analysis_result`, 2400 branches) are run on HPC by the student. The aggregated `sim_results_table` and `sim_figs_by_scenario` objects are downloaded into the local store's `objects/` directory. However, the store **metadata** (`meta/`) is not transferred, so `targets` does not recognize these objects as up-to-date.

**DO NOT run `tar_make(names = "sim_figs_by_scenario", script = "_targets_sims.R")` locally** — it will attempt to recompute all 2,400 analysis branches from scratch and appear to freeze.

Instead, rebuild simulation figures directly from the pre-computed results table:

```r
source("sim_figs.R")
store <- "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"
sim_results_table <- readRDS(file.path(store, "objects", "sim_results_table"))
sim_figs_by_scenario <- build_sim_figs_by_scenario(sim_results_table)

# Save rebuilt object back to store (REQUIRED for paper rendering)
saveRDS(sim_figs_by_scenario, file.path(store, "objects", "sim_figs_by_scenario"))

# Also save standalone PDFs
save_sim_figs_by_scenario(sim_figs_by_scenario,
                          sim_dir = "figures/sim",
                          figure_configs = list())
```

**IMPORTANT:** You must save the rebuilt object back to the store. The paper (`04_results.Rmd:25`) loads `sim_figs_by_scenario` via `readRDS()` from the store. The HPC-built ggplot objects are **not renderable locally** due to ggplot2 version mismatch (`ggplot_build` fails on deserialized objects from different versions). Always rebuild from `sim_results_table` (raw data, version-agnostic) rather than reusing serialized ggplot objects.

### Rendering Commands

```bash
# Render main paper (requires targets store)
Rscript -e 'rmarkdown::render("paper/paper.Rmd")'

# Render supplementary methods
Rscript -e 'rmarkdown::render("paper/supp_methods.Rmd")'

# The paper uses this targets store (set in paper/_targets.yaml):
# store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full
```

### Targets Store for Paper

The paper reads from a specific targets store. **All four store references must match:**

| File | Setting |
|------|---------|
| `_targets.yaml` (root) | `main: store:` |
| `paper/_targets.yaml` | `main: store:` |
| `paper/paper.Rmd` | `params: tar_store:` |
| `_targets_sims_local.sh` | `tar_config_set(store = ...)` |

**Current store:** `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full`

To render with a different store:
```r
# In paper.Rmd, the store is set via params$tar_store
rmarkdown::render("paper/paper.Rmd",
  params = list(tar_store = "your_store_name"))
```

**Run `scripts/verify_consistency.R` to check all store paths match.**

## Code Style

- 2-space indentation, snake_case for functions/files
- Configuration edits go in `targets_*_configs.R` files, not hardcoded
- Keep helpers small and composable
- Commit messages: short, imperative, lowercase

### Shell Script Best Practices

**Directory handling (NEVER hardcode paths):**
```bash
# WRONG - breaks on other machines
cd /home/naimrashid/Downloads/DeSurv-paper

# RIGHT - portable
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"
```

**External command guards (prevent `set -e` crashes):**
```bash
# WRONG - crashes if sacct unavailable
failed=$(sacct -u $USER ...)

# RIGHT - graceful fallback
if command -v sacct &> /dev/null; then
    failed=$(sacct -u $USER ... || true)
else
    echo "sacct not available"
fi
```

**Multi-line command output (prevent arithmetic errors):**
```bash
# WRONG - sinfo may return multiple lines
CPUS=$(sinfo -h -o "%C")

# RIGHT - take first line
CPUS=$(sinfo -h -o "%C" | head -1)
```

**Process elapsed time (use etime, not CPU time):**
```bash
# WRONG - $10 in ps aux is CPU time
ps aux | awk '{print $10}'

# RIGHT - etime shows wall-clock elapsed
ps -eo pid,etime,args
```

**Fallible commands with `set -e` (grep exits 1 on no match):**
```bash
# WRONG - script exits if no matches
set -e
matches=$(grep "pattern" file)

# RIGHT - suppress exit code with || true
matches=$(grep "pattern" file || true)

# Also RIGHT - use if statement
if grep -q "pattern" file; then
    # handle match
fi
```

**Environment differences (local vs HPC):**
```bash
# WRONG - assumes module system exists
module load r/4.4.0

# RIGHT - graceful fallback
module load r/4.4.0 2>/dev/null || true

# For critical dependencies, check explicitly
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found in PATH"
    exit 1
fi
```

**Parsing multi-value output (avoid word splitting issues):**
```bash
# WRONG - breaks on spaces in values
for val in $(sinfo -h -o "%P %a"); do ...

# RIGHT - use read with IFS
sinfo -h -o "%P|%a" | while IFS='|' read partition avail; do
    echo "Partition: $partition, Available: $avail"
done
```

**`local` keyword only works inside functions:**
```bash
# WRONG - 'local' at top level or in while loop causes error with set -e
while read line; do
    local var=0    # ERROR: local only valid in function
done

# RIGHT - use plain variables outside functions
while read line; do
    var=0          # OK at any scope
done
```

**Capturing exit status with `|| true` (masks failures):**
```bash
# WRONG - || true makes $? always 0
output=$(some_command || true)
if [ $? -ne 0 ]; then  # Never triggers!
    echo "Failed"
fi

# RIGHT - capture status before || true
output=$(some_command 2>&1)
status=$?
if [ $status -ne 0 ]; then
    echo "Failed with exit $status"
fi
# Now safe to continue even if failed
```

## Consistency Checking

Before committing changes, run:
```bash
Rscript scripts/verify_consistency.R
```

This checks:
- Store path consistency across config files
- Config file validation
- Pipeline validation
- Debug statement detection
- Figure target status
- Documentation completeness

**Key documents for consistency:**
- [CONSISTENCY_STANDARDS.md](CONSISTENCY_STANDARDS.md) - Full guidelines
- [DEFAULTS.md](DEFAULTS.md) - Parameter defaults registry
- [CHANGELOG.md](CHANGELOG.md) - Change history

## Optional Dependencies

ORA/KEGG enrichment requires `clusterProfiler` and `org.Hs.eg.db`.

## Local Development & Config Architecture

### Config Symlink System

The root config files are **symlinks** to files in `local_slurm/`:
```
targets_bo_configs.R     → local_slurm/targets_bo_configs.R
targets_run_configs.R    → local_slurm/targets_run_configs.R
targets_val_configs.R    → local_slurm/targets_val_configs.R
targets_figure_configs.R → local_slurm/targets_figure_configs.R
```

The `local_slurm/` directory contains:
- Root level: Current active configs (symlink targets)
- `quick/`: Reduced iteration configs for testing
- `full/`: Publication-quality configs

**WARNING:** The `local_slurm/full/` configs may diverge from root configs. Always verify which config is active before running pipelines.

### Obsolete Files (Do Not Use)

| File | Status |
|------|--------|
| `local_slurm/targets_configs.R` | Obsolete monolithic file, never sourced |
| `local_slurm/submit_targets_quick.R` | References old file structure |
| `local_slurm/_targets_sims.R` | Incomplete stub, needs main file |

### Switching Configs

```bash
# Check current mode
cat local_slurm/.current_mode

# To use quick mode (for testing only):
# 1. Backup current symlinks
# 2. Point symlinks to local_slurm/quick/
# 3. Run pipeline
# 4. Restore symlinks
```

## Current Branch Configuration

**Branch:** `naimedits0125`
**Store:** `store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full`

### Local Desktop Mode

The pipeline is configured for local desktop execution (bypassing HPC Slurm issues):
- Main pipeline uses `crew_controller_local` with 2 workers
- Simulation pipeline uses `crew_controller_local` with 2 workers
- Resource limits: 19 CPUs, ~30GB RAM available

### Key Configuration Files

| File | Purpose |
|------|---------|
| `targets_bo_configs.R` | BO hyperparameter bounds (k, alpha, lambda, nu, ngene, ntop) |
| `targets_setup.R` | Controller definitions, package loading |
| `paper/_targets.yaml` | Store path for paper rendering |
| `paper/paper.Rmd` | Paper params including `tar_store` |

### Job Monitoring Tools

**IMPORTANT:** Always use these tools to prevent orphaned processes, stale jobs, and silent failures.

#### Pre-Submission Check (REQUIRED before sbatch)

```bash
# Run before submitting any job
./scripts/preflight_check.sh

# Or use as wrapper (runs command if checks pass)
./scripts/preflight_check.sh sbatch _targets_wait.sh
```

Checks for:
- Orphaned R/crew processes in project directory
- Stale/failed Slurm jobs from this project
- Crew backend lock files and socket states
- Targets store locks and errored targets
- System resources (CPU, memory, disk)
- Slurm cluster availability

#### Post-Submission Watchdog (for long-running jobs)

```bash
# One-time status check
./scripts/watchdog.sh

# Continuous monitoring (adaptive 5-15 min interval)
./scripts/watchdog.sh --watch

# With email alerts on issues
./scripts/watchdog.sh --watch --alert
```

Monitors for:
- Jobs running but not producing output (hung)
- Orphaned workers without controller (controller died)
- Log files not updating (30+ min stale)
- Store not receiving new results (60+ min stale)
- Jobs running > 24 hours
- NEW job failures (compared to previous check)
- System resource exhaustion

Logs to: `logs/watchdog.log`

#### Quick Status Commands

```bash
squeue -u $USER                    # View job queue
tail -f logs/slurm-<JOBID>.out     # Monitor job output
sacct -j <JOBID>                   # Job accounting info
```

#### Pipeline Runtime Diagnostics

**IMPORTANT:** When asked "what's running?" or "when will it finish?", **always run `tar_progress()` first**. Do NOT rely on the pipeline log, process trees, or store timestamps — these are unreliable for in-flight status.

**1. What is actually running right now?**
```r
# THE authoritative source — shows dispatched/completed/errored targets
p <- targets::tar_progress()
print(p[p$progress != "skipped", ])
```

**2. Estimate remaining runtime:**
```r
# Get historical runtimes for all outdated targets
m <- targets::tar_meta()
od <- targets::tar_outdated()
od_meta <- m[m$name %in% od, c("name", "seconds")]
od_meta <- od_meta[order(-od_meta$seconds), ]
print(od_meta)
cat("Total sequential runtime:", sum(od_meta$seconds, na.rm=TRUE)/3600, "hours\n")
```

**3. Identify the critical path (longest sequential chain):**

The cv worker is almost always the bottleneck. These targets run **sequentially** on cv:

```
fit_std → std_nmf_selected_k → desurv_bo_results_elbowk (~60 min)
                              → fit_std_desurvk (~4s)
                              → fit_std_elbowk (~3s)
```

Everything else (seed fits on med_mem, clustering/figures on default) runs in parallel and is rarely the bottleneck.

**4. Controller-to-target mapping (which worker runs what):**

| Controller | Key Targets | Typical Runtime |
|-----------|-------------|-----------------|
| `cv` | `desurv_bo_results_*`, `fit_std_*` | 1-60 min each |
| `med_mem` | `desurv_seed_fits_*`, `desurv_consensus_init_*`, `tar_fit_desurv_*` | 1-12 min each |
| `default` | clusters, validation, figures, ORA | 1-2 min each |

**5. Process-level confirmation (only if tar_progress is ambiguous):**
```bash
# Which controllers have active workers?
ps -eo pid,ppid,etime,rss,args | grep "crew_worker" | grep -v grep
# Look for controller="cv" vs controller="med_mem" in the command line

# Check mclapply children of a specific worker
ps --ppid <WORKER_PID> -o pid,etime,%cpu
```

**Known observability gaps in targets + crew:**
- The pipeline log only writes lines on dispatch/complete/error — long-running targets produce zero output for hours
- Crew worker stdout goes to `/dev/null` — `verbose=TRUE` output from BO/NMF is invisible
- The log "dispatched" event fires when a target is sent to the controller queue, NOT when a worker starts executing it — a dispatched target may be blocked on unsatisfied dependencies
- `tar_meta()` only has data for previously completed targets — there is no "elapsed time so far" for in-flight targets
- mclapply forks from different target types (BO CV folds vs NMF parallel runs) look identical in `ps` output

### Common Pitfalls to Avoid

| Pitfall | Symptom | Prevention |
|---------|---------|------------|
| **Orphaned crew workers** | Workers running after controller dies | Run `./scripts/preflight_check.sh` before new submissions |
| **Stale store locks** | "Store is locked" errors | Check for `.lock` files in store directory |
| **SIGTERM without cleanup** | Job killed, workers keep running | Use `./scripts/watchdog.sh --watch` to detect early |
| **Hardcoded paths** | Scripts fail on other machines | Use `${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}` |
| **sacct unavailable** | Scripts crash on `set -e` | Always guard with `command -v sacct` |
| **Multi-line sinfo output** | Arithmetic errors in bash | Use `head -1` when parsing sinfo |
| **CPU time vs elapsed time** | Miss long-running orphans | Use `ps -eo pid,etime,args` not `ps aux` |
| **grep with set -e** | Script exits on no match (exit 1) | Use `grep ... \|\| true` or `if grep -q` |
| **Local vs HPC environment** | Works locally, fails on cluster | Guard `module load` with `2>/dev/null \|\| true` |
| **Word splitting in loops** | Breaks on spaces in values | Use `IFS='|'` delimiter with `read` |
| **`local` outside function** | `set -e` aborts script | Only use `local` inside functions |
| **`|| true` masks exit status** | `$?` always 0, error check never triggers | Capture `$?` before `\|\| true` |
| **Cross-project resource contention** | Controller OOM crash, orphaned workers | Run `preflight_check.sh` - check section 2b warnings |
| **No monitoring during long runs** | Silent failures, hours of lost work | Always run `watchdog.sh --watch` or use `start_pipeline.sh` |
| **Skipping preflight** | Resource conflicts, stale locks | Use `./scripts/start_pipeline.sh` wrapper (enforces preflight) |
| **Piping pipeline to `head`/`tail`** | SIGPIPE kills controller, orphans workers | NEVER pipe `start_pipeline.sh` through `head`. Run in background with `nohup` |
| **Too many crew workers** | OOM crash kills pipeline mid-run, losing in-flight BO results | Keep total workers under memory budget (see `targets_setup.R` comments). Peak = sum of all controller workers × per-worker memory |
| **Changing `desurv_parallel_grid`/`ncores_grid`** | Invalidates ALL BO targets (multi-hour reruns) | These are inside `bo_config` which is hashed. Change worker counts in `targets_setup.R` instead (safe, no invalidation) |
| **Using pipeline log to diagnose status** | Log only updates on dispatch/complete events — silent for hours during long targets | Use `tar_progress()` as first diagnostic, not the log file |
| **Trusting "dispatched" in log** | "Dispatched" means queued, not executing — target may be blocked on dependencies | Cross-reference with `tar_progress()` and `tar_network()` dependency graph |
| **Using `ps` to identify running target** | mclapply forks from BO and NMF look identical (5 R processes at 100% CPU) | Use `tar_progress()` to identify the target; only use `ps` to confirm worker health |
| **Stale dispatched progress on restart** | Crashed pipeline leaves targets in "dispatched" state; `tar_make()` skips them | Use `start_pipeline.sh` (auto-clears) or `Rscript -e "targets::tar_destroy(destroy='progress')"` |
| **Paper target kills pipeline early** | Paper rendering error before crew initializes; `error="continue"` doesn't help | Use `start_pipeline.sh` (excludes paper by default); add `--with-paper` only when ready |
| **Worker idle timeout starvation** | med_mem/default workers die during long BO; downstream targets never run | All controllers now use `seconds_idle = Inf` — workers stay alive for entire run |
| **Running `tar_make` for sim figures locally** | Appears frozen — tries to recompute 2,400 analysis branches because store metadata is missing (objects were downloaded from HPC without meta) | Use `readRDS()` to load `sim_results_table` directly and call `build_sim_figs_by_scenario()`. See "Simulation Figures" section |

### Multi-Project Coordination

When running pipelines concurrently across projects (e.g., DeSurv + baton):

**Check total resource usage FIRST:**
```bash
# See all crew workers system-wide
ps aux | grep "crew::crew_worker" | grep -v grep

# Estimate total memory used
ps aux | grep "crew::crew_worker" | awk '{sum+=$6} END {print sum/1024/1024 "GB"}'

# Check which projects they belong to
for pid in $(pgrep -f "crew::crew_worker"); do
    echo "PID $pid: $(readlink -f /proc/$pid/cwd)"
done
```

**Resource guidelines for 30GB system:**
- Max ~22 crew workers total across all projects
- Each DeSurv BO task needs ~4GB (2 workers × nested parallelism)
- Leave 8GB headroom for system + non-crew tasks

**RECOMMENDED: Use the safe launcher:**
```bash
# Run in background - NEVER pipe through head/tail/less (SIGPIPE kills controller)
nohup ./scripts/start_pipeline.sh > logs/pipeline.log 2>&1 &
tail -f logs/pipeline.log             # Monitor separately

# Include paper rendering (excluded by default to prevent early errors killing the pipeline):
nohup ./scripts/start_pipeline.sh --with-paper > logs/pipeline.log 2>&1 &

# Or for simulations:
nohup ./scripts/start_pipeline.sh --sims > logs/sims.log 2>&1 &
```

This wrapper:
1. Runs mandatory preflight checks
2. Checks cross-project resource contention
3. Auto-clears stale "dispatched"/"started" progress from crashed runs
4. Excludes paper target by default (use `--with-paper` to include)
5. Auto-starts watchdog monitoring
6. Cleans up on completion

### Slurm Job Lifecycle

```
1. PRE-SUBMISSION
   └── Run ./scripts/preflight_check.sh
       ├── Fix any errors before proceeding
       └── Review warnings

2. SUBMISSION
   └── sbatch _targets_wait.sh (or wrapper mode)
       └── Job waits for system load < threshold

3. RUNNING
   └── Controller dispatches to crew workers
       ├── Monitor with: tail -f logs/slurm-<JOB>.out
       └── Or run: ./scripts/watchdog.sh --watch

4. IF PROBLEMS
   ├── Controller dies → Workers become orphaned
   │   └── Detection: watchdog shows "ORPHANED WORKERS"
   │   └── Fix: scancel <worker_job_ids>
   │
   ├── Workers die → Controller hangs
   │   └── Detection: Log stale for 30+ min
   │   └── Fix: scancel <controller_job_id>
   │
   └── Job killed externally → Partial state
       └── Detection: sacct shows FAILED/CANCELLED
       └── Fix: Run preflight, resubmit

5. COMPLETION
   └── Check store for results
       └── tar_meta() shows completed targets
```

### Recovery from Failed Jobs

```bash
# 1. Check what happened
sacct -j <JOBID> --format=JobID,JobName,State,ExitCode,MaxRSS,Elapsed

# 2. Clean up orphaned processes
./scripts/preflight_check.sh
# If orphans found:
kill <PIDs shown>
scancel <job_ids shown>

# 3. Clear stale progress (targets stuck in dispatched/started state)
Rscript -e 'targets::tar_destroy(destroy = "progress")'
# NOTE: start_pipeline.sh does this automatically

# 4. Check store state
Rscript -e 'targets::tar_meta()' | grep error

# 5. Remove stale locks if needed
rm store_*/.lock

# 6. Resubmit (start_pipeline.sh handles steps 3-5 automatically)
./scripts/preflight_check.sh sbatch _targets_wait.sh
```

## Known Issues & Code Review

See [CODE_REVIEW.md](CODE_REVIEW.md) for comprehensive code review findings including:
- **Critical bugs** - mostly resolved (see status below)
- **Architectural issues** in pipeline configuration and resource management
- **DeSurv package issues** in C++ numerical stability and validation
- **Test coverage gaps** for critical functions

### Critical Items (Status)

1. ~~**`R/compute_metrics.R:9`** - `sdZ` incorrectly assigned from `meanZ`~~ **RESOLVED**: Bug fixed and file removed (was dead code - see BUG1_IMPACT_ANALYSIS.md)
2. ~~**`R/select_best_init.R:5`** - Filtered dataframe `keep` never used~~ **RESOLVED**: File removed (was dead code - never called in pipeline)
3. ~~**`targets_common_pipeline.R`** - Remove `browser()` debug statements~~ **RESOLVED**: All 6 browser() statements removed
4. **`targets_setup.R:110`** - `DEFAULT_NINIT` referenced before definition

### DeSurv Package Critical Items

1. ~~**`src/functions.cpp:472`** - Beta backtracking logic inverted~~ **NOT A BUG**: Analysis confirmed the logic is correct for maximization (see BUG4_ANALYSIS.md)
2. **`R/cv_helpers.R:231`** - Missing gene subset validation in CV
3. **`R/predict_methods.R:71-106`** - Non-finite validation inconsistency

## Narrative Arc Analysis

See [NARRATIVE_ARC.md](NARRATIVE_ARC.md) for detailed analysis of the paper's storytelling structure:
- Five-act story structure (Problem → Solution → Proof → Insight → Generalization)
- Figure-by-figure narrative role analysis
- The central insight: "Variance ≠ Prognosis"
- Rhetorical strategies (problem-solution framing, escalating claims, consistent baseline)
- Information disclosure hierarchy (main text vs supplement)

## PNAS Submission Review

See [PNAS_REVIEW.md](PNAS_REVIEW.md) for comprehensive publication readiness assessment including:
- **PNAS compliance checklist** (word count, figures, references, required sections)
- **Critical code bugs** that may affect reported results
- **Paper-code inconsistencies** in mathematical formulations
- **Mathematical notation issues** (symbol overloading, typos, incomplete equations)
- **Scientific content assessment** (strengths, weaknesses, PNAS suitability)
- **Supplement completeness** review
- **Action checklist** organized by priority

### Key PNAS Issues

1. **Figures:** 6 figures (PNAS max: 4) - move 2 to supplement
2. **References:** 71 references (PNAS max: 50) - reduce by 21+
3. **Placeholders:** Keywords, author contributions, conflict of interest incomplete
4. **Methods section:** Too brief (~30 lines) for PNAS
5. **Discussion section:** Too brief (6 lines) for PNAS

## Figure Analysis

See [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md) for panel-by-panel evaluation of all manuscript figures:
- **Visualization effectiveness** - whether each panel optimally supports its narrative point
- **Proposed improvements** - alternative plotting approaches and annotations
- **Aesthetic assessment** - structure, proportions, color consistency, labels
- **Narrative arc impact** - how figure improvements strengthen the paper's story

### High-Impact Figure Changes (Minimal Effort)

1. **Figure 2:** Reorder panels B-C-D above A (problem before solution)
2. **Figure 4:** ADD FACTOR LABELS to heatmaps (dramatically improves interpretability)
3. **Figures 5-6:** Add HR, CI, p-value to all KM plots (standard practice)
4. **Figure 3:** Add "true k=3, n=100 replicates" annotation
5. **All figures:** Standardize legend terminology ("Standard NMF" vs "DeSurv")

## Suggested Replacement Text

See [SUGGESTED_TEXT.md](SUGGESTED_TEXT.md) for drop-in replacement language implementing editorial recommendations:
- **Revised Introduction** - 5-paragraph replacement with proper prior work acknowledgment
- **Related Work paragraph** - Explicitly acknowledges DECODER and prior PDAC work
- **6 targeted claim fixes** - Specific sentence replacements for W-vs-H, cross-cancer, overfitting claims
- **Revised Significance Statement** - More concrete alternative

### Key Text Changes

1. **Acknowledge Huang et al. and Le Goff et al.** - prior survival-NMF work cited fairly
2. **Acknowledge DECODER** - UNC prior work explicitly recognized
3. **Temper cross-cancer claim** - "consistent with prior reports" not novel discovery
4. **Add W-vs-H gradient explanation** - preempts mathematical reviewers
5. **Add explicit out-of-sample statement** - preempts double-dipping concerns

## Manuscript Pipeline Analysis

See [MANUSCRIPT_PIPELINE_ANALYSIS.md](MANUSCRIPT_PIPELINE_ANALYSIS.md) for comprehensive documentation of:
- **Complete dependency chain** from paper.pdf back to raw data
- **Figure-to-target mapping** showing which targets produce which figures
- **Configuration system** explaining BO, Run, and Validation configs
- **Simulation pipeline** with scenarios, analysis methods, and metrics
- **Real data pipeline** with datasets, preprocessing, and validation
- **Impact of configuration choices on paper claims**
- **Quick vs Full mode differences** and their effects on results

### Key Insights for Reproducibility

1. **Paper claims require both DeSurv AND NMF baselines** - `bo` and `bo_alpha0` analysis methods must both run
2. **Quick mode may produce different results** - fewer initializations, shorter BO, missing comparisons
3. **Figure 3 boxplots need alpha=0 runs** - only shows one box if NMF baseline is missing
4. **External validation requires all cohorts** - Dijk, Moffitt, PACA_AU_array, PACA_AU_seq, Puleo

### Critical Targets for Paper

| Target | Figure | Must Have |
|--------|--------|-----------|
| `sim_figs_by_scenario` | Fig 2E, 3 | Both bo and bo_alpha0 results |
| `fig_bo_heat_tcgacptac` | Fig 2D | Sufficient bo_n_iter |
| `tar_fit_desurv_tcgacptac` | Multiple | Stable factors (high ninit_full) |
| `fig_extval_*` | Fig 6 | All validation datasets |
