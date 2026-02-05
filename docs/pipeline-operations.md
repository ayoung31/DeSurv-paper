# Pipeline Operations Reference

Monitoring, diagnostics, pitfall avoidance, and recovery procedures for the DeSurv targets + crew pipelines.

## Job Monitoring Tools

**IMPORTANT:** Always use these tools to prevent orphaned processes, stale jobs, and silent failures.

### Pre-Submission Check (REQUIRED before sbatch)

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

### Post-Submission Watchdog (for long-running jobs)

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

### Quick Status Commands

```bash
squeue -u $USER                    # View job queue
tail -f logs/slurm-<JOBID>.out     # Monitor job output
sacct -j <JOBID>                   # Job accounting info
```

## Pipeline Runtime Diagnostics

**IMPORTANT:** When asked "what's running?" or "when will it finish?", **always run `tar_progress()` first**. Do NOT rely on the pipeline log, process trees, or store timestamps — these are unreliable for in-flight status.

### 1. What is actually running right now?

```r
# THE authoritative source — shows dispatched/completed/errored targets
p <- targets::tar_progress()
print(p[p$progress != "skipped", ])
```

### 2. Estimate remaining runtime

```r
# Get historical runtimes for all outdated targets
m <- targets::tar_meta()
od <- targets::tar_outdated()
od_meta <- m[m$name %in% od, c("name", "seconds")]
od_meta <- od_meta[order(-od_meta$seconds), ]
print(od_meta)
cat("Total sequential runtime:", sum(od_meta$seconds, na.rm=TRUE)/3600, "hours\n")
```

### 3. Identify the critical path (longest sequential chain)

The cv worker is almost always the bottleneck. These targets run **sequentially** on cv:

```
fit_std → std_nmf_selected_k → desurv_bo_results_elbowk (~60 min)
                              → fit_std_desurvk (~4s)
                              → fit_std_elbowk (~3s)
```

Everything else (seed fits on med_mem, clustering/figures on default) runs in parallel and is rarely the bottleneck.

### 4. Controller-to-target mapping (which worker runs what)

| Controller | Key Targets | Typical Runtime |
|-----------|-------------|-----------------|
| `cv` | `desurv_bo_results_*`, `fit_std_*` | 1-60 min each |
| `med_mem` | `desurv_seed_fits_*`, `desurv_consensus_init_*`, `tar_fit_desurv_*` | 1-12 min each |
| `default` | clusters, validation, figures, ORA | 1-2 min each |

### 5. Process-level confirmation (only if tar_progress is ambiguous)

```bash
# Which controllers have active workers?
ps -eo pid,ppid,etime,rss,args | grep "crew_worker" | grep -v grep
# Look for controller="cv" vs controller="med_mem" in the command line

# Check mclapply children of a specific worker
ps --ppid <WORKER_PID> -o pid,etime,%cpu
```

### Known observability gaps in targets + crew

- The pipeline log only writes lines on dispatch/complete/error — long-running targets produce zero output for hours
- Crew worker stdout goes to `/dev/null` — `verbose=TRUE` output from BO/NMF is invisible
- The log "dispatched" event fires when a target is sent to the controller queue, NOT when a worker starts executing it — a dispatched target may be blocked on unsatisfied dependencies
- `tar_meta()` only has data for previously completed targets — there is no "elapsed time so far" for in-flight targets
- mclapply forks from different target types (BO CV folds vs NMF parallel runs) look identical in `ps` output

## Common Pitfalls

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
| **Running `tar_make` for sim figures locally** | Appears frozen — tries to recompute 2,400 analysis branches because store metadata is missing (objects were downloaded from HPC without meta) | Use `readRDS()` to load `sim_results_table` directly and call `build_sim_figs_by_scenario()`. See "Simulation Figures" in CLAUDE.md |

## Multi-Project Coordination

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

## Slurm Job Lifecycle

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

## Recovery from Failed Jobs

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
