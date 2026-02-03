# Incident Postmortem: Pipeline Controller Crash (2026-02-02)

## Summary

**What happened:** Two incidents on the same day:
1. **13:52** - DeSurv controller crashed silently, leaving orphaned workers for 3+ hours
2. **19:17** - After implementing fixes, a race condition in the new safe launcher allowed the pipeline to start despite resource warnings, crashing the entire system

**Impact:**
- ~3 hours of BO computation lost (first incident)
- Complete system unresponsiveness requiring forced restart (second incident)

**Root cause:**
1. First incident: Resource contention without coordination
2. Second incident: Race condition between preflight check and start_pipeline.sh re-check

## Timeline

### Incident 1: Silent Controller Crash

| Time | Event |
|------|-------|
| 09:48 | DeSurv pipeline started (`desurv_bo_results_tcgacptac` dispatched) |
| ~13:52 | Controller crashed (last meta update) |
| 13:52-16:54 | Workers continued running without controller (orphaned) |
| 16:54 | Crash detected during manual status check |
| 16:57 | Cleanup initiated |

### Incident 2: System Crash from Race Condition (after "fix")

| Time | Event |
|------|-------|
| 19:16:21 | Created `start_pipeline.sh` with preflight + watchdog |
| 19:16:21 | Preflight detected 21 baton workers using 13.9GB |
| 19:16:21 | Preflight WARNED: Combined 28.9GB > 23GB threshold |
| 19:16:2x | start_pipeline.sh Step 2 re-checked: saw only 15 workers (6.3GB) |
| 19:16:27 | Step 2 INCORRECTLY passed: 21.3GB < 23GB threshold |
| 19:16:27 | Pipeline started, spawned 22 new workers |
| ~19:17 | 43 total workers demanding ~45GB on 30GB+8GB system |
| ~19:17+ | System thrashing, load spiked to 39+ (20 CPUs) |
| ~19:22 | System unresponsive, forced restart |

**The Bug:** Preflight and start_pipeline.sh Step 2 both checked worker counts, but at different times. Workers were dynamically starting/stopping (baton pipeline active), so Step 2 got different numbers and overrode preflight's warning.

## Root Cause Analysis

### Primary Cause: Resource Contention

Two pipelines running concurrently:

| Pipeline | Workers | Memory per Worker | Total |
|----------|---------|-------------------|-------|
| DeSurv (cv controller) | 2 | ~4GB | ~8GB |
| DeSurv (nested parallel) | 2×5=10 | ~1GB | ~10GB |
| Baton (slurm_calibration) | ~20 | ~0.8GB | ~16GB |
| **Combined** | **~32** | - | **~34GB** |

System had 30GB RAM + 8GB swap = 38GB total. Combined workload exceeded safe thresholds.

### Contributing Factors

1. **Preflight check not run** - Would have detected cross-project workers
2. **Watchdog not running** - Would have detected crash within 30 min
3. **No automatic monitoring** - Required manual intervention
4. **Memory threshold too lenient** - Preflight warns at 20% free, but crash happened at ~23% used

## What Went Wrong

| Control | Status | Gap |
|---------|--------|-----|
| `preflight_check.sh` | EXISTS but not run | No enforcement |
| `watchdog.sh` | EXISTS but not run | Not automatic |
| Cross-project detection | EXISTS in preflight | Only warns, doesn't block |
| Memory monitoring | 20% threshold | Too lenient for concurrent workloads |
| Process lock | Exists | Only protects same-project, not cross-project |

## Recommended Fixes

### 1. MANDATORY Preflight Enforcement

**Location:** `_targets.R` startup section

Add to beginning of `_targets.R`:
```r
# Preflight safety check
if (interactive() && Sys.getenv("DESURV_SKIP_PREFLIGHT") != "1") {
  preflight_result <- system("./scripts/preflight_check.sh --quiet", intern = TRUE)
  if (attr(preflight_result, "status") != 0) {
    stop("Preflight check failed. Run ./scripts/preflight_check.sh for details.")
  }
}
```

### 2. Auto-Start Watchdog

**Location:** `scripts/start_pipeline.sh` (new wrapper script)

```bash
#!/bin/bash
# Start pipeline with automatic monitoring

# Run preflight (mandatory)
./scripts/preflight_check.sh || exit 1

# Start watchdog in background
nohup ./scripts/watchdog.sh --watch --alert > logs/watchdog.log 2>&1 &
WATCHDOG_PID=$!
echo "Watchdog started: PID $WATCHDOG_PID"

# Start pipeline
Rscript -e 'targets::tar_make()'
EXIT_CODE=$?

# Stop watchdog
kill $WATCHDOG_PID 2>/dev/null

exit $EXIT_CODE
```

### 3. Enhanced Resource Checking

**Location:** `scripts/preflight_check.sh` - Update Section 6

Add combined resource estimation:
```bash
# Calculate TOTAL resource usage across all crew workers
total_crew_mem=0
total_crew_count=0
for pid in $(pgrep -f "crew::crew_worker"); do
    mem=$(ps -o rss= -p $pid 2>/dev/null || echo 0)
    total_crew_mem=$((total_crew_mem + mem))
    total_crew_count=$((total_crew_count + 1))
done

total_crew_mem_gb=$(echo "scale=1; $total_crew_mem / 1048576" | bc)
echo "   Total crew workers system-wide: $total_crew_count using ${total_crew_mem_gb}GB"

# Estimate new pipeline resource needs
expected_workers=12  # 2 cv + 10 nested
expected_mem_gb=18   # Estimated

combined_usage_gb=$(echo "scale=1; $total_crew_mem_gb + $expected_mem_gb" | bc)
mem_total_gb=$(echo "scale=1; $mem_total / 1048576" | bc)

if (( $(echo "$combined_usage_gb > $mem_total_gb * 0.8" | bc -l) )); then
    check_fail "Combined workload would exceed 80% memory (${combined_usage_gb}GB / ${mem_total_gb}GB)"
fi
```

### 4. Stricter Memory Thresholds

**Location:** `scripts/preflight_check.sh` line 319

Change from:
```bash
if [ "$mem_pct" -lt 20 ]; then
```

To:
```bash
if [ "$mem_pct" -lt 40 ]; then
    check_fail "Insufficient available memory: ${mem_pct}% free (need 40%+)"
elif [ "$mem_pct" -lt 30 ]; then
    check_warn "Low available memory: ${mem_pct}% free"
```

### 5. Cross-Project Coordination Protocol

**Location:** `CLAUDE.md` - Add new section

```markdown
### Multi-Project Coordination

When running pipelines concurrently across projects:

1. **Check total resource usage FIRST:**
   ```bash
   # See all crew workers system-wide
   ps aux | grep "crew::crew_worker" | grep -v grep

   # Estimate total memory
   ps aux | grep "crew::crew_worker" | awk '{sum+=$6} END {print sum/1024/1024 "GB"}'
   ```

2. **Coordinate with other pipelines:**
   - If baton/adaptive-trial is running: Wait or reduce DeSurv workers
   - Total crew workers should not exceed: (Total RAM - 8GB) / 1GB
   - For 30GB system: max ~22 workers total across all projects

3. **Use start_pipeline.sh wrapper:**
   ```bash
   ./scripts/start_pipeline.sh  # Runs preflight + watchdog automatically
   ```
```

### 6. CLAUDE.md Documentation Updates

Add to "Common Pitfalls to Avoid" table:

| Pitfall | Symptom | Prevention |
|---------|---------|------------|
| **Cross-project resource contention** | Controller OOM crash, orphaned workers | Run `preflight_check.sh` - check section 2b warnings |
| **No monitoring during long runs** | Silent failures, hours of lost work | Always run `watchdog.sh --watch` in background |
| **Skipping preflight** | Resource conflicts, stale locks | Use `start_pipeline.sh` wrapper (enforces preflight) |

### 7. CRITICAL FIX: Race Condition in Resource Checks (Incident 2)

**Problem:** `start_pipeline.sh` Step 2 duplicated the resource check from preflight. Due to timing differences, it got different worker counts and OVERRODE preflight's warning.

**Fix 1: Make preflight resource warning a BLOCKING error**

**Location:** `scripts/preflight_check.sh` combined crew worker section

```bash
# Changed from check_warn to check_fail:
if (( $(echo "$combined_gb > $threshold_gb" | bc -l) )); then
    check_fail "CRITICAL: System-wide crew workers: $total_crew_count using ${total_crew_mem_gb}GB"
    echo "         This EXCEEDS 75% threshold - BLOCKING"
```

**Fix 2: Remove redundant check from start_pipeline.sh**

**Location:** `scripts/start_pipeline.sh` Step 2

```bash
# Removed the re-check entirely. Now just confirms preflight passed:
echo "   Preflight passed - proceeding with pipeline"
echo "   (Resource checks handled by preflight to avoid race conditions)"
```

**Principle:** Single source of truth for resource decisions. Never re-check dynamic state that could change between checks.

## Action Items

| Priority | Item | Owner | Status |
|----------|------|-------|--------|
| P0 | Kill orphaned DeSurv workers | - | DONE |
| P0 | Remove stale locks | - | DONE |
| P0 | Fix race condition in start_pipeline.sh | - | **DONE** |
| P0 | Make preflight resource check BLOCKING | - | **DONE** |
| P1 | Create `scripts/start_pipeline.sh` wrapper | - | **DONE** |
| P1 | Update preflight memory thresholds (25%/40%) | - | **DONE** |
| P1 | Add combined resource estimation to preflight | - | **DONE** |
| P2 | Update CLAUDE.md with multi-project section | - | **DONE** |
| P2 | Add preflight enforcement to `_targets.R` | - | TODO |
| P3 | Consider cgroup/resource limits | - | TODO |

## Lessons Learned

1. **Tools exist but process failed** - Preflight and watchdog were available but not used
2. **Warnings should be errors** - Cross-project workers were a warning, should block by default
3. **Automation beats discipline** - Relying on manual checks doesn't scale
4. **Monitor long-running jobs** - Any job >1 hour needs active monitoring
5. **CRITICAL: Never duplicate dynamic checks** - The race condition in start_pipeline.sh showed that re-checking dynamic state (worker counts) at different times can produce inconsistent results
6. **Single source of truth** - Resource decisions should be made once, not re-validated by multiple checks that can disagree

## Prevention Checklist

Before starting any DeSurv pipeline:

- [ ] Wait for other pipelines (baton, etc.) to complete OR verify <60% memory in use
- [ ] Run `./scripts/start_pipeline.sh` (enforces preflight, auto-starts watchdog)
- [ ] If preflight FAILS on resource check, DO NOT use --force unless you understand the risk
- [ ] Monitor pipeline: check watchdog.log or run `./scripts/watchdog.sh` manually
