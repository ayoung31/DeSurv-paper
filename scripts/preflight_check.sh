#!/bin/bash
# DeSurv Pre-Submission Preflight Check
# Run before submitting jobs: ./scripts/preflight_check.sh
# Or use as wrapper: ./scripts/preflight_check.sh sbatch _targets.sh
#
# Checks for:
# 1. Orphaned R/crew processes
# 2. Stale Slurm jobs from this project
# 3. Crew backend state issues
# 4. Targets store lock files
# 5. System resource availability
# 6. Failed jobs that may have left state

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"

errors=0
warnings=0

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║           DeSurv Pre-Submission Preflight Check              ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Project: $PROJECT_DIR"
echo "Time: $(date)"
echo ""

check_pass() {
    echo -e "   ${GREEN}✓${NC} $1"
}

check_warn() {
    echo -e "   ${YELLOW}⚠${NC} $1"
    ((warnings++)) || true
}

check_fail() {
    echo -e "   ${RED}✗${NC} $1"
    ((errors++)) || true
}

# =============================================================================
# 1. Orphaned R Processes
# =============================================================================
echo "1. Checking for orphaned R processes..."
echo "   ─────────────────────────────────────"

# Find R processes in this project directory using elapsed time (etime)
orphaned_r=""
while read pid etime cmd; do
    [ -z "$pid" ] && continue
    cwd=$(readlink /proc/$pid/cwd 2>/dev/null || echo "")
    if [[ "$cwd" == "$PROJECT_DIR"* ]]; then
        orphaned_r="$orphaned_r$pid $etime $cmd\n"
    fi
done < <(ps -eo pid,etime,args | grep -E "[R]script|[R].*--no-echo" | grep -v grep)

if [ -n "$orphaned_r" ]; then
    count=$(echo -e "$orphaned_r" | grep -c . || echo 0)
    check_warn "Found $count R process(es) in project directory"
    echo -e "$orphaned_r" | while read pid etime cmd; do
        [ -z "$pid" ] && continue
        # Parse elapsed time to detect long-running processes (format: [[DD-]HH:]MM:SS)
        hours=0
        if [[ "$etime" =~ ^([0-9]+)-([0-9]+):([0-9]+):([0-9]+)$ ]]; then
            # DD-HH:MM:SS format
            hours=$(( ${BASH_REMATCH[1]} * 24 + ${BASH_REMATCH[2]} ))
        elif [[ "$etime" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
            # HH:MM:SS format
            hours=${BASH_REMATCH[1]}
        fi

        age_note=""
        if [ "$hours" -ge 12 ]; then
            age_note=" [LONG-RUNNING]"
        fi
        echo "         PID $pid (elapsed: $etime)$age_note: ${cmd:0:50}..."
    done
    pids_only=$(echo -e "$orphaned_r" | awk '{print $1}' | tr '\n' ' ')
    echo ""
    echo "         To kill: kill $pids_only"
else
    check_pass "No orphaned R processes in project directory"
fi

# =============================================================================
# 2. Orphaned Crew Workers
# =============================================================================
echo ""
echo "2. Checking for orphaned crew workers..."
echo "   ─────────────────────────────────────"

# Find crew workers in this project
crew_workers=$(pgrep -f "crew::crew_worker" 2>/dev/null || true)
orphaned_crew=""

if [ -n "$crew_workers" ]; then
    for pid in $crew_workers; do
        cwd=$(readlink /proc/$pid/cwd 2>/dev/null || echo "")
        if [[ "$cwd" == "$PROJECT_DIR"* ]]; then
            orphaned_crew="$orphaned_crew $pid"
        fi
    done
fi

if [ -n "$orphaned_crew" ]; then
    count=$(echo $orphaned_crew | wc -w)
    check_warn "Found $count crew worker(s) in project directory"
    echo "         PIDs:$orphaned_crew"
    echo "         To kill: kill$orphaned_crew"
else
    check_pass "No orphaned crew workers in project directory"
fi

# =============================================================================
# 2b. Cross-Project Crew Workers (Resource Contention)
# =============================================================================
echo ""
echo "2b. Checking for cross-project crew workers..."
echo "    ──────────────────────────────────────────"

# Find ALL crew workers regardless of project
cross_project_pids=""
cross_project_info=""
cross_project_jobs=""

if [ -n "$crew_workers" ]; then
    for pid in $crew_workers; do
        cwd=$(readlink /proc/$pid/cwd 2>/dev/null || echo "")
        if [[ -n "$cwd" && "$cwd" != "$PROJECT_DIR"* ]]; then
            # Get elapsed time
            etime=$(ps -o etime= -p $pid 2>/dev/null | tr -d ' ' || echo "unknown")
            proj_name=$(basename "$cwd")
            cross_project_pids="$cross_project_pids $pid"
            cross_project_info="$cross_project_info\n         PID $pid ($etime): $proj_name"
        fi
    done
fi

# Also check for Slurm crew-worker jobs from other projects
if command -v squeue &> /dev/null; then
    while IFS='|' read jobid name workdir; do
        [ -z "$jobid" ] && continue
        if [[ "$name" == crew-worker* && "$workdir" != "$PROJECT_DIR"* ]]; then
            proj_name=$(basename "$workdir")
            cross_project_jobs="$cross_project_jobs $jobid"
            cross_project_info="$cross_project_info\n         Job $jobid: $proj_name ($name)"
        fi
    done < <(squeue -u $USER -h -o "%i|%j|%Z" 2>/dev/null | grep -i "crew-worker" || true)
fi

if [ -n "$cross_project_pids" ] || [ -n "$cross_project_jobs" ]; then
    check_warn "Found crew workers from OTHER projects (may cause resource contention):"
    echo -e "$cross_project_info"
    echo ""
    if [ -n "$cross_project_pids" ]; then
        pids_trimmed=$(echo $cross_project_pids | xargs)
        echo "         To kill processes: kill $pids_trimmed"
    fi
    if [ -n "$cross_project_jobs" ]; then
        jobs_trimmed=$(echo $cross_project_jobs | xargs)
        echo "         To cancel Slurm jobs: scancel $jobs_trimmed"
    fi
else
    check_pass "No cross-project crew workers detected"
fi

# =============================================================================
# 3. Stale Slurm Jobs
# =============================================================================
echo ""
echo "3. Checking for stale/failed Slurm jobs..."
echo "   ────────────────────────────────────────"

# Check for running jobs from this project
running_jobs=$(squeue -u $USER -h -o "%i %j %T" 2>/dev/null | grep -E "desurv|DeSurv" || true)

if [ -n "$running_jobs" ]; then
    check_warn "Found running DeSurv jobs:"
    echo "$running_jobs" | while read line; do
        echo "         $line"
    done
else
    check_pass "No running DeSurv Slurm jobs"
fi

# Check for recently failed jobs (last 24 hours)
if command -v sacct &> /dev/null; then
    failed_jobs=$(sacct -u $USER --starttime=$(date -d '24 hours ago' '+%Y-%m-%dT%H:%M:%S') \
        --format=JobID,JobName,State,ExitCode,End -P 2>/dev/null | \
        grep -E "desurv|DeSurv" | grep -E "FAILED|CANCELLED|TIMEOUT" | head -5 || true)

    if [ -n "$failed_jobs" ]; then
        check_warn "Recently failed/cancelled DeSurv jobs (last 24h):"
        echo "$failed_jobs" | while IFS='|' read jobid name state exit end; do
            echo "         $jobid: $state (exit $exit) at $end"
        done
    else
        check_pass "No recently failed DeSurv jobs"
    fi
else
    check_pass "sacct not available - skipping failure history check"
fi

# =============================================================================
# 4. Crew Backend State
# =============================================================================
echo ""
echo "4. Checking crew backend state..."
echo "   ───────────────────────────────"

# Check for crew lock files
crew_locks=$(find "$PROJECT_DIR" -name "*.lock" -o -name "crew-*" -type f 2>/dev/null | head -10)

if [ -n "$crew_locks" ]; then
    count=$(echo "$crew_locks" | wc -l)
    check_warn "Found $count potential crew state file(s):"
    echo "$crew_locks" | while read f; do
        echo "         $f"
    done
else
    check_pass "No crew lock files found"
fi

# Check for mirai/nanonext socket files
socket_files=$(find /tmp -name "*nanonext*" -o -name "*mirai*" -user $USER 2>/dev/null | head -5 || true)

if [ -n "$socket_files" ]; then
    count=$(echo "$socket_files" | wc -l)
    check_warn "Found $count socket file(s) in /tmp (may be from active sessions):"
    echo "$socket_files" | head -3 | while read f; do
        echo "         $f"
    done
else
    check_pass "No stale socket files in /tmp"
fi

# =============================================================================
# 5. Targets Store Lock
# =============================================================================
echo ""
echo "5. Checking targets store state..."
echo "   ────────────────────────────────"

# Find store directories
stores=$(find "$PROJECT_DIR" -maxdepth 1 -type d -name "store_*" 2>/dev/null)

if [ -n "$stores" ]; then
    for store in $stores; do
        store_name=$(basename "$store")

        # Check for lock file
        if [ -f "$store/.lock" ]; then
            # Check if lock is stale (process doesn't exist)
            lock_pid=$(cat "$store/.lock" 2>/dev/null || echo "")
            if [ -n "$lock_pid" ] && ! kill -0 "$lock_pid" 2>/dev/null; then
                check_fail "$store_name: STALE LOCK (PID $lock_pid doesn't exist)"
                echo "         To fix: rm $store/.lock"
            else
                check_warn "$store_name: locked by PID $lock_pid"
            fi
        else
            check_pass "$store_name: not locked"
        fi

        # Check for errored targets
        if [ -f "$store/meta/meta" ]; then
            # Quick check via R (extract just the number, removing R's "[1] " prefix)
            errored=$(Rscript --vanilla -e "
                suppressMessages(library(targets))
                tar_config_set(store = '$store')
                meta <- tar_meta()
                cat(sum(!is.na(meta\$error) & nzchar(meta\$error)))
            " 2>/dev/null || echo "0")

            if [ "$errored" -gt 0 ] 2>/dev/null; then
                check_warn "$store_name: $errored errored target(s)"
            fi
        fi
    done
else
    check_pass "No store directories found (will be created)"
fi

# =============================================================================
# 6. System Resources
# =============================================================================
echo ""
echo "6. Checking system resources..."
echo "   ─────────────────────────────"

# Check load average
load=$(cat /proc/loadavg | awk '{print $1}')
cpus=$(nproc)
load_pct=$(echo "$load $cpus" | awk '{printf "%.0f", ($1/$2)*100}')

if [ "$load_pct" -gt 80 ]; then
    check_warn "High system load: $load ($load_pct% of $cpus CPUs)"
elif [ "$load_pct" -gt 50 ]; then
    check_warn "Moderate system load: $load ($load_pct% of $cpus CPUs)"
else
    check_pass "System load OK: $load ($load_pct% of $cpus CPUs)"
fi

# Check available memory
mem_available=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
mem_total=$(grep MemTotal /proc/meminfo | awk '{print $2}')
mem_pct=$(echo "$mem_available $mem_total" | awk '{printf "%.0f", ($1/$2)*100}')

if [ "$mem_pct" -lt 25 ]; then
    check_fail "CRITICAL: Very low memory: ${mem_pct}% free (need 25%+)"
elif [ "$mem_pct" -lt 40 ]; then
    check_warn "Low available memory: ${mem_pct}% free (recommend 40%+ for BO)"
else
    check_pass "Memory OK: ${mem_pct}% available"
fi

# Check disk space
disk_pct=$(df "$PROJECT_DIR" | tail -1 | awk '{print $5}' | tr -d '%')

if [ "$disk_pct" -gt 90 ]; then
    check_fail "Low disk space: ${disk_pct}% used"
elif [ "$disk_pct" -gt 80 ]; then
    check_warn "Disk space warning: ${disk_pct}% used"
else
    check_pass "Disk space OK: ${disk_pct}% used"
fi

# Check combined crew worker resource usage (all projects)
echo ""
echo "   Combined crew worker analysis:"
total_crew_count=0
total_crew_mem_kb=0
for pid in $(pgrep -f "crew::crew_worker" 2>/dev/null); do
    mem=$(ps -o rss= -p $pid 2>/dev/null || echo 0)
    total_crew_mem_kb=$((total_crew_mem_kb + mem))
    total_crew_count=$((total_crew_count + 1))
done

if [ "$total_crew_count" -gt 0 ]; then
    total_crew_mem_gb=$(echo "scale=1; $total_crew_mem_kb / 1048576" | bc 2>/dev/null || echo "?")
    mem_total_gb=$(echo "scale=1; $mem_total / 1048576" | bc 2>/dev/null || echo "?")

    # Estimate if adding DeSurv would exceed threshold (assume ~15GB for BO)
    expected_new_gb=15
    combined_gb=$(echo "scale=1; $total_crew_mem_kb / 1048576 + $expected_new_gb" | bc 2>/dev/null || echo "999")
    threshold_gb=$(echo "scale=1; $mem_total / 1048576 * 0.75" | bc 2>/dev/null || echo "0")

    if (( $(echo "$combined_gb > $threshold_gb" | bc -l 2>/dev/null || echo 1) )); then
        check_fail "CRITICAL: System-wide crew workers: $total_crew_count using ${total_crew_mem_gb}GB"
        echo "         Combined with new pipeline (~${expected_new_gb}GB) would be ~${combined_gb}GB"
        echo "         This EXCEEDS 75% threshold (${threshold_gb}GB)"
        echo "         BLOCKING: Wait for other pipelines to complete or use --force"
    else
        check_pass "System-wide crew workers: $total_crew_count using ${total_crew_mem_gb}GB (headroom OK)"
    fi
else
    check_pass "No crew workers running system-wide"
fi

# =============================================================================
# 7. Slurm Cluster Status
# =============================================================================
echo ""
echo "7. Checking Slurm cluster status..."
echo "   ──────────────────────────────────"

slurm_status=$(sinfo -h -o "%P %a %D %T" 2>/dev/null | head -1 || echo "")

if [ -z "$slurm_status" ]; then
    check_fail "Cannot connect to Slurm controller"
else
    partition=$(echo "$slurm_status" | awk '{print $1}')
    avail=$(echo "$slurm_status" | awk '{print $2}')
    nodes=$(echo "$slurm_status" | awk '{print $3}')
    state=$(echo "$slurm_status" | awk '{print $4}')

    if [ "$avail" == "up" ]; then
        check_pass "Slurm partition $partition: $avail ($nodes node(s), $state)"
    else
        check_fail "Slurm partition $partition: $avail"
    fi
fi

# Pending jobs count
pending=$(squeue -u $USER -h -t PD 2>/dev/null | wc -l)
running=$(squeue -u $USER -h -t R 2>/dev/null | wc -l)

check_pass "Your jobs: $running running, $pending pending"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
printf "║  Summary: %d error(s), %d warning(s)                           ║\n" "$errors" "$warnings"
echo "╚══════════════════════════════════════════════════════════════╝"

if [ $errors -gt 0 ]; then
    echo ""
    echo -e "${RED}✗ PREFLIGHT FAILED${NC} - Fix errors before submitting jobs"
    echo ""
    exit 1
elif [ $warnings -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}⚠ PREFLIGHT PASSED WITH WARNINGS${NC} - Review before proceeding"
    echo ""
fi

# If arguments provided, execute them (wrapper mode)
if [ $# -gt 0 ]; then
    if [ $errors -eq 0 ]; then
        echo "Executing: $@"
        echo ""
        exec "$@"
    fi
else
    if [ $errors -eq 0 ]; then
        echo -e "${GREEN}✓ Ready to submit jobs${NC}"
        echo ""
    fi
fi
