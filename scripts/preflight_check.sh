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

if [ "$mem_pct" -lt 20 ]; then
    check_warn "Low available memory: ${mem_pct}% free"
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
