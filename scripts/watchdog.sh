#!/bin/bash
# =============================================================================
# WATCHDOG: Monitor running DeSurv jobs for hangs/failures
# =============================================================================
# Checks for:
#   - Jobs running but not producing output (hung)
#   - Orphaned R processes after Slurm jobs end
#   - Targets controller died while workers running
#   - Log files not updating (stalled progress)
#   - Store not receiving new results
#
# Usage:
#   ./scripts/watchdog.sh              # One-time check
#   ./scripts/watchdog.sh --watch      # Continuous monitoring (every 5 min)
#   ./scripts/watchdog.sh --watch 120  # Custom interval (seconds)
#   ./scripts/watchdog.sh --alert      # Send email on issues
#
# =============================================================================

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_NAME="$(basename "$PROJECT_DIR")"
LOG_DIR="$PROJECT_DIR/logs"
WATCHDOG_LOG="$LOG_DIR/watchdog.log"
STORE_PATTERN="store_PKG_VERSION=*"
LAST_FAILURE_FILE="/tmp/desurv_watchdog_last_failures_$USER"

# Thresholds
LOG_STALE_MINUTES=30        # Log file not updated for this long = stale
STORE_STALE_MINUTES=60      # No new store objects for this long = stale
WORKER_ORPHAN_MINUTES=10    # Workers running without controller for this long
JOB_LONG_RUNNING_HOURS=24   # Jobs running longer than this trigger warning

# Adaptive interval settings
INTERVAL_FAST=300           # 5 minutes when issues detected
INTERVAL_SLOW=900           # 15 minutes when stable
consecutive_ok=0            # Track consecutive OK checks

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# State
issues=0
alerts=""

# Parse arguments
WATCH_MODE=false
WATCH_INTERVAL=300
ALERT_MODE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --watch|-w)
            WATCH_MODE=true
            if [[ $2 =~ ^[0-9]+$ ]]; then
                WATCH_INTERVAL=$2
                shift
            fi
            shift
            ;;
        --alert|-a)
            ALERT_MODE=true
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# =============================================================================
# Helper Functions
# =============================================================================

log_to_file() {
    if [ -n "$WATCHDOG_LOG" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$WATCHDOG_LOG"
    fi
}

log_info() {
    echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $1"
    log_to_file "INFO: $1"
}

log_ok() {
    echo -e "${GREEN}[$(date '+%H:%M:%S')] ✓${NC} $1"
    log_to_file "OK: $1"
}

log_warn() {
    echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠${NC} $1"
    ((issues++)) || true
    alerts="$alerts\n⚠ $1"
    log_to_file "WARN: $1"
}

log_error() {
    echo -e "${RED}[$(date '+%H:%M:%S')] ✗${NC} $1"
    ((issues++)) || true
    alerts="$alerts\n✗ $1"
    log_to_file "ERROR: $1"
}

get_file_age_minutes() {
    local file=$1
    if [ -f "$file" ]; then
        local now=$(date +%s)
        local mtime=$(stat -c %Y "$file" 2>/dev/null || echo $now)
        echo $(( (now - mtime) / 60 ))
    else
        echo 9999
    fi
}

# =============================================================================
# Check Functions
# =============================================================================

check_slurm_jobs() {
    log_info "Checking Slurm jobs..."

    # Get DeSurv-related jobs
    local jobs=$(squeue -u $USER -h -o "%i|%j|%T|%M|%N" 2>/dev/null | grep -iE "desurv|crew" || true)

    if [ -z "$jobs" ]; then
        log_ok "No DeSurv jobs currently running"
        return
    fi

    echo "$jobs" | while IFS='|' read jobid name state time node; do
        if [ "$state" == "RUNNING" ]; then
            # Parse elapsed time (format: D-HH:MM:SS or HH:MM:SS or MM:SS)
            local hours=0
            if [[ "$time" =~ ^([0-9]+)-([0-9]+): ]]; then
                # D-HH:MM:SS format
                hours=$(( ${BASH_REMATCH[1]} * 24 + ${BASH_REMATCH[2]} ))
            elif [[ "$time" =~ ^([0-9]+):([0-9]+):([0-9]+)$ ]]; then
                # HH:MM:SS format
                hours=${BASH_REMATCH[1]}
            fi

            # Check for long-running jobs (> 24 hours)
            if [ "$hours" -ge "$JOB_LONG_RUNNING_HOURS" ]; then
                log_warn "Job $jobid ($name) running for ${hours}+ hours - may be stuck"
            fi

            # Check if this is a main controller or worker
            if [[ "$name" == *"crew"* ]]; then
                log_ok "Worker $jobid ($name) running for $time"
            else
                log_ok "Controller $jobid ($name) running for $time"

                # For controllers, check if log is updating
                local log_file="$LOG_DIR/slurm-${jobid}.out"
                if [ -f "$log_file" ]; then
                    local age=$(get_file_age_minutes "$log_file")
                    if [ "$age" -gt "$LOG_STALE_MINUTES" ]; then
                        log_warn "Job $jobid log stale for ${age}m (threshold: ${LOG_STALE_MINUTES}m)"
                    fi
                fi
            fi
        elif [ "$state" == "PENDING" ]; then
            log_info "Job $jobid ($name) pending"
        fi
    done
}

check_controller_worker_sync() {
    log_info "Checking controller-worker synchronization..."

    # Find controller jobs (main pipeline jobs, not crew workers)
    local controllers=$(squeue -u $USER -h -o "%i %j" 2>/dev/null | grep -iE "desurv_pipeline|desurv_sims" | grep -v crew || true)
    local workers=$(squeue -u $USER -h -o "%i %j" 2>/dev/null | grep -i "crew" || true)

    local controller_count=$(echo "$controllers" | grep -c . || echo 0)
    local worker_count=$(echo "$workers" | grep -c . || echo 0)

    if [ "$worker_count" -gt 0 ] && [ "$controller_count" -eq 0 ]; then
        log_error "ORPHANED WORKERS: $worker_count crew workers running without controller!"
        echo "$workers" | head -5 | while read jobid name; do
            echo "         $jobid: $name"
        done
        echo ""
        echo "         To cancel: scancel $(echo "$workers" | awk '{print $1}' | tr '\n' ' ')"
    elif [ "$controller_count" -gt 0 ] && [ "$worker_count" -eq 0 ]; then
        log_warn "Controller running but no workers (may be waiting for Slurm)"
    elif [ "$controller_count" -gt 0 ] && [ "$worker_count" -gt 0 ]; then
        log_ok "Controller(s) and workers in sync ($controller_count controllers, $worker_count workers)"
    else
        log_ok "No active controller-worker pairs"
    fi
}

check_orphaned_processes() {
    log_info "Checking for orphaned processes..."

    # Find R processes in project directory
    local orphaned=0

    for pid in $(pgrep -f "R.*--no-echo" 2>/dev/null || true); do
        local cwd=$(readlink /proc/$pid/cwd 2>/dev/null || echo "")
        if [[ "$cwd" == "$PROJECT_DIR"* ]]; then
            # Check if this process is part of an active Slurm job
            local in_slurm=false
            local job_pids=$(squeue -u $USER -h -o "%i" 2>/dev/null | while read jid; do
                scontrol listpids $jid 2>/dev/null | awk '{print $1}'
            done || true)

            if echo "$job_pids" | grep -q "^${pid}$"; then
                in_slurm=true
            fi

            if [ "$in_slurm" = false ]; then
                ((orphaned++)) || true
                local cmd=$(ps -p $pid -o args= 2>/dev/null | head -c 60)
                log_warn "Orphaned process PID $pid: $cmd..."
            fi
        fi
    done

    if [ "$orphaned" -eq 0 ]; then
        log_ok "No orphaned R processes"
    fi
}

check_log_progress() {
    log_info "Checking log file progress..."

    # Find recent log files
    local recent_logs=$(find "$LOG_DIR" -name "slurm-*.out" -mmin -60 2>/dev/null | head -5)

    if [ -z "$recent_logs" ]; then
        log_info "No recent log files (last 60 min)"
        return
    fi

    for log in $recent_logs; do
        local jobid=$(basename "$log" | sed 's/slurm-\([0-9]*\).out/\1/')
        local age=$(get_file_age_minutes "$log")
        local size=$(stat -c %s "$log" 2>/dev/null || echo 0)
        local size_kb=$((size / 1024))

        # Check if job is still running
        local job_state=$(squeue -j "$jobid" -h -o "%T" 2>/dev/null || echo "COMPLETED")

        if [ "$job_state" == "RUNNING" ]; then
            if [ "$age" -gt "$LOG_STALE_MINUTES" ]; then
                # Check last few lines for activity indicators
                local last_lines=$(tail -5 "$log" 2>/dev/null)
                if echo "$last_lines" | grep -qE "waiting|dispatched|completed|running"; then
                    log_ok "Job $jobid active (${size_kb}KB, ${age}m since update)"
                else
                    log_warn "Job $jobid may be stalled (${size_kb}KB, ${age}m since update)"
                fi
            else
                log_ok "Job $jobid progressing (${size_kb}KB, updated ${age}m ago)"
            fi
        fi
    done
}

check_store_progress() {
    log_info "Checking targets store progress..."

    # Find the active store
    local stores=$(find "$PROJECT_DIR" -maxdepth 1 -type d -name "$STORE_PATTERN" 2>/dev/null)

    for store in $stores; do
        local store_name=$(basename "$store")
        local objects_dir="$store/objects"

        if [ ! -d "$objects_dir" ]; then
            continue
        fi

        # Find most recent object
        local newest=$(find "$objects_dir" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1)

        if [ -n "$newest" ]; then
            local newest_file=$(echo "$newest" | cut -d' ' -f2-)
            local age=$(get_file_age_minutes "$newest_file")
            local count=$(find "$objects_dir" -type f 2>/dev/null | wc -l)

            if [ "$age" -gt "$STORE_STALE_MINUTES" ]; then
                # Check if any jobs are running
                local running=$(squeue -u $USER -h -o "%T" 2>/dev/null | grep -c "RUNNING" || echo 0)
                if [ "$running" -gt 0 ]; then
                    log_warn "$store_name: No new objects in ${age}m ($count total) - jobs may be stalled"
                else
                    log_ok "$store_name: $count objects (last update ${age}m ago, no jobs running)"
                fi
            else
                log_ok "$store_name: Active ($count objects, last update ${age}m ago)"
            fi
        fi
    done
}

check_failed_jobs() {
    log_info "Checking for recently failed jobs..."

    # Skip if sacct is not available (some clusters don't have accounting)
    if ! command -v sacct &> /dev/null; then
        log_info "sacct not available - skipping job failure check"
        return
    fi

    # Check last 2 hours - capture both output and exit status
    local sacct_output
    local sacct_status
    sacct_output=$(sacct -u $USER --starttime=$(date -d '2 hours ago' '+%Y-%m-%dT%H:%M:%S') \
        --format=JobID,JobName,State,ExitCode,End -P 2>&1)
    sacct_status=$?

    # Handle sacct returning error
    if [ $sacct_status -ne 0 ]; then
        log_warn "sacct query failed (exit $sacct_status) - Slurm accounting may be unavailable"
        return
    fi

    # Filter for DeSurv failures
    local failed=$(echo "$sacct_output" | grep -iE "desurv" | grep -E "FAILED|TIMEOUT|OUT_OF_ME" | head -10 || true)

    if [ -n "$failed" ]; then
        # Check for NEW failures (not seen before)
        local current_failures=$(echo "$failed" | sort)
        local previous_failures=""
        if [ -f "$LAST_FAILURE_FILE" ]; then
            previous_failures=$(cat "$LAST_FAILURE_FILE" | sort)
        fi

        # Find new failures
        local new_failures=$(comm -23 <(echo "$current_failures") <(echo "$previous_failures") 2>/dev/null || echo "$current_failures")

        if [ -n "$new_failures" ] && [ "$new_failures" != "" ]; then
            log_error "NEW job failures detected:"
            echo "$new_failures" | while IFS='|' read jobid name state exit end; do
                echo "         $jobid ($name): $state, exit $exit at $end"
            done
        else
            log_warn "Previously seen failed jobs (last 2h):"
            echo "$failed" | head -3 | while IFS='|' read jobid name state exit end; do
                echo "         $jobid ($name): $state"
            done
        fi

        # Update cache
        echo "$current_failures" > "$LAST_FAILURE_FILE"
    else
        log_ok "No recently failed jobs"
        # Clear cache if no failures
        rm -f "$LAST_FAILURE_FILE" 2>/dev/null || true
    fi
}

check_resource_usage() {
    log_info "Checking resource usage..."

    # System load
    local load=$(cat /proc/loadavg | awk '{print $1}')
    local cpus=$(nproc)
    local load_pct=$(echo "$load $cpus" | awk '{printf "%.0f", ($1/$2)*100}')

    if [ "$load_pct" -gt 95 ]; then
        log_warn "System overloaded: ${load_pct}% CPU usage"
    elif [ "$load_pct" -gt 80 ]; then
        log_info "High load: ${load_pct}% CPU usage"
    else
        log_ok "CPU load OK: ${load_pct}%"
    fi

    # Memory
    local mem_available=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
    local mem_total=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    local mem_pct=$(echo "$mem_available $mem_total" | awk '{printf "%.0f", ($1/$2)*100}')

    if [ "$mem_pct" -lt 10 ]; then
        log_error "Critical memory: only ${mem_pct}% available"
    elif [ "$mem_pct" -lt 20 ]; then
        log_warn "Low memory: ${mem_pct}% available"
    else
        log_ok "Memory OK: ${mem_pct}% available"
    fi
}

# =============================================================================
# Alert Function
# =============================================================================

send_alert() {
    if [ "$ALERT_MODE" = true ] && [ "$issues" -gt 0 ]; then
        local subject="[DeSurv Watchdog] $issues issue(s) detected"
        local body="Watchdog detected issues at $(date):\n$alerts\n\nProject: $PROJECT_DIR"

        # Try to send email via mail command or msmtp
        if command -v mail &> /dev/null; then
            echo -e "$body" | mail -s "$subject" "$USER"
            log_info "Alert email sent"
        elif command -v msmtp &> /dev/null; then
            echo -e "Subject: $subject\n\n$body" | msmtp "$USER@localhost"
            log_info "Alert sent via msmtp"
        else
            log_warn "No mail command available for alerts"
        fi
    fi
}

# =============================================================================
# Main Execution
# =============================================================================

run_checks() {
    issues=0
    alerts=""

    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║           DeSurv Watchdog - $(date '+%Y-%m-%d %H:%M:%S')            ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Project: $PROJECT_DIR"
    echo ""

    check_slurm_jobs
    echo ""
    check_controller_worker_sync
    echo ""
    check_orphaned_processes
    echo ""
    check_log_progress
    echo ""
    check_store_progress
    echo ""
    check_failed_jobs
    echo ""
    check_resource_usage

    echo ""
    echo "═══════════════════════════════════════════════════════════════"
    if [ "$issues" -eq 0 ]; then
        echo -e "${GREEN}All checks passed - system healthy${NC}"
    else
        echo -e "${YELLOW}$issues issue(s) detected - review above${NC}"
    fi
    echo "═══════════════════════════════════════════════════════════════"
    echo ""

    send_alert
}

# Ensure log directory exists
mkdir -p "$LOG_DIR"

# Run once or in watch mode
if [ "$WATCH_MODE" = true ]; then
    log_info "Starting watchdog in continuous mode (adaptive interval: ${INTERVAL_FAST}s-${INTERVAL_SLOW}s)"
    log_info "Log file: $WATCHDOG_LOG"
    log_info "Press Ctrl+C to stop"
    echo ""

    while true; do
        run_checks

        # Adaptive interval: slow down if no issues, speed up if problems
        if [ "$issues" -eq 0 ]; then
            ((consecutive_ok++)) || true
            if [ "$consecutive_ok" -ge 3 ]; then
                WATCH_INTERVAL=$INTERVAL_SLOW
            fi
        else
            consecutive_ok=0
            WATCH_INTERVAL=$INTERVAL_FAST
        fi

        echo ""
        log_info "Next check in ${WATCH_INTERVAL}s ($([ $consecutive_ok -ge 3 ] && echo 'stable' || echo 'monitoring'))..."
        sleep $WATCH_INTERVAL
        clear
    done
else
    run_checks
fi
