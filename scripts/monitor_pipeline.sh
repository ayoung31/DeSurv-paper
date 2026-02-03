#!/usr/bin/env bash
# Monitor pipeline resource usage
# Phase 1: every 5 min for first 30 min
# Phase 2: every 15 min after that
# Usage: ./scripts/monitor_pipeline.sh

set -euo pipefail

LOG_DIR="logs"
MONITOR_LOG="${LOG_DIR}/monitor_$(date +%Y%m%d_%H%M%S).log"
PIPELINE_LOG="${LOG_DIR}/pipeline_restart.log"
PHASE1_INTERVAL=300   # 5 min in seconds
PHASE2_INTERVAL=900   # 15 min in seconds
PHASE1_DURATION=1800  # 30 min in seconds

START_TIME=$(date +%s)

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MONITOR_LOG"
}

check_status() {
  log_msg "========== STATUS CHECK =========="

  # Memory
  mem_info=$(free -m | awk '/^Mem:/ {printf "Used: %dMB / %dMB (%.0f%%), Available: %dMB", $3, $2, $3/$2*100, $7}')
  swap_info=$(free -m | awk '/^Swap:/ {printf "Swap: %dMB / %dMB", $3, $2}')
  log_msg "MEMORY: $mem_info | $swap_info"

  # CPU load
  load=$(uptime | awk -F'load average:' '{print $2}')
  log_msg "CPU LOAD:$load"

  # Pipeline process
  pipe_procs=$(pgrep -f "tar_make|crew::crew_worker" 2>/dev/null | wc -l)
  log_msg "PROCESSES: $pipe_procs (tar_make + crew workers)"

  # Crew workers by controller (from cmdline args)
  if pgrep -f "crew::crew_worker" > /dev/null 2>&1; then
    worker_count=$(pgrep -cf "crew::crew_worker" 2>/dev/null || echo 0)
    worker_rss=$(ps -C R -C Rscript -o rss= 2>/dev/null | awk '{sum+=$1} END {printf "%.0f", sum/1024}')
    log_msg "CREW WORKERS: $worker_count active, ${worker_rss}MB total RSS"
  fi

  # Latest pipeline log lines (last completed/dispatched target)
  if [ -f "$PIPELINE_LOG" ]; then
    last_target=$(grep -E "^[✔✖+]" "$PIPELINE_LOG" 2>/dev/null | tail -3)
    if [ -n "$last_target" ]; then
      log_msg "LAST TARGETS:"
      echo "$last_target" | while read -r line; do
        log_msg "  $line"
      done
    fi

    # Count completed vs errored
    completed=$(grep -c "^✔" "$PIPELINE_LOG" 2>/dev/null || echo 0)
    errored=$(grep -c "^✖" "$PIPELINE_LOG" 2>/dev/null || echo 0)
    dispatched=$(grep -c "^+" "$PIPELINE_LOG" 2>/dev/null || echo 0)
    log_msg "PROGRESS: $completed completed, $errored errored, $dispatched dispatched"
  fi

  # Check if pipeline is still running
  if ! pgrep -f "tar_make" > /dev/null 2>&1; then
    log_msg "WARNING: Pipeline process not found - may have finished or crashed"
    if [ -f "$PIPELINE_LOG" ]; then
      log_msg "FINAL LOG LINES:"
      tail -5 "$PIPELINE_LOG" | while read -r line; do
        log_msg "  $line"
      done
    fi
    return 1
  fi

  log_msg "---"
  return 0
}

log_msg "Pipeline monitor started"
log_msg "Phase 1: every ${PHASE1_INTERVAL}s for ${PHASE1_DURATION}s"
log_msg "Phase 2: every ${PHASE2_INTERVAL}s after that"
log_msg "Monitor log: $MONITOR_LOG"

# Initial check
check_status || true

while true; do
  elapsed=$(( $(date +%s) - START_TIME ))

  if [ "$elapsed" -lt "$PHASE1_DURATION" ]; then
    sleep "$PHASE1_INTERVAL"
  else
    sleep "$PHASE2_INTERVAL"
  fi

  if ! check_status; then
    log_msg "Pipeline appears to have stopped. Exiting monitor."
    exit 0
  fi
done
