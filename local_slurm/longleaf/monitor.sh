#!/bin/bash
# monitor.sh - Pipeline monitoring helper for DeSurv on Longleaf
#
# Usage:
#   bash local_slurm/longleaf/monitor.sh              # Summary (queue + progress counts)
#   bash local_slurm/longleaf/monitor.sh --progress    # Full tar_progress() output
#   bash local_slurm/longleaf/monitor.sh --log         # Tail pipeline log
#   bash local_slurm/longleaf/monitor.sh --workers     # Show active crew worker jobs
#   bash local_slurm/longleaf/monitor.sh --errors      # Show errored targets
#   bash local_slurm/longleaf/monitor.sh --help        # This help

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

# Ensure R is available
load_r() {
    if ! command -v Rscript &> /dev/null; then
        module load r/4.4.0 2>/dev/null || true
    fi
    if ! command -v Rscript &> /dev/null; then
        echo -e "${RED}Error: Rscript not available. Run: module load r/4.4.0${NC}"
        exit 1
    fi
}

show_help() {
    cat <<'EOF'
Usage: bash local_slurm/longleaf/monitor.sh [OPTION]

Monitor DeSurv pipeline progress on Longleaf.

Options:
  (no args)     Summary: Slurm queue + tar_progress() counts
  --progress    Full tar_progress() (all non-skipped targets)
  --log         Tail logs/_targets.out (last 50 lines)
  --workers     Show active crew worker Slurm jobs
  --errors      Show errored targets with messages
  --help        Show this help

Examples:
  monitor.sh                # Quick status check
  monitor.sh --progress     # Detailed target status
  watch -n 30 monitor.sh    # Auto-refresh every 30 seconds
EOF
}

# =========================================================================
# Default: summary view
# =========================================================================
cmd_summary() {
    echo ""
    echo -e "${BOLD}=== DeSurv Pipeline Monitor ===${NC}"
    echo "$(date)"
    echo ""

    # Slurm queue
    echo -e "${BOLD}Slurm Queue${NC}"
    echo "─────────────────────────────────────────"
    my_jobs=$(squeue -u "$USER" -h -o "%-10i %-20j %-8T %-10M %-6C %R" 2>/dev/null || true)
    if [ -n "$my_jobs" ]; then
        echo "JOBID      NAME                 STATE    TIME       CPUS   REASON"
        echo "$my_jobs"
    else
        echo "  (no jobs in queue)"
    fi
    echo ""

    # tar_progress summary
    echo -e "${BOLD}Target Progress${NC}"
    echo "─────────────────────────────────────────"
    load_r
    Rscript --vanilla -e "
suppressMessages(library(targets))
tryCatch({
  p <- tar_progress()
  counts <- table(p\$progress)
  for (status in names(counts)) {
    cat(sprintf('  %-12s %d\n', status, counts[status]))
  }
  total <- nrow(p)
  active <- sum(p\$progress %in% c('dispatched', 'started'))
  done <- sum(p\$progress == 'completed')
  cat(sprintf('\n  Total: %d | Active: %d | Done: %d\n', total, active, done))
}, error = function(e) {
  cat('  (no targets store found or no progress yet)\n')
})
" 2>/dev/null
    echo ""
}

# =========================================================================
# --progress: full non-skipped target listing
# =========================================================================
cmd_progress() {
    echo ""
    echo -e "${BOLD}=== Target Progress (non-skipped) ===${NC}"
    echo ""
    load_r
    Rscript --vanilla -e "
suppressMessages(library(targets))
tryCatch({
  p <- tar_progress()
  active <- p[p\$progress != 'skipped', ]
  if (nrow(active) == 0) {
    cat('  All targets skipped (up to date).\n')
  } else {
    # Sort: dispatched/started first, then completed, then others
    priority <- c('started' = 1, 'dispatched' = 2, 'errored' = 3,
                   'canceled' = 4, 'completed' = 5)
    active\$sort_key <- ifelse(active\$progress %in% names(priority),
                              priority[active\$progress], 6)
    active <- active[order(active\$sort_key, active\$name), ]
    for (i in seq_len(nrow(active))) {
      status <- active\$progress[i]
      marker <- switch(status,
        started = '\033[0;32m>>',
        dispatched = '\033[1;33m..',
        errored = '\033[0;31m!!',
        completed = '\033[0;32mOK',
        '\033[0m--'
      )
      cat(sprintf('  %s\033[0m %-12s %s\n', marker, status, active\$name[i]))
    }
    cat(sprintf('\n  Showing %d non-skipped target(s).\n', nrow(active)))
  }
}, error = function(e) {
  cat('  Error:', conditionMessage(e), '\n')
})
" 2>/dev/null
    echo ""
}

# =========================================================================
# --log: tail pipeline log
# =========================================================================
cmd_log() {
    local logfile="$PROJECT_DIR/logs/_targets.out"
    echo ""
    echo -e "${BOLD}=== Pipeline Log (last 50 lines) ===${NC}"
    echo "File: $logfile"
    echo "─────────────────────────────────────────"
    if [ -f "$logfile" ]; then
        tail -n 50 "$logfile"
    else
        echo "  (log file not found — pipeline may not have started yet)"
    fi
    echo ""
}

# =========================================================================
# --workers: crew worker Slurm jobs
# =========================================================================
cmd_workers() {
    echo ""
    echo -e "${BOLD}=== Active Crew Workers ===${NC}"
    echo ""

    crew_jobs=$(squeue -u "$USER" -h -o "%-10i %-25j %-8T %-10M %-6C %-8m %R" 2>/dev/null | grep -i "crew" || true)
    if [ -n "$crew_jobs" ]; then
        echo "JOBID      NAME                      STATE    TIME       CPUS   MEM      REASON"
        echo "$crew_jobs"
        echo ""
        n_workers=$(echo "$crew_jobs" | wc -l)
        n_running=$(echo "$crew_jobs" | grep -c "RUNNING" || true)
        n_pending=$(echo "$crew_jobs" | grep -c "PENDING" || true)
        echo "  Total: $n_workers | Running: $n_running | Pending: $n_pending"
    else
        echo "  (no crew worker jobs found in Slurm queue)"
    fi
    echo ""
}

# =========================================================================
# --errors: errored targets with messages
# =========================================================================
cmd_errors() {
    echo ""
    echo -e "${BOLD}=== Errored Targets ===${NC}"
    echo ""
    load_r
    Rscript --vanilla -e "
suppressMessages(library(targets))
tryCatch({
  m <- tar_meta(fields = c('error', 'warnings'))
  errored <- m[!is.na(m\$error) & nzchar(m\$error), ]
  warned <- m[!is.na(m\$warnings) & nzchar(m\$warnings), ]

  if (nrow(errored) == 0 && nrow(warned) == 0) {
    cat('  No errors or warnings found.\n')
  } else {
    if (nrow(errored) > 0) {
      cat(sprintf('Errors (%d):\n', nrow(errored)))
      for (i in seq_len(nrow(errored))) {
        cat(sprintf('  \033[0;31m!!\033[0m %s\n     %s\n',
            errored\$name[i], errored\$error[i]))
      }
      cat('\n')
    }
    if (nrow(warned) > 0) {
      cat(sprintf('Warnings (%d):\n', nrow(warned)))
      for (i in seq_len(min(nrow(warned), 10))) {
        cat(sprintf('  \033[1;33m!!\033[0m %s\n     %s\n',
            warned\$name[i], warned\$warnings[i]))
      }
      if (nrow(warned) > 10) {
        cat(sprintf('  ... and %d more\n', nrow(warned) - 10))
      }
    }
  }
}, error = function(e) {
  cat('  Error reading store:', conditionMessage(e), '\n')
})
" 2>/dev/null
    echo ""
}

# =========================================================================
# Main dispatch
# =========================================================================
case "${1:-}" in
    --progress|-p) cmd_progress ;;
    --log|-l)      cmd_log ;;
    --workers|-w)  cmd_workers ;;
    --errors|-e)   cmd_errors ;;
    --help|-h)     show_help ;;
    "")            cmd_summary ;;
    *)
        echo "Unknown option: $1"
        show_help
        exit 1
        ;;
esac
