#!/bin/bash
# =============================================================================
# start_pipeline.sh - Safe pipeline launcher with mandatory checks
# =============================================================================
# This wrapper ENFORCES preflight checks and AUTO-STARTS watchdog monitoring.
# Use this instead of directly running tar_make() for long-running pipelines.
#
# Usage:
#   ./scripts/start_pipeline.sh                    # Run main pipeline (excludes paper)
#   ./scripts/start_pipeline.sh --with-paper       # Run pipeline including paper rendering
#   ./scripts/start_pipeline.sh --sims             # Run simulation pipeline
#   ./scripts/start_pipeline.sh --script FILE.R   # Run custom script
#   ./scripts/start_pipeline.sh --no-watchdog     # Skip watchdog (not recommended)
#   ./scripts/start_pipeline.sh --force           # Skip preflight (DANGEROUS)
#
# Created: 2026-02-02 (after controller crash incident)
# =============================================================================

set -e

# Configuration
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PIPELINE_LOG="$LOG_DIR/pipeline_${TIMESTAMP}.log"
WATCHDOG_LOG="$LOG_DIR/watchdog.log"

# Defaults
SCRIPT="_targets.R"
RUN_WATCHDOG=true
FORCE_MODE=false
ALERT_MODE=false
INCLUDE_PAPER=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --sims)
            SCRIPT="_targets_sims_local.R"
            shift
            ;;
        --script)
            SCRIPT="$2"
            shift 2
            ;;
        --no-watchdog)
            RUN_WATCHDOG=false
            shift
            ;;
        --force)
            FORCE_MODE=true
            shift
            ;;
        --with-paper)
            INCLUDE_PAPER=true
            shift
            ;;
        --alert)
            ALERT_MODE=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --sims          Run simulation pipeline"
            echo "  --script FILE   Run custom targets script"
            echo "  --with-paper    Include paper rendering target (excluded by default)"
            echo "  --no-watchdog   Skip watchdog monitoring (not recommended)"
            echo "  --force         Skip preflight checks (DANGEROUS)"
            echo "  --alert         Enable email alerts on issues"
            echo "  -h, --help      Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           DeSurv Safe Pipeline Launcher                      ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Project: $PROJECT_DIR"
echo "Script:  $SCRIPT"
echo "Log:     $PIPELINE_LOG"
echo "Time:    $(date)"
echo ""

# =============================================================================
# Step 1: Preflight Checks (MANDATORY unless --force)
# =============================================================================
echo -e "${YELLOW}Step 1: Running preflight checks...${NC}"
echo ""

if [ "$FORCE_MODE" = true ]; then
    echo -e "${RED}⚠ WARNING: Skipping preflight checks (--force mode)${NC}"
    echo -e "${RED}  This is DANGEROUS and may cause resource conflicts!${NC}"
    echo ""
    sleep 3
else
    if ! ./scripts/preflight_check.sh; then
        echo ""
        echo -e "${RED}✗ Preflight checks FAILED${NC}"
        echo "  Fix the issues above before starting the pipeline."
        echo ""
        echo "  If you understand the risks, use --force to skip checks."
        exit 1
    fi
fi

# =============================================================================
# Step 2: Confirm Preflight Decision (no redundant re-check)
# =============================================================================
echo ""
echo -e "${YELLOW}Step 2: Confirming preflight resource decision...${NC}"
echo ""
# NOTE: We trust preflight's resource check. DO NOT re-check here!
# A race condition caused a crash on 2026-02-02 when Step 2 got different
# numbers than preflight and incorrectly overrode a warning.
# See: docs/INCIDENT_POSTMORTEM_20260202.md

echo -e "   ${GREEN}✓${NC} Preflight passed - proceeding with pipeline"
echo "   (Resource checks handled by preflight to avoid race conditions)"

# =============================================================================
# Step 3: Start Watchdog (unless --no-watchdog)
# =============================================================================
echo ""
WATCHDOG_PID=""

if [ "$RUN_WATCHDOG" = true ]; then
    echo -e "${YELLOW}Step 3: Starting watchdog monitor...${NC}"

    WATCHDOG_ARGS="--watch"
    if [ "$ALERT_MODE" = true ]; then
        WATCHDOG_ARGS="$WATCHDOG_ARGS --alert"
    fi

    nohup ./scripts/watchdog.sh $WATCHDOG_ARGS >> "$WATCHDOG_LOG" 2>&1 &
    WATCHDOG_PID=$!

    sleep 2
    if kill -0 $WATCHDOG_PID 2>/dev/null; then
        echo -e "   ${GREEN}✓${NC} Watchdog started (PID $WATCHDOG_PID)"
        echo "   Log: $WATCHDOG_LOG"
    else
        echo -e "   ${YELLOW}⚠${NC} Watchdog failed to start - continuing anyway"
        WATCHDOG_PID=""
    fi
else
    echo -e "${YELLOW}Step 3: Watchdog SKIPPED (--no-watchdog)${NC}"
    echo -e "   ${YELLOW}⚠${NC} No automatic monitoring - check manually!"
fi

# =============================================================================
# Step 3.5: Clear stale progress from crashed runs
# =============================================================================
echo ""
echo -e "${YELLOW}Step 3.5: Checking for stale progress entries...${NC}"

Rscript --vanilla -e "
  suppressMessages(library(targets))
  prog <- targets::tar_progress()
  stuck <- sum(prog\$progress %in% c('dispatched', 'started'))
  if (stuck > 0) {
    cat('Found', stuck, 'targets stuck in dispatched/started state. Clearing progress...\n')
    targets::tar_destroy(destroy = 'progress')
    cat('Progress cleared. These targets will be re-dispatched.\n')
  } else {
    cat('No stale progress entries.\n')
  }
" 2>&1

# =============================================================================
# Step 4: Start Pipeline
# =============================================================================
echo ""
echo -e "${YELLOW}Step 4: Starting pipeline...${NC}"
echo ""
echo "   Script: $SCRIPT"
if [ "$INCLUDE_PAPER" = true ]; then
    echo "   Paper:  included (--with-paper)"
else
    echo "   Paper:  excluded (use --with-paper to include)"
fi
echo "   Log: $PIPELINE_LOG"
echo ""
echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
echo ""

# Record start info
START_TIME=$(date +%s)
echo "Pipeline started: $(date)" >> "$PIPELINE_LOG"
echo "Script: $SCRIPT" >> "$PIPELINE_LOG"
echo "Watchdog PID: ${WATCHDOG_PID:-none}" >> "$PIPELINE_LOG"
echo "---" >> "$PIPELINE_LOG"

# Run the pipeline
set +e  # Don't exit on error - we need to cleanup
if [ "$INCLUDE_PAPER" = true ]; then
    Rscript -e "targets::tar_make(script = '$SCRIPT')" 2>&1 | tee -a "$PIPELINE_LOG"
else
    Rscript -e "targets::tar_make(names = !tidyselect::any_of('paper'), script = '$SCRIPT')" 2>&1 | tee -a "$PIPELINE_LOG"
fi
EXIT_CODE=${PIPESTATUS[0]}
set -e

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo -e "${BLUE}════════════════════════════════════════════════════════════════${NC}"
echo ""

# =============================================================================
# Step 5: Cleanup
# =============================================================================
echo -e "${YELLOW}Step 5: Cleanup...${NC}"

# Stop watchdog
if [ -n "$WATCHDOG_PID" ] && kill -0 $WATCHDOG_PID 2>/dev/null; then
    kill $WATCHDOG_PID 2>/dev/null || true
    echo -e "   ${GREEN}✓${NC} Watchdog stopped"
fi

# Report
echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  Pipeline completed successfully                             ║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════╝${NC}"
else
    echo -e "${RED}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║  Pipeline FAILED (exit code $EXIT_CODE)                              ║${NC}"
    echo -e "${RED}╚══════════════════════════════════════════════════════════════╝${NC}"
fi

echo ""
echo "Duration: ${DURATION_MIN} minutes"
echo "Log: $PIPELINE_LOG"
echo ""

exit $EXIT_CODE
