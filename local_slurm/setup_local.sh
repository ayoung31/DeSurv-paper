#!/bin/bash
# setup_local.sh - Set up local configuration for DeSurv pipeline
#
# This script creates symlinks from local_slurm configs to the main directory,
# allowing you to run the pipeline locally without modifying the student's files.
#
# Supports three modes:
#   --quick (default): Fast testing with reduced iterations
#   --full: HPC-quality iterations, capped at 19 CPUs
#   --longleaf: UNC Longleaf HPC with large worker pools
#
# Usage:
#   ./local_slurm/setup_local.sh              # Quick mode (default)
#   ./local_slurm/setup_local.sh --quick      # Quick mode (explicit)
#   ./local_slurm/setup_local.sh --full       # Full quality mode
#   ./local_slurm/setup_local.sh --longleaf   # Longleaf HPC mode
#   ./local_slurm/setup_local.sh --restore    # Restore original configs
#   ./local_slurm/setup_local.sh --status     # Show current config status

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Config files to manage
CONFIG_FILES=(
    "targets_setup.R"
    "targets_bo_configs.R"
    "targets_run_configs.R"
    "targets_val_configs.R"
    "targets_figure_configs.R"
)

# Backup directory
BACKUP_DIR="$PROJECT_DIR/.config_backup"

# Mode file to track current mode
MODE_FILE="$SCRIPT_DIR/.current_mode"

setup_local() {
    local mode="${1:-quick}"
    local mode_dir="$SCRIPT_DIR/$mode"

    # Validate mode
    if [[ ! -d "$mode_dir" ]]; then
        echo "Error: Mode '$mode' not found. Expected directory: $mode_dir"
        echo "Available modes: quick, full, longleaf"
        exit 1
    fi

    if [[ "$mode" == "longleaf" ]]; then
        echo "Setting up LONGLEAF HPC configuration..."
    else
        echo "Setting up LOCAL configuration for 20-CPU desktop..."
    fi
    echo "Mode: $mode"
    echo ""

    # Create backup directory
    mkdir -p "$BACKUP_DIR"

    # Create logs directory if it doesn't exist
    mkdir -p "$PROJECT_DIR/logs"

    for config in "${CONFIG_FILES[@]}"; do
        local_config="$mode_dir/$config"
        main_config="$PROJECT_DIR/$config"
        backup_config="$BACKUP_DIR/$config"

        # Check if local config exists
        if [[ ! -f "$local_config" ]]; then
            echo "  Warning: $local_config not found, skipping"
            continue
        fi

        # Backup original if it exists and hasn't been backed up
        if [[ -f "$main_config" && ! -L "$main_config" && ! -f "$backup_config" ]]; then
            echo "  Backing up: $config"
            cp "$main_config" "$backup_config"
        fi

        # Remove existing file/symlink
        if [[ -e "$main_config" || -L "$main_config" ]]; then
            rm "$main_config"
        fi

        # Create symlink to local config
        ln -s "$local_config" "$main_config"
        echo "  Linked: $config -> local_slurm/$mode/$config"
    done

    # Save current mode
    echo "$mode" > "$MODE_FILE"

    echo ""
    echo "=== Configuration Summary ==="
    if [[ "$mode" == "quick" ]]; then
        echo "Mode:          QUICK (fast testing)"
        echo "ninit:         4   (vs 30 in full)"
        echo "bo_n_iter:     4   (vs 50 in full)"
        echo "ninit_full:    19  (vs 100 in full)"
        echo "CPUs per task: 19  (local desktop limit)"
    elif [[ "$mode" == "longleaf" ]]; then
        echo "Mode:          LONGLEAF (UNC HPC)"
        echo "ninit:         50  (50 CPUs per BO worker)"
        echo "bo_n_iter:     100 (tcgacptac) / 50 (bladder)"
        echo "ninit_full:    100"
        echo "ncores_grid:   50  (matches CPUs per BO worker)"
        echo "CPU limit:     200 (via DESURV_CPU_LIMIT env var)"
    else
        echo "Mode:          FULL (HPC-quality)"
        echo "ninit:         30"
        echo "bo_n_iter:     50"
        echo "ninit_full:    100"
        echo "CPUs per task: 19  (local desktop limit)"
    fi
    echo ""
    if [[ "$mode" == "longleaf" ]]; then
        echo "To run the pipeline on Longleaf:"
        echo "  sbatch local_slurm/longleaf/_targets.sh"
        echo ""
        echo "Or for simulations:"
        echo "  sbatch local_slurm/longleaf/_targets_sims.sh"
    else
        echo "To run the pipeline:"
        echo "  sbatch local_slurm/_targets.sh"
        echo ""
        echo "Or for simulations:"
        echo "  sbatch local_slurm/_targets_sims.sh"
        echo ""
        echo "Or for bladder cancer:"
        echo "  sbatch local_slurm/_targets_bladder.sh"
    fi
    echo ""
    echo "To switch modes:"
    echo "  ./local_slurm/setup_local.sh --quick      # Fast testing"
    echo "  ./local_slurm/setup_local.sh --full       # Full quality"
    echo "  ./local_slurm/setup_local.sh --longleaf   # Longleaf HPC"
    echo ""
    echo "To restore original configs:"
    echo "  ./local_slurm/setup_local.sh --restore"
}

restore_original() {
    echo "Restoring original configuration..."

    for config in "${CONFIG_FILES[@]}"; do
        main_config="$PROJECT_DIR/$config"
        backup_config="$BACKUP_DIR/$config"

        # Remove symlink if it exists
        if [[ -L "$main_config" ]]; then
            rm "$main_config"
        fi

        # Restore from backup if available
        if [[ -f "$backup_config" ]]; then
            cp "$backup_config" "$main_config"
            echo "  Restored: $config"
        else
            echo "  No backup for: $config (may not have existed originally)"
        fi
    done

    # Remove mode file
    rm -f "$MODE_FILE"

    echo ""
    echo "Original configuration restored."
}

check_status() {
    echo "Configuration status:"
    echo ""

    # Check current mode
    if [[ -f "$MODE_FILE" ]]; then
        current_mode=$(cat "$MODE_FILE")
        echo "Current mode: $current_mode"
    else
        echo "Current mode: (not set - using original configs)"
    fi
    echo ""

    for config in "${CONFIG_FILES[@]}"; do
        main_config="$PROJECT_DIR/$config"

        if [[ -L "$main_config" ]]; then
            target=$(readlink "$main_config")
            # Extract mode from path
            if [[ "$target" == *"/quick/"* ]]; then
                echo "  $config -> QUICK"
            elif [[ "$target" == *"/full/"* ]]; then
                echo "  $config -> FULL"
            elif [[ "$target" == *"/longleaf/"* ]]; then
                echo "  $config -> LONGLEAF"
            else
                echo "  $config -> $target"
            fi
        elif [[ -f "$main_config" ]]; then
            echo "  $config (ORIGINAL)"
        else
            echo "  $config (MISSING)"
        fi
    done

    echo ""
    echo "Mode comparison:"
    echo "  Parameter       QUICK    FULL     LONGLEAF"
    echo "  ----------      -----    ----     --------"
    echo "  ninit             4       30         50"
    echo "  bo_n_init         4       20         50"
    echo "  bo_n_iter         4       50        100"
    echo "  bo_candidate    200     4000       4000"
    echo "  ninit_full       19      100        100"
    echo "  ncores_grid       5        5         50"
    echo "  CPUs/task        19       19        200"
}

show_help() {
    echo "Usage: $0 [--quick|--full|--longleaf|--restore|--status|--help]"
    echo ""
    echo "Modes:"
    echo "  --quick     Quick testing mode (default)"
    echo "              - ninit=4, bo_n_iter=4, ninit_full=19"
    echo "              - Good for testing pipeline runs quickly"
    echo ""
    echo "  --full      Full quality mode"
    echo "              - ninit=30, bo_n_iter=50, ninit_full=100"
    echo "              - Matches HPC quality, but takes much longer"
    echo "              - CPUs still capped at 19 for local desktop"
    echo ""
    echo "  --longleaf  UNC Longleaf HPC mode"
    echo "              - ninit=50, bo_n_iter=100, ninit_full=100"
    echo "              - ncores_grid=50 (50 CPUs per BO worker)"
    echo "              - 202 CV workers, 120 med_mem workers"
    echo "              - Submit via: sbatch local_slurm/longleaf/_targets.sh"
    echo ""
    echo "Options:"
    echo "  --restore   Restore original HPC configuration"
    echo "  --status    Show current configuration status"
    echo "  --help      Show this help"
    echo ""
    echo "Examples:"
    echo "  $0                  # Set up quick mode (default)"
    echo "  $0 --full           # Set up full quality mode"
    echo "  $0 --longleaf       # Set up Longleaf HPC mode"
    echo "  $0 --status         # Check current configuration"
    echo "  $0 --restore        # Restore student's original configs"
}

# Main
case "${1:-}" in
    --restore|-r)
        restore_original
        ;;
    --status|-s)
        check_status
        ;;
    --help|-h)
        show_help
        ;;
    --full|-f)
        setup_local "full"
        ;;
    --longleaf|-l)
        setup_local "longleaf"
        ;;
    --quick|-q|"")
        setup_local "quick"
        ;;
    *)
        echo "Unknown option: $1"
        show_help
        exit 1
        ;;
esac
