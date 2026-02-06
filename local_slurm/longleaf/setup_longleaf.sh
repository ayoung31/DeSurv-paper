#!/bin/bash
# setup_longleaf.sh - One-time environment setup for DeSurv on UNC Longleaf
#
# Run once after cloning the repo on Longleaf to:
#   1. Verify R module and prerequisites
#   2. Install all required R packages
#   3. Install the DeSurv package
#   4. Verify data files exist
#   5. Create required directories
#   6. Activate Longleaf configs
#   7. Run preflight check
#
# Usage:
#   bash local_slurm/longleaf/setup_longleaf.sh
#   bash local_slurm/longleaf/setup_longleaf.sh --skip-packages
#   bash local_slurm/longleaf/setup_longleaf.sh --help

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

SKIP_PACKAGES=false

show_help() {
    cat <<'EOF'
Usage: bash local_slurm/longleaf/setup_longleaf.sh [OPTIONS]

One-time setup for running DeSurv pipelines on UNC Longleaf.

Options:
  --skip-packages   Skip R package installation (if already installed)
  --help            Show this help

Steps performed:
  1. Check prerequisites (R module, Rscript)
  2. Install R dependencies (CRAN + Bioconductor)
  3. Install DeSurv package from ../DeSurv
  4. Verify data files exist
  5. Create output directories (logs/, figures/)
  6. Activate Longleaf pipeline configs
  7. Run preflight check

After setup, submit pipelines with:
  sbatch local_slurm/longleaf/_targets.sh
  sbatch local_slurm/longleaf/_targets_sims.sh
EOF
}

pass() { echo -e "  ${GREEN}[OK]${NC} $1"; }
warn() { echo -e "  ${YELLOW}[WARN]${NC} $1"; }
fail() { echo -e "  ${RED}[FAIL]${NC} $1"; }

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --skip-packages) SKIP_PACKAGES=true ;;
        --help|-h) show_help; exit 0 ;;
        *) echo "Unknown option: $arg"; show_help; exit 1 ;;
    esac
done

echo ""
echo "========================================================"
echo "  DeSurv Longleaf Setup"
echo "========================================================"
echo ""
echo "Project: $PROJECT_DIR"
echo "Time:    $(date)"
echo ""

# =========================================================================
# 1. Check prerequisites
# =========================================================================
echo "${BOLD}1. Checking prerequisites...${NC}"
echo "   ────────────────────────"

# Check module command exists
if ! command -v module &> /dev/null; then
    fail "module command not found (are you on Longleaf?)"
    exit 1
fi
pass "module command available"

# Load R module
if module load r/4.4.0 2>/dev/null; then
    pass "module load r/4.4.0 succeeded"
else
    fail "could not load r/4.4.0 module"
    echo "       Available R modules:"
    module avail r/ 2>&1 | head -5
    exit 1
fi

# Check Rscript
if command -v Rscript &> /dev/null; then
    r_version=$(Rscript --version 2>&1 | head -1)
    pass "Rscript on PATH ($r_version)"
else
    fail "Rscript not found after module load"
    exit 1
fi

# Check git
if command -v git &> /dev/null; then
    pass "git available"
else
    warn "git not found (some targets provenance features may not work)"
fi

echo ""

# =========================================================================
# 2. Install R dependencies
# =========================================================================
echo "${BOLD}2. Installing R packages...${NC}"
echo "   ────────────────────────"

if [ "$SKIP_PACKAGES" = true ]; then
    echo "  Skipping (--skip-packages)"
else
    Rscript --vanilla -e '
# CRAN packages required by the pipeline
cran_pkgs <- c(
  # Pipeline infrastructure
  "targets", "tarchetypes", "crew", "crew.cluster", "gert",
  # Data manipulation
  "dplyr", "purrr", "tibble", "tidyverse", "tidyselect", "digest",
  # ML / statistics
  "survival", "glmnet", "caret", "pec", "cvwrapr",
  # NMF / visualization
  "NMF", "pheatmap", "ggrepel", "survminer",
  # Parallel
  "parallel", "foreach", "doParallel", "doMC",
  # Other
  "rmarkdown", "webshot2", "devtools"
)

# Bioconductor packages
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db")

# Install missing CRAN packages
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
  cat("All Bioconductor packages already installed.\n")
}

# Verify all packages load
all_pkgs <- c(cran_pkgs, bioc_pkgs)
failed <- character(0)
for (pkg in all_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
  }
}

if (length(failed) > 0) {
  cat("ERROR: Failed to install:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("All", length(all_pkgs), "packages verified.\n")
}
'

    if [ $? -eq 0 ]; then
        pass "R package installation complete"
    else
        fail "R package installation failed (see errors above)"
        exit 1
    fi
fi

echo ""

# =========================================================================
# 3. Install DeSurv package
# =========================================================================
echo "${BOLD}3. Installing DeSurv package...${NC}"
echo "   ────────────────────────────"

DESURV_DIR="$(cd "$PROJECT_DIR/.." && pwd)/DeSurv"

if [ ! -d "$DESURV_DIR" ]; then
    fail "DeSurv package not found at $DESURV_DIR"
    echo "       Clone it first: git clone <desurv-repo-url> $(dirname "$DESURV_DIR")/DeSurv"
    exit 1
fi
pass "DeSurv source found at $DESURV_DIR"

Rscript --vanilla -e "devtools::install_local('$DESURV_DIR', upgrade = 'never', force = TRUE)"

if [ $? -eq 0 ]; then
    pass "DeSurv installed successfully"
else
    fail "DeSurv installation failed"
    exit 1
fi

echo ""

# =========================================================================
# 4. Check data files
# =========================================================================
echo "${BOLD}4. Checking data files...${NC}"
echo "   ──────────────────────"

data_ok=true

# Key original data files (TCGA/CPTAC pancreatic)
for f in CPTAC.rds CPTAC.survival_data.rds; do
    if [ -f "$PROJECT_DIR/data/original/$f" ]; then
        pass "data/original/$f"
    else
        fail "data/original/$f not found"
        data_ok=false
    fi
done

# Key derived data files
for f in cmbSubtypes_formatted.RData; do
    if [ -f "$PROJECT_DIR/data/derv/$f" ]; then
        pass "data/derv/$f"
    else
        warn "data/derv/$f not found (some analyses may be unavailable)"
    fi
done

# Count total data files
n_original=$(ls "$PROJECT_DIR/data/original/"*.rds 2>/dev/null | wc -l)
n_original_rdata=$(ls "$PROJECT_DIR/data/original/"*.RData 2>/dev/null | wc -l)
echo "  Found $n_original RDS + $n_original_rdata RData files in data/original/"

if [ "$data_ok" = false ]; then
    echo ""
    warn "Missing data files. Ensure data is symlinked from /proj/rashidlab/ or copied."
fi

echo ""

# =========================================================================
# 5. Create directories
# =========================================================================
echo "${BOLD}5. Creating directories...${NC}"
echo "   ────────────────────────"

for dir in logs figures/panels figures/sim; do
    target_dir="$PROJECT_DIR/$dir"
    if [ -d "$target_dir" ]; then
        pass "$dir/ (exists)"
    else
        mkdir -p "$target_dir"
        pass "$dir/ (created)"
    fi
done

echo ""

# =========================================================================
# 6. Activate Longleaf configs
# =========================================================================
echo "${BOLD}6. Activating Longleaf configs...${NC}"
echo "   ──────────────────────────────"

"$PROJECT_DIR/local_slurm/setup_local.sh" --longleaf

echo ""

# =========================================================================
# 7. Run preflight check
# =========================================================================
echo "${BOLD}7. Running preflight check...${NC}"
echo "   ──────────────────────────"

# Preflight may fail on warnings (expected on fresh setup), so don't exit on error
set +e
"$PROJECT_DIR/scripts/preflight_check.sh"
preflight_exit=$?
set -e

echo ""

# =========================================================================
# Summary
# =========================================================================
echo "========================================================"
echo "  Setup Complete"
echo "========================================================"
echo ""
echo "Next steps:"
echo ""
echo "  1. Submit the main pipeline:"
echo "     cd $PROJECT_DIR"
echo "     sbatch local_slurm/longleaf/_targets.sh"
echo ""
echo "  2. Submit simulations (if needed):"
echo "     sbatch local_slurm/longleaf/_targets_sims.sh"
echo ""
echo "  3. Monitor progress:"
echo "     bash local_slurm/longleaf/monitor.sh"
echo "     bash local_slurm/longleaf/monitor.sh --progress"
echo "     bash local_slurm/longleaf/monitor.sh --workers"
echo ""
if [ $preflight_exit -ne 0 ]; then
    echo -e "  ${YELLOW}Note: Preflight had warnings (normal for fresh setup).${NC}"
    echo "  Review above and re-run: ./scripts/preflight_check.sh"
    echo ""
fi
