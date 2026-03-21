#!/bin/bash
# assemble_clean_repo.sh — Build the clean rashidlab/DeSurv-paper repo
#
# Creates a directory with only the files needed for reproduction.
# Run from the current DeSurv-paper repo root.
#
# Usage: bash assemble_clean_repo.sh [TARGET_DIR]
#   Default target: ../DeSurv-paper-clean

set -euo pipefail

SRC="$(cd "$(dirname "$0")" && pwd)"
TARGET="${1:-$(dirname "$SRC")/DeSurv-paper-clean}"

if [ -d "$TARGET" ]; then
  echo "ERROR: Target directory already exists: $TARGET"
  echo "Remove it first or specify a different path."
  exit 1
fi

echo "Assembling clean repo from: $SRC"
echo "Target directory: $TARGET"
echo ""

mkdir -p "$TARGET"

# ═══════════════════════════════════════════════════════════════════════════
# 1. Pipeline scripts
# ═══════════════════════════════════════════════════════════════════════════
echo "── Pipeline scripts ──"
mkdir -p "$TARGET/code"
cp "$SRC/code/00_helpers.R" "$TARGET/code/"
cp "$SRC/code/01_install.R" "$TARGET/code/"
cp "$SRC/code/02_load_data.R" "$TARGET/code/"
cp "$SRC/code/03_bayesian_optimization.R" "$TARGET/code/"
cp "$SRC/code/04_fit_models.R" "$TARGET/code/"
cp "$SRC/code/05_external_validation.R" "$TARGET/code/"
cp "$SRC/code/06_sensitivity_analysis.R" "$TARGET/code/"
cp "$SRC/code/07_simulations.R" "$TARGET/code/"
cp "$SRC/code/08_figures.R" "$TARGET/code/"
cp "$SRC/code/09_render_paper.R" "$TARGET/code/"
echo "  10 scripts"

# ═══════════════════════════════════════════════════════════════════════════
# 2. R helper functions (sourced by code/*.R)
# ═══════════════════════════════════════════════════════════════════════════
echo "── R helper functions ──"
mkdir -p "$TARGET/R"

# Directly sourced by pipeline scripts
for f in \
  load_data.R \
  load_data_internal.R \
  bo_helpers.R \
  pick_k_elbow.R \
  targets_config.R \
  fit_cox_model.R \
  predict_validation_scores.R \
  cv_grid_helpers.R \
  get_top_genes.R \
  figure_targets.R \
  cluster_alignment.R \
  compare_models.R \
  plot_survival.R \
  enrichment_map.R \
  plot_heatmap.R \
  heatmap.3.R \
; do
  cp "$SRC/R/$f" "$TARGET/R/"
done
echo "  16 R files"

# Simulation functions
mkdir -p "$TARGET/R/simulation_functions"
cp "$SRC/R/simulation_functions/"*.R "$TARGET/R/simulation_functions/"
echo "  $(ls "$TARGET/R/simulation_functions/" | wc -l) simulation function files"

# Top-level sim figs script
cp "$SRC/sim_figs.R" "$TARGET/"

# ═══════════════════════════════════════════════════════════════════════════
# 3. Top-level entry points
# ═══════════════════════════════════════════════════════════════════════════
echo "── Entry points ──"
cp "$SRC/Makefile" "$TARGET/"
cp "$SRC/run_pipeline.R" "$TARGET/"
cp "$SRC/export_precomputed.R" "$TARGET/"

# ═══════════════════════════════════════════════════════════════════════════
# 4. Pre-computed results
# ═══════════════════════════════════════════════════════════════════════════
echo "── Pre-computed results ──"
mkdir -p "$TARGET/results/precomputed"
cp "$SRC/results/precomputed/"*.rds "$TARGET/results/precomputed/"
echo "  $(ls "$TARGET/results/precomputed/"*.rds | wc -l) precomputed RDS files"

mkdir -p "$TARGET/results/cv_grid"
cp "$SRC/results/cv_grid/"*.rds "$TARGET/results/cv_grid/" 2>/dev/null || true
cp "$SRC/results/cv_grid/"*.csv "$TARGET/results/cv_grid/" 2>/dev/null || true
echo "  $(ls "$TARGET/results/cv_grid/" | wc -l) cv_grid files"

# ═══════════════════════════════════════════════════════════════════════════
# 5. Static figures referenced by paper
# ═══════════════════════════════════════════════════════════════════════════
echo "── Static figures ──"
mkdir -p "$TARGET/figures/cv_grid" "$TARGET/figures/km_dichot"

cp "$SRC/figures/model_schematic_final.pdf" "$TARGET/figures/"
cp "$SRC/figures/cutpoint_curve_logrank_tcgacptac.pdf" "$TARGET/figures/"

cp "$SRC/figures/cv_grid/cv_cindex_by_k_primary.pdf" "$TARGET/figures/cv_grid/"
cp "$SRC/figures/cv_grid/k3_k7_factor_correlation_heatmap.pdf" "$TARGET/figures/cv_grid/"

for f in \
  km_val_pooled_logrank_tcgacptac.pdf \
  km_val_Dijk_logrank_tcgacptac.pdf \
  km_val_Moffitt_GEO_array_logrank_tcgacptac.pdf \
  km_val_PACA_AU_logrank_tcgacptac.pdf \
  km_val_Puleo_array_logrank_tcgacptac.pdf \
  subtype_overlap_pooled_logrank_tcgacptac.pdf \
; do
  cp "$SRC/figures/km_dichot/$f" "$TARGET/figures/km_dichot/"
done
echo "  10 static PDF figures"

# ═══════════════════════════════════════════════════════════════════════════
# 6. Paper source files
# ═══════════════════════════════════════════════════════════════════════════
echo "── Paper source ──"
mkdir -p "$TARGET/paper"

# Main paper + children
cp "$SRC/paper/paper.Rmd" "$TARGET/paper/"
cp "$SRC/paper/02_introduction_REVISED.Rmd" "$TARGET/paper/"
cp "$SRC/paper/03_methods_REVISED.Rmd" "$TARGET/paper/"
cp "$SRC/paper/04_results_REVISED.Rmd" "$TARGET/paper/"
cp "$SRC/paper/05_discussion_REVISED.Rmd" "$TARGET/paper/"

# SI Appendix
cp "$SRC/paper/si_appendix.Rmd" "$TARGET/paper/"

# Precomputed loader
cp "$SRC/paper/load_precomputed.R" "$TARGET/paper/"

# Bibliography and styles
cp "$SRC/paper/references_30102025.bib" "$TARGET/paper/"
cp "$SRC/paper/pnas.csl" "$TARGET/paper/"
cp "$SRC/paper/pnas-new.cls" "$TARGET/paper/"
cp "$SRC/paper/pnasresearcharticle.sty" "$TARGET/paper/"

# LaTeX algorithm packages
cp "$SRC/paper/algorithm.sty" "$TARGET/paper/" 2>/dev/null || true
cp "$SRC/paper/algpseudocode.sty" "$TARGET/paper/" 2>/dev/null || true
cp "$SRC/paper/algorithmicx.sty" "$TARGET/paper/" 2>/dev/null || true

# Gene lists CSV (referenced in paper)
cp "$SRC/paper/gene_lists_top270_k3.csv" "$TARGET/paper/" 2>/dev/null || true

# Top-level gene list
cp "$SRC/top_genes_desurv_k3_tcgacptac.csv" "$TARGET/" 2>/dev/null || true

echo "  $(ls "$TARGET/paper/" | wc -l) paper files"

# ═══════════════════════════════════════════════════════════════════════════
# 7. Data directory (empty with README — data downloaded separately)
# ═══════════════════════════════════════════════════════════════════════════
echo "── Data placeholder ──"
mkdir -p "$TARGET/data/original"
cat > "$TARGET/data/README.md" << 'DATA_EOF'
# Data

This directory should contain the input datasets for the DeSurv analysis.

## Required files

The following datasets are needed in `data/original/`:

### Training cohorts
- `TCGA_PAAD.rds`, `TCGA_PAAD.survival_data.rds`, `TCGA_PAAD_subtype.csv`
- `CPTAC.rds`, `CPTAC.survival_data.rds`, `CPTAC_subtype.csv`

### Validation cohorts
- `Dijk.rds`, `Dijk.survival_data.rds`, `Dijk_subtype.csv`
- `Moffitt_GEO_array.rds`, `Moffitt_GEO_array.survival_data.rds`, `Moffitt_GEO_array_subtype.csv`
- `PACA_AU_array.rds`, `PACA_AU_array.survival_data.rds`, `PACA_AU_array_subtype.csv`
- `PACA_AU_seq.rds`, `PACA_AU_seq.survival_data.rds`, `PACA_AU_seq_subtype.csv`
- `Puleo_array.rds`, `Puleo_array.survival_data.rds`, `Puleo_array_subtype.csv`

### Reference data
- `cmbSubtypes.RData` (combined molecular subtypes for sensitivity analysis)

## Data sources

- **TCGA_PAAD:** The Cancer Genome Atlas (https://www.cancer.gov/tcga)
- **CPTAC:** Clinical Proteomic Tumor Analysis Consortium
- **Moffitt:** GEO accession GSE71729
- **Puleo:** ArrayExpress E-MTAB-6134
- **Dijk:** ArrayExpress E-MTAB-6830
- **PACA-AU:** ICGC data portal, EGA study EGAS00001000154

## Download

Data files are available from [Zenodo DOI: TBD].
DATA_EOF
echo "  Data README created (data files not copied — will be on Zenodo)"

# ═══════════════════════════════════════════════════════════════════════════
# 8. Slurm scripts
# ═══════════════════════════════════════════════════════════════════════════
echo "── Slurm scripts ──"
mkdir -p "$TARGET/slurm"
cat > "$TARGET/slurm/run_full_pipeline.sh" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=desurv-pipeline
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/desurv_%j.log

module load r/4.3.1
cd $SLURM_SUBMIT_DIR

make all NCORES=8
SLURM_EOF
echo "  1 slurm script"

# ═══════════════════════════════════════════════════════════════════════════
# 9. Top-level files
# ═══════════════════════════════════════════════════════════════════════════
echo "── Top-level files ──"

cat > "$TARGET/.gitignore" << 'GI_EOF'
# Data (hosted on Zenodo)
data/original/*.rds
data/original/*.RData
data/original/*.csv

# R artifacts
.Rhistory
.RData
.Rproj.user/
Rplots.pdf

# renv library (generated by renv::restore())
renv/library/
renv/staging/
renv/cellar/

# Build artifacts
paper/*.log
paper/*.tex
paper/*.pdf
paper/*_files/
logs/

# OS
.DS_Store
*~
GI_EOF

cat > "$TARGET/LICENSE" << 'LIC_EOF'
MIT License

Copyright (c) 2026 Amber M. Young, Naim U. Rashid

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
LIC_EOF

echo "  .gitignore, LICENSE created"

# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "Clean repo assembled at: $TARGET"
echo ""
echo "Directory structure:"
find "$TARGET" -type f | sed "s|$TARGET/||" | sort | head -80
echo ""
TOTAL_FILES=$(find "$TARGET" -type f | wc -l)
TOTAL_SIZE=$(du -sh "$TARGET" | cut -f1)
echo "Total: $TOTAL_FILES files, $TOTAL_SIZE"
echo ""
echo "Next steps:"
echo "  1. Copy data files to $TARGET/data/original/ (or symlink)"
echo "  2. cd $TARGET && git init && git add -A && git commit -m 'initial clean repo'"
echo "  3. Write README.md"
echo "  4. Set up renv: Rscript -e 'renv::init()'"
echo "  5. Test: make paper (from precomputed)"
echo "  6. Test: make quick (with data symlinked)"
echo "════════════════════════════════════════════════════════════════"
