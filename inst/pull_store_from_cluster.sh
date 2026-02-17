#!/bin/bash
# Pull paper's targets store objects from Longleaf cluster
# Uses SSH control socket for single authentication
# Usage: bash inst/pull_store_from_cluster.sh
#
# Pulls core K=3 model targets, K=5 (elbow K) fit,
# BO diagnostics, and CV grid search summaries (k × alpha)
# for K-sensitivity comparison analyses.

set -euo pipefail

REMOTE="nur2@longleaf.unc.edu"
REMOTE_DIR="/work/users/a/y/ayoung31/DeSurv-paper"
LOCAL_DIR="$HOME/Downloads/DeSurv-paper"
STORE="store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"
SOCK="/tmp/longleaf-sock"
SSH_OPTS="-o ControlPath=${SOCK}"

echo "Pulling targets from ${REMOTE}:${REMOTE_DIR}/${STORE}"
echo "To: ${LOCAL_DIR}/${STORE}"
echo ""

# Open persistent SSH connection (single password prompt)
echo "--- Opening SSH control socket (authenticate once) ---"
ssh -fNM -o ControlPath="${SOCK}" "${REMOTE}"
trap 'echo "Closing SSH socket..."; ssh -O exit -o ControlPath="${SOCK}" "${REMOTE}" 2>/dev/null' EXIT

# Helper: pull a target, warn on failure instead of aborting
pull_target() {
    local target="$1"
    local label="${2:-}"
    if scp ${SSH_OPTS} "${REMOTE}:${REMOTE_DIR}/${STORE}/objects/${target}" \
        "${LOCAL_DIR}/${STORE}/objects/" 2>/dev/null; then
        echo "  OK  ${target}"
    else
        echo "  SKIP ${target} (not found on remote)"
    fi
}

# Create local directories
mkdir -p "${LOCAL_DIR}/${STORE}/objects"
mkdir -p "${LOCAL_DIR}/${STORE}/meta"

# Metadata
echo "--- Pulling metadata ---"
scp ${SSH_OPTS} "${REMOTE}:${REMOTE_DIR}/${STORE}/meta/meta" \
    "${LOCAL_DIR}/${STORE}/meta/"

# ── Core K=3 model targets ──────────────────────────────────
echo ""
echo "--- Core K=3 model targets ---"
for t in \
    val_run_bundle_tcgacptac \
    tar_data_tcgacptac \
    bo_bundle_selected_tcgacptac \
    data_val_filtered_tcgacptac \
    data_val_filtered_surv_tcgacptac \
    desurv_lp_stats_tcgacptac \
    desurv_optimal_z_cutpoint_tcgacptac \
; do
    pull_target "$t"
done

# ── BO diagnostics (GP surface, K selection) ────────────────
echo ""
echo "--- BO diagnostics ---"
for t in \
    tar_k_selection_tcgacptac \
    tar_params_best_tcgacptac \
    desurv_bo_history_tcgacptac \
; do
    pull_target "$t"
done

# ── K=5 / elbow-K DeSurv fit and validation ─────────────────
echo ""
echo "--- K=5 (elbow K) DeSurv targets ---"
for t in \
    std_nmf_selected_k_tcgacptac \
    tar_params_best_elbowk_tcgacptac \
    tar_fit_desurv_elbowk_tcgacptac \
    tar_tops_desurv_elbowk_tcgacptac \
    tar_data_filtered_elbowk_tcgacptac \
    data_val_filtered_elbowk_tcgacptac \
    val_latent_desurv_elbowk_tcgacptac \
    val_cindex_desurv_elbowk_tcgacptac \
    val_predictions_desurv_elbowk_tcgacptac \
; do
    pull_target "$t"
done

# ── Standard NMF at elbow K (for comparison) ────────────────
echo ""
echo "--- Standard NMF at elbow K ---"
for t in \
    fit_std_elbowk_tcgacptac \
    tar_tops_std_elbowk_tcgacptac \
    val_latent_std_elbowk_tcgacptac \
    val_cindex_std_elbowk_tcgacptac \
; do
    pull_target "$t"
done

# ── CV grid search summaries (k × alpha exhaustive grid) ──
echo ""
echo "--- CV grid search summaries ---"
for t in \
    cv_grid_summary \
    cv_grid_best_alpha \
    cv_grid_val_summary \
    cv_grid_fit_summary \
    cv_grid_optimal_cutpoint \
    cv_grid_alpha0 \
    cv_grid_data \
; do
    pull_target "$t"
done

# ── CV grid full fits: K=2,3,5 at key alphas (factor structure analysis) ──
# These 9 branches enable apples-to-apples factor structure comparison
# across K values with identical lambda=0.3, nu=0.05, 100 inits.
# Mapping from (K, alpha, ntop=NULL) to branch hash determined via
# expand.grid order matching cv_grid_fit_summary rows to metadata children.
echo ""
echo "--- CV grid fits: K=2,3,5 at key alphas ---"
for t in \
    cv_grid_fit_083dae3cd3bba1c5 \
    cv_grid_fit_22ac9cbd8337ebdf \
    cv_grid_fit_128a06d24f80875e \
    cv_grid_fit_e675ee353d14ef89 \
    cv_grid_fit_7ced5aaef3b594d7 \
    cv_grid_fit_1555a2f1fcc58ab2 \
    cv_grid_fit_ce943309f0b3e212 \
    cv_grid_fit_5c725504aa99734e \
    cv_grid_fit_ed2d68b41575585f \
; do
    pull_target "$t"
done
# Hash mapping:
#   K=2 alpha=0.00  cv_grid_fit_083dae3cd3bba1c5  (stdNMF baseline)
#   K=3 alpha=0.00  cv_grid_fit_22ac9cbd8337ebdf  (stdNMF baseline)
#   K=5 alpha=0.00  cv_grid_fit_128a06d24f80875e  (stdNMF baseline)
#   K=2 alpha=0.35  cv_grid_fit_e675ee353d14ef89  (val-optimal)
#   K=2 alpha=0.55  cv_grid_fit_7ced5aaef3b594d7  (matched to K=3 CV-optimal)
#   K=3 alpha=0.35  cv_grid_fit_1555a2f1fcc58ab2  (val-optimal)
#   K=3 alpha=0.55  cv_grid_fit_ce943309f0b3e212  (CV-optimal)
#   K=5 alpha=0.25  cv_grid_fit_5c725504aa99734e  (val-optimal)
#   K=5 alpha=0.55  cv_grid_fit_ed2d68b41575585f  (matched to K=3 CV-optimal)
#
# To pull ALL 440 fit branches (large!), uncomment below:
# echo ""
# echo "--- CV grid full fits (all 440 branches) ---"
# ssh ${SSH_OPTS} "${REMOTE}" "ls ${REMOTE_DIR}/${STORE}/objects/cv_grid_fit_*" 2>/dev/null | while read remote_path; do
#     fname=$(basename "$remote_path")
#     scp ${SSH_OPTS} "${REMOTE}:${remote_path}" "${LOCAL_DIR}/${STORE}/objects/" 2>/dev/null && echo "  OK  ${fname}" || echo "  SKIP ${fname}"
# done

echo ""
echo "Done. Verify with:"
echo "  ls -lh ${LOCAL_DIR}/${STORE}/objects/"
