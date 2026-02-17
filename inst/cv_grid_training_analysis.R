##############################################################################
# cv_grid_training_analysis.R
#
# Training-only analysis of 9 cv_grid DeSurv fits (K=2,3,5 x alpha=0,0.25/0.35,0.55).
# All fits use fixed lambda=0.3, nu=0.05, ntop=NULL, 100 initializations.
#
# Run from the DeSurv-paper repo root:
#   Rscript inst/cv_grid_training_analysis.R
#
# Produces: per-fit H-score correlations with original K=3, reference
# signature overlaps, beta structure, and per-factor Cox models (unadjusted).
#
# See: docs/plans/2026-02-17-k-sensitivity-synthesis.md for interpretation.
##############################################################################

suppressMessages({
  library(survival)
  library(data.table)
})

STORE <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"

# Load reference gene signatures (Elyada, SCISSORS, DECODER, etc.)
load("data/derv/cmbSubtypes_formatted.RData")

# Load original production fit (BO-selected, alpha=0.334, ntop=270)
orig_bundle <- readRDS(file.path(STORE, "objects/val_run_bundle_tcgacptac"))
orig_fit <- orig_bundle$fit_desurv
orig_W <- orig_fit$W
orig_beta <- orig_fit$beta

# Load training data (TCGA + CPTAC, 1970 genes x 273 samples)
bo_bundle <- readRDS(file.path(STORE, "objects/bo_bundle_selected_tcgacptac"))
X_train <- bo_bundle$data_filtered$ex
surv_info <- bo_bundle$data_filtered$sampInfo

source("R/get_top_genes.R")
orig_tops <- get_top_genes(orig_W, 270)

# Build reference signature lists
elyada_icaf <- top_genes$Elyada_CAF$iCAF
elyada_mycaf <- top_genes$Elyada_CAF$myCAF
scissors_icaf <- top_genes$SCISSORS_CAF_top25$iCAF
scissors_mycaf <- top_genes$SCISSORS_CAF_top25$myCAF
decaf_procaf <- c("IGFL2", "NOX4", "VSNL1", "BICD1", "NPR3", "ETV1", "ITGA11", "CNIH3", "COL11A1")
decaf_restcaf <- c("CHRDL1", "OGN", "PI16", "ANK2", "ABCA8", "TGFBR3", "FBLN5", "SCARA5", "KIAA1217")
decoder_immune <- top_genes$DECODER$Immune
decoder_actstroma <- top_genes$DECODER$ActivatedStroma
decoder_normstroma <- top_genes$DECODER$NormalStroma
decoder_basaltumor <- top_genes$DECODER$BasalTumor
decoder_classtumor <- top_genes$DECODER$ClassicalTumor

cat("Gene universe:", nrow(X_train), "genes x", ncol(X_train), "samples\n\n")

# Load all 9 cv_grid fits
fits <- list(
  K2_a0   = readRDS(file.path(STORE, "objects/cv_grid_fit_083dae3cd3bba1c5")),
  K3_a0   = readRDS(file.path(STORE, "objects/cv_grid_fit_22ac9cbd8337ebdf")),
  K5_a0   = readRDS(file.path(STORE, "objects/cv_grid_fit_128a06d24f80875e")),
  K2_a35  = readRDS(file.path(STORE, "objects/cv_grid_fit_e675ee353d14ef89")),
  K2_a55  = readRDS(file.path(STORE, "objects/cv_grid_fit_7ced5aaef3b594d7")),
  K3_a35  = readRDS(file.path(STORE, "objects/cv_grid_fit_1555a2f1fcc58ab2")),
  K3_a55  = readRDS(file.path(STORE, "objects/cv_grid_fit_ce943309f0b3e212")),
  K5_a25  = readRDS(file.path(STORE, "objects/cv_grid_fit_5c725504aa99734e")),
  K5_a55  = readRDS(file.path(STORE, "objects/cv_grid_fit_ed2d68b41575585f"))
)

ref_sigs <- list(
  "Elyada_iCAF" = elyada_icaf,
  "Elyada_myCAF" = elyada_mycaf,
  "SCISSORS_iCAF" = scissors_icaf,
  "SCISSORS_myCAF" = scissors_mycaf,
  "DeCAF_restCAF" = decaf_restcaf,
  "DeCAF_proCAF" = decaf_procaf,
  "DECODER_Immune" = decoder_immune,
  "DECODER_ActStroma" = decoder_actstroma,
  "DECODER_NormStroma" = decoder_normstroma,
  "DECODER_BasalTumor" = decoder_basaltumor,
  "DECODER_ClassTumor" = decoder_classtumor
)

cat("All 9 fits loaded.\n\n")

# ============================================================
# Per-fit analysis
# ============================================================
for (nm in names(fits)) {
  fit <- fits[[nm]]
  W <- fit$fit$W
  beta <- fit$fit$beta
  k <- fit$k
  alpha <- fit$alpha

  cat("========================================\n")
  cat(sprintf("%s (K=%d, alpha=%.2f)\n", nm, k, alpha))
  cat("========================================\n")

  # H-score correlations with original K=3
  # H = X^T W, computed on the training data
  genes <- intersect(rownames(W), rownames(orig_W))
  genes <- intersect(genes, rownames(X_train))
  H_orig <- t(X_train[genes, ]) %*% orig_W[genes, ]
  H_new <- t(X_train[genes, ]) %*% W[genes, ]

  cor_mat <- cor(H_new, H_orig, method = "spearman")
  cat("\nH-score correlation with original K=3:\n")
  cat(sprintf("%-6s %8s %8s %8s\n", "", "OrigF1", "OrigF2", "OrigF3"))
  for (j in seq_len(nrow(cor_mat))) {
    cat(sprintf("F%-5d %8.3f %8.3f %8.3f\n", j, cor_mat[j, 1], cor_mat[j, 2], cor_mat[j, 3]))
  }

  # Top-270 gene overlap with orig F1 (the iCAF factor)
  orig_f1_genes <- orig_tops$top_genes[, 1]
  tops_new <- get_top_genes(W, 270)
  cat("\nTop-270 gene overlap with original F1 (iCAF):\n")
  for (j in seq_len(k)) {
    new_genes_j <- tops_new$top_genes[, j]
    ov <- length(intersect(new_genes_j, orig_f1_genes))
    cat(sprintf("  F%d: %d/270\n", j, ov))
  }

  # Identify best iCAF factor (highest H-cor with original F1)
  best_f <- which.max(cor_mat[, 1])
  cat(sprintf("\nBest iCAF match: F%d (rho=%.3f)\n", best_f, cor_mat[best_f, 1]))

  # Reference signature overlaps for best iCAF factor's top-270 genes
  icaf_genes <- tops_new$top_genes[, best_f]
  cat(sprintf("\nRef signature overlap for F%d (top-270):\n", best_f))
  for (sig_name in names(ref_sigs)) {
    ref <- ref_sigs[[sig_name]]
    ref_u <- intersect(ref, rownames(W))  # genes present in our universe
    ov <- intersect(icaf_genes, ref_u)
    pct <- if (length(ref_u) > 0) round(100 * length(ov) / length(ref_u), 1) else 0
    cat(sprintf("  %-22s %3d/%3d (%5.1f%%)\n", sig_name, length(ov), length(ref_u), pct))
  }

  # Reference overlaps for ALL factors (Elyada iCAF, DECODER Immune, DECODER Basal)
  cat("\nRef overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor):\n")
  for (j in seq_len(k)) {
    fg <- tops_new$top_genes[, j]
    e_ov <- length(intersect(fg, intersect(elyada_icaf, rownames(W))))
    d_im <- length(intersect(fg, intersect(decoder_immune, rownames(W))))
    d_bt <- length(intersect(fg, intersect(decoder_basaltumor, rownames(W))))
    cat(sprintf("  F%d: Elyada_iCAF=%d/%d  DECODER_Immune=%d/%d  DECODER_Basal=%d/%d\n",
      j, e_ov, length(intersect(elyada_icaf, rownames(W))),
      d_im, length(intersect(decoder_immune, rownames(W))),
      d_bt, length(intersect(decoder_basaltumor, rownames(W)))))
  }

  # Beta structure
  cat(sprintf("\nBeta: %s\n", paste(sprintf("F%d=%.4e", seq_along(beta), beta), collapse=", ")))

  # Per-factor Cox (unadjusted, standardized H-scores)
  H_all <- t(X_train[genes, ]) %*% W[genes, ]
  surv_df <- data.frame(time = surv_info$time, event = surv_info$event, scale(H_all))
  cat("\nPer-factor Cox (unadjusted, standardized H-scores):\n")
  cat(sprintf("%-6s %8s %8s %10s %4s\n", "Factor", "HR", "z", "p", ""))
  for (j in seq_len(k)) {
    fv <- colnames(surv_df)[j + 2]
    cx <- coxph(as.formula(paste0("Surv(time, event) ~ ", fv)), data = surv_df)
    s <- summary(cx)
    hr <- s$conf.int[1, 1]
    z_val <- s$coefficients[1, "z"]
    p_val <- s$coefficients[1, "Pr(>|z|)"]
    sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "")))
    cat(sprintf("F%-5d %8.3f %8.3f %10.4f %4s\n", j, hr, z_val, p_val, sig))
  }
  cat("\n")
}

# ============================================================
# Summary table
# ============================================================
cat("\n========== SUMMARY ==========\n")
cat(sprintf("%-10s %3s %5s %6s %8s %8s\n", "Fit", "K", "Alpha", "iCAF_F", "H_cor", "Overlap"))
for (nm in names(fits)) {
  W <- fits[[nm]]$fit$W
  genes <- intersect(rownames(W), rownames(orig_W))
  genes <- intersect(genes, rownames(X_train))
  H_orig <- t(X_train[genes, ]) %*% orig_W[genes, ]
  H_new <- t(X_train[genes, ]) %*% W[genes, ]
  cors <- cor(H_orig[, 1], H_new, method = "spearman")
  best_f <- which.max(cors)

  tops_n <- get_top_genes(W, 270)
  ov <- length(intersect(tops_n$top_genes[, best_f], orig_tops$top_genes[, 1]))

  cat(sprintf("%-10s %3d %5.2f %6d %8.3f %5d/270\n",
    nm, fits[[nm]]$k, fits[[nm]]$alpha, best_f, max(cors), ov))
}
