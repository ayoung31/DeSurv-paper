##############################################################################
# cv_grid_training_analysis.R
#
# Training-only analysis of cv_grid DeSurv fits over
# K=2,3,5,7,9 x alpha=0,0.25,0.35,0.55,0.75,0.85 (30 combos).
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

# Load all cv_grid_fit branches from store, then filter by (k, alpha)
cat("Loading all cv_grid_fit branches from store...\n")
fit_files <- list.files(
  file.path(STORE, "objects"),
  pattern = "^cv_grid_fit_[0-9a-f]+$",
  full.names = TRUE
)
all_fits <- lapply(fit_files, readRDS)
cat(sprintf("  Loaded %d branches.\n", length(all_fits)))

find_fit <- function(all_fits, k, alpha, ntop = NA) {
  for (fe in all_fits) {
    if (is.null(fe) || is.null(fe$k)) next
    fe_ntop <- if (is.null(fe$ntop)) NA_real_ else as.numeric(fe$ntop)
    ntop_match <- (is.na(ntop) && is.na(fe_ntop)) ||
                  (!is.na(ntop) && !is.na(fe_ntop) && fe_ntop == ntop)
    if (fe$k == k && abs(fe$alpha - alpha) < 1e-6 && ntop_match) {
      return(fe)
    }
  }
  stop(sprintf("No fit found for k=%d, alpha=%.2f, ntop=%s",
               k, alpha, ifelse(is.na(ntop), "NULL", ntop)))
}

# Define desired parameter combos (ntop=NA means ntop=NULL / all genes)
GRID_K     <- c(2, 3, 5, 7, 9)
GRID_ALPHA <- c(0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95)

fit_specs <- list()
for (ki in GRID_K) {
  for (ai in GRID_ALPHA) {
    nm <- sprintf("K%d_a%s", ki, sub("\\.", "", sprintf("%.2f", ai)))
    fit_specs[[nm]] <- list(k = ki, alpha = ai)
  }
}

fits <- lapply(fit_specs, function(spec) {
  find_fit(all_fits, k = spec$k, alpha = spec$alpha)
})

# Verify loaded params match expectations
for (nm in names(fits)) {
  exp <- fit_specs[[nm]]
  got <- fits[[nm]]
  if (got$k != exp$k || abs(got$alpha - exp$alpha) > 1e-6) {
    stop(sprintf("Mismatch for %s: expected k=%d alpha=%.2f, got k=%d alpha=%.2f",
                 nm, exp$k, exp$alpha, got$k, got$alpha))
  }
  cat(sprintf("  %s: k=%d alpha=%.2f OK\n", nm, got$k, got$alpha))
}

# Free memory from unused fits
rm(all_fits)
gc(verbose = FALSE)

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

cat(sprintf("All %d fits loaded.\n\n", length(fits)))

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
# Summary table (long format)
# ============================================================
cat("\n========== SUMMARY ==========\n")
cat(sprintf("%-14s %3s %5s %6s %8s %8s\n", "Fit", "K", "Alpha", "iCAF_F", "H_cor", "Overlap"))

# Collect H_cor values for the k x alpha matrix
hcor_matrix <- matrix(NA_real_, nrow = length(GRID_K), ncol = length(GRID_ALPHA),
                       dimnames = list(paste0("K=", GRID_K),
                                       paste0("a=", sprintf("%.2f", GRID_ALPHA))))

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

  hcor_val <- max(cors)
  cat(sprintf("%-14s %3d %5.2f %6d %8.3f %5d/270\n",
    nm, fits[[nm]]$k, fits[[nm]]$alpha, best_f, hcor_val, ov))

  # Store in matrix (round to avoid floating-point mismatch from seq())
  row_idx <- match(fits[[nm]]$k, GRID_K)
  col_idx <- match(round(fits[[nm]]$alpha, 2), round(GRID_ALPHA, 2))
  hcor_matrix[row_idx, col_idx] <- hcor_val
}

# ============================================================
# H_cor matrix: k (rows) x alpha (columns)
# ============================================================
cat("\n========== H_cor MATRIX (best iCAF factor rho with original F1) ==========\n")
cat(sprintf("%-6s", "K"))
for (ai in GRID_ALPHA) cat(sprintf(" %8s", sprintf("a=%.2f", ai)))
cat("\n")
for (i in seq_along(GRID_K)) {
  cat(sprintf("%-6d", GRID_K[i]))
  for (j in seq_along(GRID_ALPHA)) {
    val <- hcor_matrix[i, j]
    if (is.na(val)) {
      cat(sprintf(" %8s", "NA"))
    } else {
      cat(sprintf(" %8.3f", val))
    }
  }
  cat("\n")
}
