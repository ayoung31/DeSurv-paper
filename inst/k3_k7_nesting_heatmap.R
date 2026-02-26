##############################################################################
# k3_k7_nesting_heatmap.R
#
# Computes Spearman correlations between K=3 (production fit) and K=7
# (cv_grid fit, alpha=0.35, ntop=270) W matrices to assess whether K=3
# factors are nested within the 7-factor solution.
#
# Addresses GitHub issue #11, question 3.
#
# Run from the DeSurv-paper repo root:
#   Rscript inst/k3_k7_nesting_heatmap.R
##############################################################################

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

STORE         <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"
TARGET_NTOP   <- 270
TARGET_LAMBDA <- 0.349
TARGET_NU     <- 0.056
OUT_DIR       <- "figures/cv_grid"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load production K=3 W matrix ----
cat("Loading production K=3 fit...\n")
orig_bundle <- readRDS(file.path(STORE, "objects/val_run_bundle_tcgacptac"))
W_k3 <- orig_bundle$fit_desurv$W
colnames(W_k3) <- paste0("K3_F", seq_len(ncol(W_k3)))
cat(sprintf("  K=3 W: %d genes x %d factors\n", nrow(W_k3), ncol(W_k3)))

# ---- Load K=7 cv_grid fit (alpha=0.35, ntop=270) ----
cat("Loading cv_grid K=7 fit (alpha=0.35, ntop=270)...\n")
fit_files <- list.files(
  file.path(STORE, "objects"),
  pattern = "^cv_grid_fit_[0-9a-f]+$",
  full.names = TRUE
)
all_fits <- lapply(fit_files, readRDS)
cat(sprintf("  Loaded %d cv_grid branches.\n", length(all_fits)))

find_fit <- function(all_fits, k, alpha, ntop, lambda, nu) {
  for (fe in all_fits) {
    if (is.null(fe) || is.null(fe$k)) next
    fe_ntop   <- if (is.null(fe$ntop))   NA_real_ else as.numeric(fe$ntop)
    fe_lambda <- if (is.null(fe$lambda)) NA_real_ else as.numeric(fe$lambda)
    fe_nu     <- if (is.null(fe$nu))     NA_real_ else as.numeric(fe$nu)
    ntop_ok   <- !is.na(fe_ntop)   && fe_ntop   == ntop
    lam_ok    <- !is.na(fe_lambda) && abs(fe_lambda - lambda) < 1e-9
    nu_ok     <- !is.na(fe_nu)     && abs(fe_nu  - nu)     < 1e-9
    if (fe$k == k && abs(fe$alpha - alpha) < 1e-6 && ntop_ok && lam_ok && nu_ok) {
      return(fe)
    }
  }
  stop(sprintf("No fit found for k=%d, alpha=%.2f, ntop=%d, lambda=%.3f, nu=%.3f",
               k, alpha, ntop, lambda, nu))
}

fit_k7 <- find_fit(all_fits, k = 7, alpha = 0.35,
                   ntop = TARGET_NTOP, lambda = TARGET_LAMBDA, nu = TARGET_NU)
W_k7 <- fit_k7$fit$W
colnames(W_k7) <- paste0("K7_F", seq_len(ncol(W_k7)))
cat(sprintf("  K=7 W: %d genes x %d factors\n", nrow(W_k7), ncol(W_k7)))
rm(all_fits); gc(verbose = FALSE)

# ---- Align genes ----
common_genes <- intersect(rownames(W_k3), rownames(W_k7))
cat(sprintf("  Common genes: %d\n", length(common_genes)))
W_k3 <- W_k3[common_genes, , drop = FALSE]
W_k7 <- W_k7[common_genes, , drop = FALSE]

# ---- Compute Spearman correlations: K3 factors x K7 factors ----
cat("Computing Spearman correlations...\n")
cor_mat <- cor(W_k3, W_k7, method = "spearman")
# cor_mat is 3 x 7: rows = K3 factors, cols = K7 factors

# ---- Reshape to long for ggplot ----
cor_long <- as.data.frame(cor_mat) |>
  tibble::rownames_to_column("K3_factor") |>
  tidyr::pivot_longer(cols = -K3_factor,
                      names_to = "K7_factor",
                      values_to = "spearman_r") |>
  dplyr::mutate(
    K3_factor = factor(K3_factor, levels = rev(colnames(W_k3))),
    K7_factor = factor(K7_factor, levels = colnames(W_k7))
  )

# ---- Print correlation matrix ----
cat("\nSpearman correlation matrix (K=3 rows x K=7 cols):\n")
print(round(cor_mat, 3))

# ---- Find best match for each K=3 factor ----
cat("\nBest K=7 match for each K=3 factor:\n")
for (f3 in colnames(W_k3)) {
  row_cors <- cor_mat[f3, ]
  best_f7  <- names(which.max(abs(row_cors)))
  cat(sprintf("  %s -> %s (r = %.3f)\n", f3, best_f7, row_cors[best_f7]))
}

# ---- Heatmap ----
p <- ggplot(cor_long, aes(x = K7_factor, y = K3_factor, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", spearman_r)),
            size = 3.5, color = "black") +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#d6604d",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman r"
  ) +
  labs(
    title    = "Factor correspondence: K=3 (production) vs K=7 (DeSurv, \u03b1=0.35)",
    subtitle = "Spearman correlation of W-matrix gene loadings across all genes",
    x        = "K=7 factors",
    y        = "K=3 factors"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9, color = "gray40"),
    legend.position  = "right"
  )

out_path <- file.path(OUT_DIR, "k3_k7_factor_correlation_heatmap.pdf")
ggsave(out_path, p, width = 7, height = 3.5)
cat(sprintf("\nSaved heatmap to: %s\n", out_path))
