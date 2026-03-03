##############################################################################
# factor2_purity_correlation.R
#
# Checks correlation between Factor 2 score (k=5 DeSurv model) and tumor
# purity in the TCGA_PAAD + CPTAC training cohort.
#
# Purity sources:
#   A) ESTIMATE (Yoshihara et al. 2013) — computed here from unfiltered
#      expression; available for both TCGA_PAAD and CPTAC.
#   B) ABSOLUTE (Taylor et al. 2018, Nat Genet) — TCGA_PAAD only.
#      Download file before running:
#        curl -o data/original/TCGA_absolute_purity.txt \
#          "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
#
# Factor scores = W[:,2]^T %*% X  (projection onto loading vector, not H).
# W comes from tar_fit_desurv_elbowk_tcgacptac in the targets store.
# X is the pipeline-filtered training expression (tar_data_filtered_elbowk).
# ESTIMATE is run on the full unfiltered expression to retain the ~2500
# immune/stromal signature genes that may be excluded by pipeline gene
# filtering.
#
# Run from repo root:
#   Rscript inst/factor2_purity_correlation.R
##############################################################################

suppressMessages({
  library(targets)
  library(tidyestimate)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

STORE         <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"
OUT_DIR       <- "figures"
ABSOLUTE_PATH <- "data/original/TCGA_absolute_purity.txt"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

##############################################################################
# 1.  Load pipeline targets
##############################################################################
cat("Loading pipeline targets...\n")

tar_fit  <- tar_read(tar_fit_desurv_elbowk_tcgacptac,  store = STORE)
tar_data <- tar_read(tar_data_filtered_elbowk_tcgacptac, store = STORE)

sinfo  <- tar_data$sampInfo
X_filt <- tar_data$ex   # genes x samples (pipeline-filtered gene set)

cat(sprintf("  Training samples: %d  (TCGA_PAAD: %d, CPTAC: %d)\n",
    nrow(sinfo),
    sum(sinfo$dataset == "TCGA_PAAD"),
    sum(sinfo$dataset == "CPTAC")))

##############################################################################
# 2.  Compute Factor 2 scores via W[:,2]^T %*% X
##############################################################################
cat("Computing Factor 2 scores (W^T X)...\n")

w2           <- tar_fit$W[, 2, drop = FALSE]   # genes x 1
common_genes <- intersect(rownames(w2), rownames(X_filt))
cat(sprintf("  Genes in W ∩ filtered expression: %d / %d\n",
    length(common_genes), nrow(w2)))

scores_f2        <- as.numeric(t(w2[common_genes, , drop = FALSE]) %*%
                                  X_filt[common_genes, ])
names(scores_f2) <- colnames(X_filt)

sinfo$ID = rownames(sinfo)

score_df <- data.frame(
  ID      = names(scores_f2),
  factor2 = scores_f2,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(sinfo |> dplyr::select(ID, dataset),
                   by = "ID")

##############################################################################
# 3.  Load raw unfiltered expression for ESTIMATE
##############################################################################
cat("Loading raw expression for ESTIMATE...\n")
source("R/load_data_internal.R")

# load_data_internal applies dataset-specific transforms (log2 for TCGA, etc.)
# and returns the full expression matrix with sampInfo$keep flagging training
# samples.  We filter to keep == 1 to match the training cohort.
load_training_expr <- function(dataname) {
  d    <- load_data_internal(dataname)
  keep <- d$sampInfo$keep == 1
  list(
    ex       = d$ex[, d$sampInfo$ID[keep], drop = FALSE],
    sampInfo = d$sampInfo[keep, ]
  )
}

tcga  <- load_training_expr("TCGA_PAAD")
cptac <- load_training_expr("CPTAC")

cat(sprintf("  TCGA_PAAD : %d samples, %d genes\n", ncol(tcga$ex),  nrow(tcga$ex)))
cat(sprintf("  CPTAC     : %d samples, %d genes\n", ncol(cptac$ex), nrow(cptac$ex)))

##############################################################################
# 4.  Run ESTIMATE
##############################################################################
cat("Running ESTIMATE...\n")

# Helper: collapse duplicate gene symbols then run ESTIMATE
run_estimate <- function(ex_mat, dataset_name, is_affymetrix = FALSE) {
  df         <- as.data.frame(ex_mat)
  df$symbol  <- rownames(df)

  # Collapse duplicate symbols by median (mirrors pipeline convention)
  df <- dplyr::as_tibble(df) |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(dplyr::across(where(is.numeric), median), .groups = "drop") |>
    as.data.frame()
  df <- df[!is.na(df$symbol) & df$symbol != "?", ]
  rownames(df) <- df$symbol
  df$symbol    <- NULL

  est <- df |>
    filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE,
                        find_alias = TRUE) |>
    estimate_score(is_affymetrix = is_affymetrix)

  dplyr::mutate(est, dataset = dataset_name)
}

est_tcga  <- run_estimate(tcga$ex,  "TCGA_PAAD", is_affymetrix = FALSE)
est_cptac <- run_estimate(cptac$ex, "CPTAC",     is_affymetrix = FALSE)

# Yoshihara et al. 2013: tumor purity from ESTIMATE score
estimate_to_purity <- function(df) {
  dplyr::mutate(df,
    purity_estimate = pmin(pmax(cos(0.6049872 + 0.0001467884 * estimate), 0), 1)
  )
}

est_df <- dplyr::bind_rows(
  estimate_to_purity(est_tcga),
  estimate_to_purity(est_cptac)
) |>
  dplyr::rename(ID = sample) |>
  dplyr::select(ID, dataset, immune, stromal, estimate, purity_estimate)

cat(sprintf("  ESTIMATE complete: %d / %d samples\n",
    nrow(est_df), nrow(score_df)))

##############################################################################
# 5.  Load ABSOLUTE purity (TCGA only)
##############################################################################
abs_df <- NULL

if (file.exists(ABSOLUTE_PATH)) {
  cat("Loading ABSOLUTE purity...\n")

  abs_raw        <- read.delim(ABSOLUTE_PATH, stringsAsFactors = FALSE,
                                check.names = FALSE)
  names(abs_raw) <- tolower(trimws(names(abs_raw)))

  # Column name varies across ABSOLUTE file versions — find sample/array column
  id_col  <- intersect(c("sample", "array", "tcga_sample_id"), names(abs_raw))[1]
  pur_col <- intersect(c("purity", "cancer_dna_fraction"),     names(abs_raw))[1]

  if (is.na(id_col) || is.na(pur_col)) {
    warning("Could not identify ID or purity column in ABSOLUTE file. ",
            "Found: ", paste(names(abs_raw), collapse = ", "))
  } else {
    abs_df <- abs_raw |>
      dplyr::select(ID_abs = !!id_col, purity_absolute = !!pur_col) |>
      dplyr::filter(!is.na(purity_absolute)) |>
      # ABSOLUTE uses 15-char IDs (TCGA-XX-XXXX-01); pipeline uses 16-char
      # (TCGA-XX-XXXX-01A).  Trim to 15 for matching.
      dplyr::mutate(ID15 = substr(as.character(ID_abs), 1, 15))

    cat(sprintf("  ABSOLUTE loaded: %d TCGA samples\n", nrow(abs_df)))
  }
} else {
  cat(sprintf(paste0(
    "  ABSOLUTE file not found at '%s'.\n",
    "  Skipping ABSOLUTE panel.\n",
    "  Download with:\n",
    "    curl -o %s \\\n",
    "      \"https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5\"\n"
  ), ABSOLUTE_PATH, ABSOLUTE_PATH))
}

##############################################################################
# 6.  Merge all sources
##############################################################################
cat("Merging...\n")

df <- score_df |>
  dplyr::left_join(
    est_df |> dplyr::select(ID, purity_estimate, stromal, immune),
    by = "ID"
  )

if (!is.null(abs_df)) {
  df <- df |>
    dplyr::mutate(
      ID15 = dplyr::if_else(dataset == "TCGA_PAAD",
                            substr(ID, 1, 15), NA_character_)
    ) |>
    dplyr::left_join(abs_df |> dplyr::select(ID15, purity_absolute),
                     by = "ID15") |>
    dplyr::select(-ID15)
}

cat(sprintf("  Rows: %d | ESTIMATE: %d covered | ABSOLUTE: %s\n",
    nrow(df),
    sum(!is.na(df$purity_estimate)),
    if (!is.null(abs_df))
      sprintf("%d / %d TCGA covered",
              sum(!is.na(df$purity_absolute)),
              sum(df$dataset == "TCGA_PAAD"))
    else "not loaded"))

##############################################################################
# 7.  Spearman correlations
##############################################################################
cat("\n--- Spearman correlations: Factor 2 vs purity ---\n")

spearman_summary <- function(x, y, label) {
  ok  <- !is.na(x) & !is.na(y)
  cr  <- cor.test(x[ok], y[ok], method = "spearman")
  cat(sprintf("  %-35s  rho = %+.3f  p = %.3g  n = %d\n",
              label, cr$estimate, cr$p.value, sum(ok)))
  invisible(list(rho = cr$estimate, pval = cr$p.value, n = sum(ok), label = label))
}

cat("ESTIMATE purity:\n")
cor_est_pooled <- spearman_summary(df$factor2, df$purity_estimate, "Pooled (TCGA + CPTAC)")
for (ds in c("TCGA_PAAD", "CPTAC")) {
  sub <- df[df$dataset == ds, ]
  spearman_summary(sub$factor2, sub$purity_estimate, ds)
}

if (!is.null(abs_df)) {
  cat("ABSOLUTE purity (TCGA only):\n")
  tcga_sub <- df[df$dataset == "TCGA_PAAD", ]
  cor_abs  <- spearman_summary(tcga_sub$factor2, tcga_sub$purity_absolute, "TCGA_PAAD")
}

##############################################################################
# 8.  Figures
##############################################################################
cat("\nGenerating figures...\n")

cohort_colors <- c("TCGA_PAAD" = "#1f78b4", "CPTAC" = "#e31a1c")

make_scatter <- function(data, purity_col, y_lab) {
  sub <- data |> dplyr::filter(!is.na(.data[[purity_col]]))
  cr  <- cor.test(sub$factor2, sub[[purity_col]], method = "spearman")
  lbl <- sprintf("Spearman \u03c1 = %+.2f\np = %.2g\nn = %d",
                 cr$estimate, cr$p.value, nrow(sub))

  ggplot(sub, aes(x = factor2, y = .data[[purity_col]], color = dataset)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
             label = lbl, size = 3.0) +
    scale_color_manual(values = cohort_colors, drop = FALSE) +
    labs(
      x     = "Factor 2 score  (DeSurv k = 5,  W\u1D40X)",
      y     = y_lab,
      color = "Cohort"
    ) +
    theme_classic(base_size = 10)
}

p_est <- make_scatter(df, "purity_estimate",
                      "Tumor purity (ESTIMATE)")

if (!is.null(abs_df) && "purity_absolute" %in% names(df)) {
  p_abs <- make_scatter(
    df |> dplyr::filter(dataset == "TCGA_PAAD"),
    "purity_absolute",
    "Tumor purity (ABSOLUTE)"
  )
  fig       <- cowplot::plot_grid(p_est, p_abs, ncol = 2, labels = c("A", "B"))
  fig_width <- 9
} else {
  fig       <- p_est
  fig_width <- 5
}

out_path <- file.path(OUT_DIR, "factor2_purity_correlation.pdf")
ggsave(out_path, fig, width = fig_width, height = 4)
cat(sprintf("Saved: %s\n", out_path))
