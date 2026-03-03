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
  library(survival)
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
  dplyr::left_join(sinfo,
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

# Extract PurIST and DeCAF for joining — available via load_data_internal subtype CSV
purist_df <- dplyr::bind_rows(
  tcga$sampInfo  |> dplyr::select(ID, PurIST, DeCAF),
  cptac$sampInfo |> dplyr::select(ID, PurIST, DeCAF)
)

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
  ) |>
  dplyr::left_join(purist_df, by = "ID")

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

# Diagnostic: confirm join populated subtype columns
cat(sprintf("  Columns in df: %s\n", paste(names(df), collapse = ", ")))
cat(sprintf("  PurIST non-NA: %d / %d\n", sum(!is.na(df$PurIST)), nrow(df)))
cat(sprintf("  DeCAF  non-NA: %d / %d\n", sum(!is.na(df$DeCAF)),  nrow(df)))

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

make_scatter(df%>%filter(dataset=="TCGA_PAAD"), "purity_estimate",
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

##############################################################################
# 9.  Boxplots: Factor 2 score by PurIST subtype
##############################################################################
cat("\nGenerating PurIST boxplots...\n")

df_purist <- df |>
  dplyr::filter(!is.na(PurIST), PurIST %in% c("Classical", "Basal-like"))

# Wilcoxon test per dataset and pooled
for (grp in c("Pooled", unique(df_purist$dataset))) {
  sub <- if (grp == "Pooled") df_purist else df_purist |> dplyr::filter(dataset == grp)
  wt  <- wilcox.test(factor2 ~ PurIST, data = sub)
  cat(sprintf("  %-12s  W = %.0f  p = %.3g  n(Classical) = %d  n(Basal-like) = %d\n",
      grp, wt$statistic, wt$p.value,
      sum(sub$PurIST == "Classical"),
      sum(sub$PurIST == "Basal-like")))
}

subtype_colors <- c("Classical" = "#2166ac", "Basal-like" = "#d6604d")

make_boxplot <- function(data, title = NULL) {
  wt  <- wilcox.test(factor2 ~ PurIST, data = data)
  lbl <- sprintf("Wilcoxon p = %.2g\nn = %d", wt$p.value, nrow(data))

  ggplot(data, aes(x = PurIST, y = factor2, fill = PurIST)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
    annotate("text", x = 1.5, y = Inf, vjust = 1.4,
             label = lbl, size = 3.0) +
    scale_fill_manual(values = subtype_colors) +
    labs(
      x     = "PurIST subtype",
      y     = "Factor 2 score  (DeSurv k = 5,  W\u1D40X)",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
}

p_box_pooled <- make_boxplot(df_purist, title = "Pooled")
p_box_tcga   <- make_boxplot(df_purist |> dplyr::filter(dataset == "TCGA_PAAD"),
                             title = "TCGA_PAAD")
p_box_cptac  <- make_boxplot(df_purist |> dplyr::filter(dataset == "CPTAC"),
                             title = "CPTAC")

fig_box <- cowplot::plot_grid(p_box_pooled, p_box_tcga, p_box_cptac,
                              ncol = 3, labels = c("A", "B", "C"))

box_path <- file.path(OUT_DIR, "factor2_purist_boxplot.pdf")
ggsave(box_path, fig_box, width = 9, height = 4)
cat(sprintf("Saved: %s\n", box_path))

##############################################################################
# 10.  Boxplots: Factor 2 score by DeCAF subtype
##############################################################################
cat("\nGenerating DeCAF boxplots...\n")

if (!"DeCAF" %in% names(df)) {
  stop("DeCAF column missing from df — check that purist_df was created with ",
       "dplyr::select(ID, PurIST, DeCAF) at step 3 and that the join succeeded.")
}
if (all(is.na(df$DeCAF))) {
  warning("DeCAF is all-NA in df — the subtype CSV join likely failed (ID mismatch).\n",
          "Skipping DeCAF boxplots.")
  quit(status = 0)
}

df_decaf <- df |> dplyr::filter(!is.na(DeCAF))
decaf_levels <- sort(unique(df_decaf$DeCAF))
cat(sprintf("  DeCAF categories: %s\n", paste(decaf_levels, collapse = ", ")))

# Wilcoxon test (pairwise if >2 groups, otherwise single test)
for (grp in c("Pooled", unique(df_decaf$dataset))) {
  sub <- if (grp == "Pooled") df_decaf else df_decaf |> dplyr::filter(dataset == grp)
  if (nrow(sub) == 0 || length(unique(sub$DeCAF)) < 2) {
    cat(sprintf("  %-12s  skipped (n=%d, levels=%d)\n",
                grp, nrow(sub), length(unique(sub$DeCAF))))
    next
  }
  if (length(unique(sub$DeCAF)) == 2) {
    wt  <- wilcox.test(factor2 ~ DeCAF, data = sub)
    cat(sprintf("  %-12s  W = %.0f  p = %.3g", grp, wt$statistic, wt$p.value))
  } else {
    kt <- kruskal.test(factor2 ~ DeCAF, data = sub)
    cat(sprintf("  %-12s  Kruskal-Wallis p = %.3g", grp, kt$p.value))
  }
  counts <- table(sub$DeCAF)
  cat(sprintf("  n = %s\n", paste(names(counts), counts, sep = ":", collapse = "  ")))
}

# Use a colorblind-friendly palette scaled to however many DeCAF groups exist
decaf_colors <- setNames(
  scales::hue_pal()(length(decaf_levels)),
  decaf_levels
)

make_boxplot_cat <- function(data, fill_var, title = NULL) {
  lvls <- sort(unique(data[[fill_var]]))
  if (length(lvls) == 2) {
    wt  <- wilcox.test(as.formula(paste("factor2 ~", fill_var)), data = data)
    lbl <- sprintf("Wilcoxon p = %.2g\nn = %d", wt$p.value, nrow(data))
  } else {
    kt  <- kruskal.test(as.formula(paste("factor2 ~", fill_var)), data = data)
    lbl <- sprintf("Kruskal-Wallis p = %.2g\nn = %d", kt$p.value, nrow(data))
  }

  ggplot(data, aes(x = .data[[fill_var]], y = factor2, fill = .data[[fill_var]])) +
    geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
    annotate("text", x = (length(lvls) + 1) / 2, y = Inf, vjust = 1.4,
             label = lbl, size = 3.0) +
    scale_fill_manual(values = decaf_colors) +
    labs(
      x     = "DeCAF subtype",
      y     = "Factor 2 score  (DeSurv k = 5,  W\u1D40X)",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

p_decaf_pooled <- make_boxplot_cat(df_decaf, "DeCAF", title = "Pooled")
p_decaf_tcga   <- make_boxplot_cat(df_decaf |> dplyr::filter(dataset == "TCGA_PAAD"),
                                   "DeCAF", title = "TCGA_PAAD")
p_decaf_cptac  <- make_boxplot_cat(df_decaf |> dplyr::filter(dataset == "CPTAC"),
                                   "DeCAF", title = "CPTAC")

fig_decaf <- cowplot::plot_grid(p_decaf_pooled, p_decaf_tcga, p_decaf_cptac,
                                ncol = 3, labels = c("A", "B", "C"))

decaf_path <- file.path(OUT_DIR, "factor2_decaf_boxplot.pdf")
ggsave(decaf_path, fig_decaf, width = 9, height = 4)
cat(sprintf("Saved: %s\n", decaf_path))

##############################################################################
# 11.  Load pre-processed validation data from targets store
#
#      data_val_filtered_elbowk_tcgacptac is a named list of datasets that
#      have already been:
#        - filtered to the same gene set as training (tar_data_filtered_elbowk)
#        - transformed via DeSurv::preprocess_data (same method_trans_train and
#          transform_target as training) so the expression scale matches W
#        - restricted to valid-survival samples via samp_keeps
##############################################################################
cat("\n\n==== VALIDATION DATASETS ====\n")
cat("Loading pre-processed validation data from targets store...\n")

val_data_list <- tar_read(data_val_filtered_elbowk_tcgacptac, store = STORE)
cat(sprintf("  Datasets: %s\n", paste(names(val_data_list), collapse = ", ")))
for (ds in names(val_data_list)) {
  cat(sprintf("  %-22s: %d samples, %d genes\n",
              ds, ncol(val_data_list[[ds]]$ex), nrow(val_data_list[[ds]]$ex)))
}

# CPTAC is also part of the tcgacptac training set — flag it so results are
# interpreted correctly (training performance, not external validation).
TRAINING_DATASETS <- c("TCGA_PAAD", "CPTAC")

##############################################################################
# 12.  Apply full W matrix to each validation dataset
#
#      Scores = W^T X  (all k factors).  Genes absent from a dataset are
#      zero-filled, consistent with preprocess_validation_data(zero_fill=TRUE).
##############################################################################
cat("\nApplying W matrix (all factors) to validation datasets...\n")

W_all <- tar_fit$W          # genes x k
k_all <- ncol(W_all)
W_gene_order <- rownames(W_all)

project_W <- function(W, X) {
  # Returns W^T %*% X after aligning genes; missing genes zero-filled.
  wg      <- rownames(W)
  common  <- intersect(wg, rownames(X))
  missing <- setdiff(wg, rownames(X))
  if (length(missing)) {
    Z <- matrix(0, nrow = length(missing), ncol = ncol(X),
                dimnames = list(missing, colnames(X)))
    X <- rbind(X[common, , drop = FALSE], Z)
  }
  t(W[wg, , drop = FALSE]) %*% X[wg, , drop = FALSE]  # k x n
}

val_h_list <- lapply(names(val_data_list), function(ds) {
  X_val  <- val_data_list[[ds]]$ex
  common <- length(intersect(W_gene_order, rownames(X_val)))
  n_miss <- length(W_gene_order) - common
  cat(sprintf("  %-22s  %d samples  %d/%d genes matched  %d zero-filled\n",
              ds, ncol(X_val), common, length(W_gene_order), n_miss))

  H    <- project_W(W_all, X_val)           # k x n
  df_h <- as.data.frame(t(H), stringsAsFactors = FALSE)
  colnames(df_h) <- paste0("factor", seq_len(k_all))
  data.frame(ID = colnames(X_val), dataset = ds, df_h,
             stringsAsFactors = FALSE)
})

val_h_df <- dplyr::bind_rows(val_h_list)
cat(sprintf("  Total: %d validation samples\n", nrow(val_h_df)))

# Save full factor score matrix for downstream use
results_dir <- file.path("results", "factor2_purity")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(val_h_df, file.path(results_dir, "val_factor_scores.rds"))
cat(sprintf("  Saved all-factor scores: %s/val_factor_scores.rds\n", results_dir))

##############################################################################
# 13.  Factor 2 scores for validation + ESTIMATE purity
#
#      ESTIMATE needs unfiltered raw expression.  We reload each dataset via
#      load_data_internal(), restrict to the exact sample IDs used in
#      val_data_list (which went through samp_keeps + survival filtering), and
#      run the same run_estimate() helper defined in Step 4.
##############################################################################
cat("\nComputing Factor 2 purity correlation (validation datasets)...\n")

val_f2_df <- val_h_df |> dplyr::select(ID, dataset, factor2)

cat("  Running ESTIMATE on validation datasets (raw expression)...\n")

est_val_list <- lapply(names(val_data_list), function(ds) {
  tryCatch({
    d        <- load_data_internal(ds)
    # Filter to the exact IDs present in the pipeline-processed val data
    val_ids  <- colnames(val_data_list[[ds]]$ex)
    avail    <- intersect(val_ids, colnames(d$ex))
    if (!length(avail)) {
      message(sprintf("    %s: no matching IDs between raw data and val target", ds))
      return(NULL)
    }
    ex_use <- d$ex[, avail, drop = FALSE]
    est    <- run_estimate(ex_use, ds)
    estimate_to_purity(est)
  }, error = function(e) {
    message(sprintf("    ESTIMATE failed for %s: %s", ds, conditionMessage(e)))
    NULL
  })
})

est_val_df <- dplyr::bind_rows(Filter(Negate(is.null), est_val_list)) |>
  dplyr::rename(ID = sample) |>
  dplyr::select(ID, dataset, purity_estimate)

cat(sprintf("  ESTIMATE: %d / %d validation samples covered\n",
            nrow(est_val_df), nrow(val_f2_df)))

# Extract PurIST / DeCAF from val sampInfo (available for datasets that have
# a *_subtype.csv with those columns; NA for datasets that don't).
purist_val_list <- lapply(names(val_data_list), function(ds) {
  si <- val_data_list[[ds]]$sampInfo
  if (is.null(si) || !nrow(si)) return(NULL)
  if (!"ID" %in% names(si)) {
    if (!is.null(rownames(si))) si$ID <- rownames(si) else return(NULL)
  }
  cols <- intersect(c("ID", "PurIST", "DeCAF"), names(si))
  if (length(cols) < 2) return(NULL)
  si[, cols, drop = FALSE]
})
purist_val_df <- dplyr::bind_rows(Filter(Negate(is.null), purist_val_list))

# Merge
val_df2 <- val_f2_df |>
  dplyr::left_join(est_val_df |> dplyr::select(ID, purity_estimate), by = "ID")
if (nrow(purist_val_df) > 0) {
  val_df2 <- dplyr::left_join(val_df2, purist_val_df, by = "ID")
} else {
  val_df2$PurIST <- NA_character_
  val_df2$DeCAF  <- NA_character_
}

cat(sprintf("  Merged: %d rows | ESTIMATE non-NA: %d | PurIST non-NA: %d | DeCAF non-NA: %d\n",
    nrow(val_df2),
    sum(!is.na(val_df2$purity_estimate)),
    sum(!is.na(val_df2$PurIST)),
    sum(!is.na(val_df2$DeCAF))))

##############################################################################
# 14.  Spearman correlations: Factor 2 vs ESTIMATE purity (validation)
##############################################################################
cat("\n--- Spearman correlations: Factor 2 vs purity (validation) ---\n")

cat("ESTIMATE purity (pooled validation):\n")
spearman_summary(val_df2$factor2, val_df2$purity_estimate, "All validation (pooled)")
for (ds in names(val_data_list)) {
  sub <- val_df2[val_df2$dataset == ds, ]
  tag <- if (ds %in% TRAINING_DATASETS) paste0(ds, " [TRAINING]") else ds
  spearman_summary(sub$factor2, sub$purity_estimate, tag)
}

##############################################################################
# 15.  Scatter plots: Factor 2 vs purity (validation)
##############################################################################
cat("\nGenerating validation scatter plots...\n")

# Build a per-dataset colour palette (extend cohort_colors for all val datasets)
all_ds      <- unique(val_df2$dataset)
extra_ds    <- setdiff(all_ds, names(cohort_colors))
extra_cols  <- setNames(
  scales::hue_pal()(length(extra_ds)),
  extra_ds
)
val_colors  <- c(cohort_colors, extra_cols)

# Helper: scatter for validation (same structure as make_scatter, but labels
# each dataset separately and uses val_colors).
make_scatter_val <- function(data, purity_col, y_lab, title = NULL) {
  sub <- data |> dplyr::filter(!is.na(.data[[purity_col]]))
  if (!nrow(sub)) {
    return(ggplot() + labs(title = paste(title, "(no data)")) + theme_void())
  }
  cr  <- cor.test(sub$factor2, sub[[purity_col]], method = "spearman")
  lbl <- sprintf("Spearman \u03c1 = %+.2f\np = %.2g\nn = %d",
                 cr$estimate, cr$p.value, nrow(sub))

  ggplot(sub, aes(x = factor2, y = .data[[purity_col]], color = dataset)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
             label = lbl, size = 3.0) +
    scale_color_manual(values = val_colors, drop = TRUE) +
    labs(
      x     = "Factor 2 score  (DeSurv k = 5,  W\u1D40X)",
      y     = y_lab,
      color = "Dataset",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "right")
}

# Pooled scatter
p_val_pooled <- make_scatter_val(val_df2, "purity_estimate",
                                 "Tumor purity (ESTIMATE)",
                                 "All validation (pooled)")

# Per-dataset scatters
p_val_each <- lapply(all_ds, function(ds) {
  make_scatter_val(val_df2 |> dplyr::filter(dataset == ds),
                   "purity_estimate", "Tumor purity (ESTIMATE)", ds)
})
names(p_val_each) <- all_ds

# Combined figure: pooled on top row, per-dataset on subsequent rows
n_each  <- length(p_val_each)
n_cols  <- min(3L, n_each)
fig_val <- cowplot::plot_grid(
  p_val_pooled,
  cowplot::plot_grid(plotlist = p_val_each, ncol = n_cols),
  ncol = 1, rel_heights = c(1, ceiling(n_each / n_cols))
)

val_scatter_path <- file.path(OUT_DIR, "factor2_val_purity_scatter.pdf")
ggsave(val_scatter_path, fig_val,
       width = 4 * n_cols, height = 4 * (1 + ceiling(n_each / n_cols)))
cat(sprintf("Saved: %s\n", val_scatter_path))

##############################################################################
# 16.  Boxplots: Factor 2 by PurIST (validation)
##############################################################################
cat("\nGenerating validation PurIST boxplots...\n")

val_purist <- val_df2 |>
  dplyr::filter(!is.na(PurIST), PurIST %in% c("Classical", "Basal-like"))

if (nrow(val_purist) == 0) {
  cat("  No PurIST annotations available in validation datasets — skipping.\n")
} else {
  ds_with_purist <- unique(val_purist$dataset)
  cat(sprintf("  Datasets with PurIST: %s\n",
              paste(ds_with_purist, collapse = ", ")))

  for (grp in c("Pooled", ds_with_purist)) {
    sub <- if (grp == "Pooled") val_purist else
             val_purist |> dplyr::filter(dataset == grp)
    if (nrow(sub) < 2 || length(unique(sub$PurIST)) < 2) next
    wt  <- wilcox.test(factor2 ~ PurIST, data = sub)
    cat(sprintf("  %-22s  W = %.0f  p = %.3g  n(Classical) = %d  n(Basal-like) = %d\n",
        grp, wt$statistic, wt$p.value,
        sum(sub$PurIST == "Classical"), sum(sub$PurIST == "Basal-like")))
  }

  p_vbox_pooled <- make_boxplot(val_purist, title = "Pooled (validation)")
  p_vbox_each   <- lapply(ds_with_purist, function(ds)
    make_boxplot(val_purist |> dplyr::filter(dataset == ds), title = ds))

  n_vp      <- length(ds_with_purist)
  fig_vpbox <- cowplot::plot_grid(
    plotlist = c(list(p_vbox_pooled), p_vbox_each),
    ncol     = min(3L, n_vp + 1L)
  )
  vbox_path <- file.path(OUT_DIR, "factor2_val_purist_boxplot.pdf")
  ggsave(vbox_path, fig_vpbox,
         width = 3 * min(3L, n_vp + 1L), height = 4)
  cat(sprintf("Saved: %s\n", vbox_path))
}

##############################################################################
# 17.  Boxplots: Factor 2 by DeCAF (validation)
##############################################################################
cat("\nGenerating validation DeCAF boxplots...\n")

if (!"DeCAF" %in% names(val_df2) || all(is.na(val_df2$DeCAF))) {
  cat("  No DeCAF annotations in validation datasets — skipping.\n")
} else {
  val_decaf      <- val_df2 |> dplyr::filter(!is.na(DeCAF))
  ds_with_decaf  <- unique(val_decaf$dataset)
  vdecaf_levels  <- sort(unique(val_decaf$DeCAF))
  cat(sprintf("  Datasets with DeCAF: %s\n",
              paste(ds_with_decaf, collapse = ", ")))
  cat(sprintf("  DeCAF categories: %s\n",
              paste(vdecaf_levels, collapse = ", ")))

  vdecaf_colors <- setNames(scales::hue_pal()(length(vdecaf_levels)), vdecaf_levels)

  for (grp in c("Pooled", ds_with_decaf)) {
    sub <- if (grp == "Pooled") val_decaf else
             val_decaf |> dplyr::filter(dataset == grp)
    if (nrow(sub) == 0 || length(unique(sub$DeCAF)) < 2) {
      cat(sprintf("  %-22s  skipped (n=%d, levels=%d)\n",
                  grp, nrow(sub), length(unique(sub$DeCAF))))
      next
    }
    if (length(unique(sub$DeCAF)) == 2) {
      wt  <- wilcox.test(factor2 ~ DeCAF, data = sub)
      cat(sprintf("  %-22s  W = %.0f  p = %.3g", grp, wt$statistic, wt$p.value))
    } else {
      kt <- kruskal.test(factor2 ~ DeCAF, data = sub)
      cat(sprintf("  %-22s  Kruskal-Wallis p = %.3g", grp, kt$p.value))
    }
    counts <- table(sub$DeCAF)
    cat(sprintf("  n = %s\n", paste(names(counts), counts, sep = ":", collapse = "  ")))
  }

  make_boxplot_decaf_val <- function(data, title = NULL) {
    lvls <- sort(unique(data$DeCAF))
    if (length(lvls) == 2) {
      wt  <- wilcox.test(factor2 ~ DeCAF, data = data)
      lbl <- sprintf("Wilcoxon p = %.2g\nn = %d", wt$p.value, nrow(data))
    } else {
      kt  <- kruskal.test(factor2 ~ DeCAF, data = data)
      lbl <- sprintf("Kruskal-Wallis p = %.2g\nn = %d", kt$p.value, nrow(data))
    }
    ggplot(data, aes(x = DeCAF, y = factor2, fill = DeCAF)) +
      geom_boxplot(outlier.size = 0.8, width = 0.5, alpha = 0.8) +
      geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
      annotate("text", x = (length(lvls) + 1) / 2, y = Inf, vjust = 1.4,
               label = lbl, size = 3.0) +
      scale_fill_manual(values = vdecaf_colors) +
      labs(x = "DeCAF subtype",
           y = "Factor 2 score  (DeSurv k = 5,  W\u1D40X)",
           title = title) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1))
  }

  p_vdecaf_pooled <- make_boxplot_decaf_val(val_decaf, title = "Pooled (validation)")
  p_vdecaf_each   <- lapply(ds_with_decaf, function(ds)
    make_boxplot_decaf_val(val_decaf |> dplyr::filter(dataset == ds), title = ds))

  n_vd        <- length(ds_with_decaf)
  fig_vdecaf  <- cowplot::plot_grid(
    plotlist = c(list(p_vdecaf_pooled), p_vdecaf_each),
    ncol     = min(3L, n_vd + 1L)
  )
  vdecaf_path <- file.path(OUT_DIR, "factor2_val_decaf_boxplot.pdf")
  ggsave(vdecaf_path, fig_vdecaf,
         width = 3 * min(3L, n_vd + 1L), height = 4)
  cat(sprintf("Saved: %s\n", vdecaf_path))
}

##############################################################################
# 18.  CoxPH: Factor 2 score in validation data, unadjusted and adjusted
#      for PurIST / DeCAF
#
#      Survival time and event come from val_data_list[[ds]]$sampInfo, which
#      is preserved through preprocess_validation_data() and contains the same
#      time / event columns as the original survival RDS.
#
#      Models fit per-dataset and pooled:
#        (a) Univariate:          Surv(time, event) ~ factor2
#        (b) Adjusted for PurIST: Surv(time, event) ~ factor2 + PurIST
#        (c) Adjusted for DeCAF:  Surv(time, event) ~ factor2 + DeCAF
#        (d) Adjusted for both:   Surv(time, event) ~ factor2 + PurIST + DeCAF
#      Pooled models stratify on dataset to allow separate baseline hazards.
##############################################################################
cat("\n\n==== COXPH: Factor 2 in validation data ====\n")

# ---- 18a.  Collect survival info from sampInfo ----
surv_val_list <- lapply(names(val_data_list), function(ds) {
  vd <- val_data_list[[ds]]
  si <- vd$sampInfo
  if (is.null(si) || !nrow(si)) {
    # Fall back to $y / $d slots if sampInfo is absent
    return(data.frame(ID = colnames(vd$ex), time = vd$y, event = vd$d,
                      stringsAsFactors = FALSE))
  }
  if (!"ID" %in% names(si)) si$ID <- rownames(si)
  # time / event are set by DeSurv::preprocess_data; grab whichever name present
  t_col <- intersect(c("time", "y"), names(si))[1]
  e_col <- intersect(c("event", "d"), names(si))[1]
  if (is.na(t_col) || is.na(e_col)) {
    if (!is.null(vd$y) && !is.null(vd$d)) {
      return(data.frame(ID = si$ID, time = vd$y, event = vd$d,
                        stringsAsFactors = FALSE))
    }
    message(sprintf("  %s: cannot locate time/event — skipping", ds))
    return(NULL)
  }
  data.frame(ID = si$ID, time = si[[t_col]], event = si[[e_col]],
             stringsAsFactors = FALSE)
})

surv_val_df <- dplyr::bind_rows(Filter(Negate(is.null), surv_val_list))

# Join survival onto the val_df2 frame (which already has factor2, PurIST, DeCAF)
val_cox_df <- val_df2 |>
  dplyr::left_join(surv_val_df, by = "ID") |>
  dplyr::filter(!is.na(time), !is.na(event), time > 0)

cat(sprintf("  Samples with valid survival data: %d / %d\n",
            nrow(val_cox_df), nrow(val_df2)))
cat(sprintf("  PurIST available: %d  |  DeCAF available: %d\n",
            sum(!is.na(val_cox_df$PurIST)), sum(!is.na(val_cox_df$DeCAF))))

# ---- 18b.  Fitting helper ----
# Returns a one-row data frame with HR / 95% CI / p for the factor2 term.
fit_cox_f2 <- function(data, formula_str, label) {
  tryCatch({
    fit   <- survival::coxph(as.formula(formula_str), data = data)
    coefs <- summary(fit)$coefficients
    row   <- coefs[grep("^factor2$", rownames(coefs)), , drop = FALSE]
    if (!nrow(row)) {
      cat(sprintf("  %-50s  factor2 term not found\n", label))
      return(NULL)
    }
    hr  <- exp(row[1, "coef"])
    se  <- row[1, "se(coef)"]
    lo  <- exp(row[1, "coef"] - 1.96 * se)
    hi  <- exp(row[1, "coef"] + 1.96 * se)
    pv  <- row[1, "Pr(>|z|)"]
    cat(sprintf("  %-50s  HR = %.3f (%.3f-%.3f)  p = %.3g  n = %d\n",
                label, hr, lo, hi, pv, fit$n))
    data.frame(label = label, HR = hr, CI_lo = lo, CI_hi = hi,
               pval = pv, n = fit$n, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat(sprintf("  %-50s  FAILED: %s\n", label, conditionMessage(e)))
    NULL
  })
}

# ---- 18c.  Pooled models (strata(dataset) for separate baseline hazards) ----
cat("\nPooled validation (strata by dataset):\n")

pooled_models <- list(
  list(adj = "univariate",
       form = "survival::Surv(time, event) ~ factor2 + strata(dataset)"),
  list(adj = "adjusted PurIST",
       form = "survival::Surv(time, event) ~ factor2 + PurIST + strata(dataset)"),
  list(adj = "adjusted DeCAF",
       form = "survival::Surv(time, event) ~ factor2 + DeCAF + strata(dataset)"),
  list(adj = "adjusted PurIST + DeCAF",
       form = "survival::Surv(time, event) ~ factor2 + PurIST + DeCAF + strata(dataset)")
)

cox_pooled <- dplyr::bind_rows(lapply(pooled_models, function(m) {
  # Subset to samples with the required covariates non-NA
  required <- character(0)
  if (grepl("PurIST", m$form)) required <- c(required, "PurIST")
  if (grepl("DeCAF",  m$form)) required <- c(required, "DeCAF")
  sub <- val_cox_df
  for (col in required) sub <- sub[!is.na(sub[[col]]), , drop = FALSE]

  res <- fit_cox_f2(sub, m$form, paste0("Pooled — ", m$adj))
  if (!is.null(res)) {
    res$dataset   <- "Pooled"
    res$adj_model <- m$adj
  }
  res
}))

# ---- 18d.  Per-dataset models ----
ds_models <- list(
  list(adj = "univariate",
       form = "survival::Surv(time, event) ~ factor2"),
  list(adj = "adjusted PurIST",
       form = "survival::Surv(time, event) ~ factor2 + PurIST"),
  list(adj = "adjusted DeCAF",
       form = "survival::Surv(time, event) ~ factor2 + DeCAF"),
  list(adj = "adjusted PurIST + DeCAF",
       form = "survival::Surv(time, event) ~ factor2 + PurIST + DeCAF")
)

cox_per_ds <- dplyr::bind_rows(lapply(names(val_data_list), function(ds) {
  sub_ds <- val_cox_df |> dplyr::filter(dataset == ds)
  cat(sprintf("\n%s (n = %d):\n", ds, nrow(sub_ds)))

  dplyr::bind_rows(lapply(ds_models, function(m) {
    required <- character(0)
    if (grepl("PurIST", m$form)) required <- c(required, "PurIST")
    if (grepl("DeCAF",  m$form)) required <- c(required, "DeCAF")
    sub <- sub_ds
    for (col in required) sub <- sub[!is.na(sub[[col]]), , drop = FALSE]
    if (nrow(sub) < 5) {
      cat(sprintf("  %-50s  skipped (n = %d)\n",
                  paste0(ds, " — ", m$adj), nrow(sub)))
      return(NULL)
    }

    res <- fit_cox_f2(sub, m$form, paste0(ds, " — ", m$adj))
    if (!is.null(res)) {
      res$dataset   <- ds
      res$adj_model <- m$adj
    }
    res
  }))
}))

# ---- 18e.  Save combined results table ----
cox_all <- dplyr::bind_rows(cox_pooled, cox_per_ds) |>
  dplyr::select(dataset, adj_model, HR, CI_lo, CI_hi, pval, n, label)

saveRDS(cox_all, file.path(results_dir, "val_cox_factor2.rds"))
cat(sprintf("\nSaved Cox results: %s/val_cox_factor2.rds\n", results_dir))

# ---- 18f.  Forest plot: univariate HR across datasets ----
cat("\nGenerating Cox forest plot...\n")

cox_forest_df <- cox_all |>
  dplyr::filter(adj_model == "univariate") |>
  dplyr::mutate(
    dataset = factor(dataset, levels = rev(c("Pooled", names(val_data_list)))),
    sig     = ifelse(pval < 0.05, "p < 0.05", "p \u2265 0.05")
  )

p_forest_uni <- ggplot(cox_forest_df,
                       aes(x = HR, y = dataset, color = sig, shape = sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("p < 0.05" = "#d73027", "p \u2265 0.05" = "#4393c3")) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p \u2265 0.05" = 1)) +
  labs(x     = "Hazard ratio (Factor 2 score, per unit)",
       y     = NULL,
       title = "CoxPH — Factor 2 (univariate, validation datasets)",
       color = NULL, shape = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

# ---- 18g.  Multi-panel forest: unadjusted vs adjusted ----
make_forest_panel <- function(adj) {
  sub <- cox_all |>
    dplyr::filter(adj_model == adj) |>
    dplyr::mutate(
      dataset = factor(dataset, levels = rev(c("Pooled", names(val_data_list)))),
      sig     = ifelse(pval < 0.05, "p < 0.05", "p \u2265 0.05")
    )
  if (!nrow(sub)) return(ggplot() + labs(title = paste(adj, "(no data)")) + theme_void())

  ggplot(sub, aes(x = HR, y = dataset, color = sig, shape = sig)) +
    geom_point(size = 2.5) +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("p < 0.05" = "#d73027", "p \u2265 0.05" = "#4393c3"),
                       drop = FALSE) +
    scale_shape_manual(values = c("p < 0.05" = 16, "p \u2265 0.05" = 1),
                       drop = FALSE) +
    labs(x = "HR (Factor 2)", y = NULL,
         title = adj, color = NULL, shape = NULL) +
    theme_classic(base_size = 9) +
    theme(legend.position = "none")
}

adj_labels <- c("univariate", "adjusted PurIST", "adjusted DeCAF",
                "adjusted PurIST + DeCAF")
forest_panels <- lapply(adj_labels, make_forest_panel)
legend_grob   <- cowplot::get_legend(
  p_forest_uni + theme(legend.position = "bottom")
)

fig_forest <- cowplot::plot_grid(
  cowplot::plot_grid(plotlist = forest_panels, ncol = 2),
  legend_grob,
  ncol = 1, rel_heights = c(1, 0.08)
)

n_ds       <- length(val_data_list) + 1L   # +1 for Pooled
forest_path <- file.path(OUT_DIR, "factor2_val_cox_forest.pdf")
ggsave(forest_path, fig_forest, width = 10,
       height = 1.5 + 0.45 * n_ds)
cat(sprintf("Saved: %s\n", forest_path))
