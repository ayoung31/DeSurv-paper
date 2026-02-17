##############################################################################
# cv_grid_validation_analysis.R
#
# Training + external validation analysis of 9 cv_grid DeSurv fits.
# Extends cv_grid_training_analysis.R with:
#   - Validation Cox models (unadjusted + adjusted for PurIST + DeCAF)
#   - Likelihood ratio tests (LRT) vs PurIST + DeCAF baseline
#   - Median-split Kaplan-Meier analysis
#
# Run from the DeSurv-paper repo root:
#   Rscript inst/cv_grid_validation_analysis.R
#
# Validation cohort: 4 independent datasets (n=570 pooled)
#   - Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69)
# All validation Cox models use strata(dataset) to account for cohort-specific
# baseline hazards.
#
# See: docs/plans/2026-02-17-k-sensitivity-synthesis.md for interpretation.
##############################################################################

suppressMessages({
  library(survival)
  library(data.table)
})

STORE <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"

# Load reference data
load("data/derv/cmbSubtypes_formatted.RData")

orig_bundle <- readRDS(file.path(STORE, "objects/val_run_bundle_tcgacptac"))
orig_fit <- orig_bundle$fit_desurv
orig_W <- orig_fit$W
orig_beta <- orig_fit$beta

bo_bundle <- readRDS(file.path(STORE, "objects/bo_bundle_selected_tcgacptac"))
X_train <- bo_bundle$data_filtered$ex
surv_train <- bo_bundle$data_filtered$sampInfo

# Load validation data (named list: Dijk, Moffitt_GEO_array, Puleo_array, PACA_AU)
val_data <- readRDS(file.path(STORE, "objects/data_val_filtered_surv_tcgacptac"))
cat("Validation datasets:", paste(names(val_data), collapse=", "), "\n")
total_val <- sum(sapply(val_data, function(x) ncol(x$ex)))
cat("Total validation n:", total_val, "\n\n")

# Load training subtypes from BOTH TCGA and CPTAC
# Note: bo_bundle$data_filtered$sampInfo only has time/event/dataset,
# not PurIST/DeCAF. We need to load those separately.
tcga_sub <- read.csv("data/original/TCGA_PAAD_subtype.csv")
cptac_sub <- read.csv("data/original/CPTAC_subtype.csv")
all_train_sub <- rbind(
  tcga_sub[, c("ID", "PurIST", "DeCAF")],
  cptac_sub[, c("ID", "PurIST", "DeCAF")]
)
all_train_sub$DeCAF[all_train_sub$DeCAF == "permCAF"] <- "proCAF"
cat("Training subtype labels: n=", nrow(all_train_sub), "\n")
cat("PurIST:", paste(names(table(all_train_sub$PurIST)), table(all_train_sub$PurIST), collapse=", "), "\n")
cat("DeCAF:", paste(names(table(all_train_sub$DeCAF)), table(all_train_sub$DeCAF), collapse=", "), "\n\n")

source("R/get_top_genes.R")

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

cat("All 9 fits loaded.\n\n")

# ============================================================
# Helper: Build pooled validation dataframe for a given W matrix
# Computes H-scores as H = X^T W using all genes in common,
# then standardizes across the pooled validation set.
# ============================================================
build_val_df <- function(W, val_data) {
  pooled <- list()
  for (ds_name in names(val_data)) {
    vd <- val_data[[ds_name]]
    X_val <- vd$ex
    si <- vd$sampInfo

    genes <- intersect(rownames(W), rownames(X_val))
    H_val <- t(X_val[genes, ]) %*% W[genes, ]
    colnames(H_val) <- paste0("Factor", seq_len(ncol(H_val)))

    decaf <- as.character(si$DeCAF)
    decaf[decaf == "permCAF"] <- "proCAF"

    ds_df <- data.frame(
      sample_id = colnames(X_val),
      time = si$time,
      event = si$event,
      PurIST = as.character(si$PurIST),
      DeCAF = decaf,
      dataset = ds_name,
      H_val,
      stringsAsFactors = FALSE
    )

    keep <- complete.cases(ds_df[, c("time", "event")])
    ds_df <- ds_df[keep, ]
    pooled[[ds_name]] <- ds_df
  }
  pooled_df <- do.call(rbind, pooled)
  rownames(pooled_df) <- NULL

  # Standardize H-scores across pooled validation
  k <- sum(grepl("^Factor", names(pooled_df)))
  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)
    pooled_df[[fn]] <- as.numeric(scale(pooled_df[[fn]]))
  }
  pooled_df
}

# ============================================================
# Helper: Build training dataframe for a given W matrix
# ============================================================
build_train_df <- function(W, X_train, surv_train, subtypes) {
  genes <- intersect(rownames(W), rownames(X_train))
  H_train <- t(X_train[genes, ]) %*% W[genes, ]
  colnames(H_train) <- paste0("Factor", seq_len(ncol(H_train)))

  train_ids <- colnames(X_train)
  idx <- match(train_ids, subtypes$ID)

  train_df <- data.frame(
    sample_id = train_ids,
    time = surv_train$time,
    event = surv_train$event,
    PurIST = as.character(subtypes$PurIST[idx]),
    DeCAF = as.character(subtypes$DeCAF[idx]),
    H_train,
    stringsAsFactors = FALSE
  )

  # Standardize H-scores
  k <- sum(grepl("^Factor", names(train_df)))
  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)
    train_df[[fn]] <- as.numeric(scale(train_df[[fn]]))
  }
  train_df
}

# ============================================================
# Run per-fit analysis (training + validation)
# ============================================================
for (nm in names(fits)) {
  fit <- fits[[nm]]
  W <- fit$fit$W
  beta <- fit$fit$beta
  k <- fit$k
  alpha <- fit$alpha

  cat("================================================================\n")
  cat(sprintf("%s (K=%d, alpha=%.2f)\n", nm, k, alpha))
  cat("================================================================\n\n")

  # Identify best iCAF factor (by H-cor with original F1 on training data)
  genes_t <- intersect(rownames(W), rownames(orig_W))
  genes_t <- intersect(genes_t, rownames(X_train))
  H_orig <- t(X_train[genes_t, ]) %*% orig_W[genes_t, ]
  H_new <- t(X_train[genes_t, ]) %*% W[genes_t, ]
  cors <- cor(H_orig[, 1], H_new, method = "spearman")
  best_f <- which.max(cors)
  cat(sprintf("Best iCAF factor: F%d (H-cor=%.3f)\n\n", best_f, max(cors)))

  # Build training and validation dataframes
  train_df <- build_train_df(W, X_train, surv_train, all_train_sub)
  val_df <- build_val_df(W, val_data)

  train_n_full <- sum(complete.cases(train_df[, c("time", "event")]))
  train_n_adj <- sum(complete.cases(train_df[, c("time", "event", "PurIST", "DeCAF")]))
  val_n <- nrow(val_df)
  val_n_adj <- sum(complete.cases(val_df[, c("time", "event", "PurIST", "DeCAF")]))

  cat(sprintf("Training: n=%d (full), n=%d (with PurIST+DeCAF)\n", train_n_full, train_n_adj))
  cat(sprintf("Validation: n=%d (pooled), n=%d (with PurIST+DeCAF)\n", val_n, val_n_adj))
  for (ds in unique(val_df$dataset)) {
    cat(sprintf("    %s: n=%d\n", ds, sum(val_df$dataset == ds)))
  }
  cat("\n")

  # === A. Per-factor Cox: UNADJUSTED ===
  cat("--- A. Per-Factor Cox (Unadjusted) ---\n")
  cat(sprintf("%-8s | %10s %10s %10s | %10s %10s %10s\n",
    "Factor", "Train HR", "Train z", "Train p",
    "Val HR", "Val z", "Val p"))
  cat(paste(rep("-", 85), collapse=""), "\n")

  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)

    # Training (no strata)
    tdf_full <- train_df[complete.cases(train_df[, c("time", "event", fn)]), ]
    cx_t <- coxph(as.formula(paste0("Surv(time, event) ~ ", fn)), data = tdf_full)
    s_t <- summary(cx_t)
    hr_t <- s_t$conf.int[1, 1]
    z_t <- s_t$coefficients[1, "z"]
    p_t <- s_t$coefficients[1, "Pr(>|z|)"]

    # Validation (with strata(dataset))
    cx_v <- coxph(as.formula(paste0("Surv(time, event) ~ ", fn, " + strata(dataset)")),
      data = val_df)
    s_v <- summary(cx_v)
    hr_v <- s_v$conf.int[1, 1]
    z_v <- s_v$coefficients[1, "z"]
    p_v <- s_v$coefficients[1, "Pr(>|z|)"]

    sig_t <- ifelse(p_t < 0.001, "***", ifelse(p_t < 0.01, "**", ifelse(p_t < 0.05, "*", "")))
    sig_v <- ifelse(p_v < 0.001, "***", ifelse(p_v < 0.01, "**", ifelse(p_v < 0.05, "*", "")))

    marker <- ifelse(j == best_f, " <-- iCAF", "")
    cat(sprintf("F%-7d | %10.3f %10.3f %10.4f%s | %10.3f %10.3f %10.4f%s%s\n",
      j, hr_t, z_t, p_t, sig_t, hr_v, z_v, p_v, sig_v, marker))
  }
  cat("\n")

  # === B. Per-factor Cox: ADJUSTED for PurIST + DeCAF ===
  cat("--- B. Per-Factor Cox (Adjusted for PurIST + DeCAF) ---\n")

  adj_cols <- c("time", "event", "PurIST", "DeCAF", paste0("Factor", seq_len(k)))
  train_adj <- train_df[complete.cases(train_df[, adj_cols]), ]
  val_adj <- val_df[complete.cases(val_df[, c("time", "event", "PurIST", "DeCAF")]), ]
  cat(sprintf("  (Training n=%d, Validation n=%d for adjusted models)\n", nrow(train_adj), nrow(val_adj)))

  cat(sprintf("%-8s | %10s %10s %10s | %10s %10s %10s\n",
    "Factor", "Train HR", "Train z", "Train p",
    "Val HR", "Val z", "Val p"))
  cat(paste(rep("-", 85), collapse=""), "\n")

  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)

    fml_t <- as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", fn))
    fml_v <- as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", fn, " + strata(dataset)"))

    cx_t <- tryCatch(coxph(fml_t, data = train_adj), error = function(e) NULL)
    cx_v <- tryCatch(coxph(fml_v, data = val_adj), error = function(e) NULL)

    if (!is.null(cx_t) && !is.null(cx_v)) {
      s_t <- summary(cx_t)
      row_t <- grep(paste0("^", fn, "$"), rownames(s_t$coefficients))
      hr_t <- exp(s_t$coefficients[row_t, "coef"])
      z_t <- s_t$coefficients[row_t, "z"]
      p_t <- s_t$coefficients[row_t, "Pr(>|z|)"]

      s_v <- summary(cx_v)
      row_v <- grep(paste0("^", fn, "$"), rownames(s_v$coefficients))
      hr_v <- exp(s_v$coefficients[row_v, "coef"])
      z_v <- s_v$coefficients[row_v, "z"]
      p_v <- s_v$coefficients[row_v, "Pr(>|z|)"]

      sig_t <- ifelse(p_t < 0.001, "***", ifelse(p_t < 0.01, "**", ifelse(p_t < 0.05, "*", "")))
      sig_v <- ifelse(p_v < 0.001, "***", ifelse(p_v < 0.01, "**", ifelse(p_v < 0.05, "*", "")))

      marker <- ifelse(j == best_f, " <-- iCAF", "")
      cat(sprintf("F%-7d | %10.3f %10.3f %10.4f%s | %10.3f %10.3f %10.4f%s%s\n",
        j, hr_t, z_t, p_t, sig_t, hr_v, z_v, p_v, sig_v, marker))
    }
  }
  cat("\n")

  # === C. LRT: Adding iCAF factor to PurIST + DeCAF baseline ===
  cat("--- C. LRT: Adding iCAF Factor to PurIST + DeCAF ---\n")
  icaf_fn <- paste0("Factor", best_f)

  base_t <- coxph(Surv(time, event) ~ PurIST + DeCAF, data = train_adj)
  full_t <- coxph(as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", icaf_fn)),
    data = train_adj)
  lrt_t <- anova(base_t, full_t)
  lrt_p_t <- lrt_t[["Pr(>|Chi|)"]][2]

  base_v <- coxph(Surv(time, event) ~ PurIST + DeCAF + strata(dataset), data = val_adj)
  full_v <- coxph(as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", icaf_fn,
    " + strata(dataset)")), data = val_adj)
  lrt_v <- anova(base_v, full_v)
  lrt_p_v <- lrt_v[["Pr(>|Chi|)"]][2]

  cat(sprintf("  Training LRT p  = %.4e (n=%d)\n", lrt_p_t, nrow(train_adj)))
  cat(sprintf("  Validation LRT p = %.4e (n=%d)\n", lrt_p_v, nrow(val_adj)))
  cat("\n")

  # === D. Median-split KM for iCAF factor ===
  cat("--- D. Median-Split KM for iCAF Factor ---\n")

  # Training KM
  tdf_km <- train_df[complete.cases(train_df[, c("time", "event", icaf_fn)]), ]
  tdf_km$icaf_group <- ifelse(tdf_km[[icaf_fn]] > median(tdf_km[[icaf_fn]]), "High", "Low")
  sf_t <- survfit(Surv(time, event) ~ icaf_group, data = tdf_km)
  med_t <- summary(sf_t)$table[, "median"]
  lr_t <- survdiff(Surv(time, event) ~ icaf_group, data = tdf_km)
  lr_p_t <- 1 - pchisq(lr_t$chisq, 1)
  cat(sprintf("  Training (n=%d):\n", nrow(tdf_km)))
  cat(sprintf("    High iCAF: median OS = %.1f months\n", med_t["icaf_group=High"]))
  cat(sprintf("    Low iCAF:  median OS = %.1f months\n", med_t["icaf_group=Low"]))
  cat(sprintf("    Log-rank p = %.4e\n", lr_p_t))

  # Validation KM
  val_df$icaf_group <- ifelse(val_df[[icaf_fn]] > median(val_df[[icaf_fn]]), "High", "Low")
  sf_v <- survfit(Surv(time, event) ~ icaf_group, data = val_df)
  med_v <- summary(sf_v)$table[, "median"]
  lr_v <- survdiff(Surv(time, event) ~ icaf_group, data = val_df)
  lr_p_v <- 1 - pchisq(lr_v$chisq, 1)
  cat(sprintf("  Validation (n=%d, pooled):\n", nrow(val_df)))
  cat(sprintf("    High iCAF: median OS = %.1f months\n", med_v["icaf_group=High"]))
  cat(sprintf("    Low iCAF:  median OS = %.1f months\n", med_v["icaf_group=Low"]))
  cat(sprintf("    Log-rank p = %.4e\n\n", lr_p_v))
}

# ============================================================
# SUMMARY TABLES
# ============================================================
cat("\n\n==================== UNADJUSTED SUMMARY ====================\n\n")
cat(sprintf("%-10s %3s %5s %6s | %8s %10s | %8s %10s | %10s %10s %10s\n",
  "Fit", "K", "Alpha", "iCAF_F",
  "Tr HR", "Tr p(unadj)",
  "Val HR", "Val p(unadj)",
  "Val KM Hi", "Val KM Lo", "Val KM p"))
cat(paste(rep("-", 115), collapse=""), "\n")

for (nm in names(fits)) {
  fit <- fits[[nm]]
  W <- fit$fit$W
  k <- fit$k
  alpha <- fit$alpha

  genes_t <- intersect(rownames(W), rownames(orig_W))
  genes_t <- intersect(genes_t, rownames(X_train))
  H_orig <- t(X_train[genes_t, ]) %*% orig_W[genes_t, ]
  H_new <- t(X_train[genes_t, ]) %*% W[genes_t, ]
  cors <- cor(H_orig[, 1], H_new, method = "spearman")
  best_f <- which.max(cors)
  icaf_fn <- paste0("Factor", best_f)

  train_df <- build_train_df(W, X_train, surv_train, all_train_sub)
  val_df <- build_val_df(W, val_data)

  cx_t <- coxph(as.formula(paste0("Surv(time, event) ~ ", icaf_fn)), data = train_df)
  s_t <- summary(cx_t)
  hr_t <- s_t$conf.int[1, 1]
  p_t <- s_t$coefficients[1, "Pr(>|z|)"]

  cx_v <- coxph(as.formula(paste0("Surv(time, event) ~ ", icaf_fn, " + strata(dataset)")),
    data = val_df)
  s_v <- summary(cx_v)
  hr_v <- s_v$conf.int[1, 1]
  p_v <- s_v$coefficients[1, "Pr(>|z|)"]

  val_df$icaf_group <- ifelse(val_df[[icaf_fn]] > median(val_df[[icaf_fn]]), "High", "Low")
  sf_v <- survfit(Surv(time, event) ~ icaf_group, data = val_df)
  med_v <- summary(sf_v)$table[, "median"]
  lr_v <- survdiff(Surv(time, event) ~ icaf_group, data = val_df)
  lr_p_v <- 1 - pchisq(lr_v$chisq, 1)

  cat(sprintf("%-10s %3d %5.2f %6d | %8.3f %11.4e | %8.3f %12.4e | %10.1f %10.1f %10.4e\n",
    nm, k, alpha, best_f,
    hr_t, p_t, hr_v, p_v,
    med_v["icaf_group=High"], med_v["icaf_group=Low"], lr_p_v))
}

cat("\n\n==================== ADJUSTED SUMMARY ====================\n\n")
cat(sprintf("%-10s %3s %5s %6s | %8s %10s %10s | %8s %10s %10s\n",
  "Fit", "K", "Alpha", "iCAF_F",
  "Tr HR", "Tr p(adj)", "Tr LRT",
  "Val HR", "Val p(adj)", "Val LRT"))
cat(paste(rep("-", 105), collapse=""), "\n")

for (nm in names(fits)) {
  fit <- fits[[nm]]
  W <- fit$fit$W
  k <- fit$k
  alpha <- fit$alpha

  genes_t <- intersect(rownames(W), rownames(orig_W))
  genes_t <- intersect(genes_t, rownames(X_train))
  H_orig <- t(X_train[genes_t, ]) %*% orig_W[genes_t, ]
  H_new <- t(X_train[genes_t, ]) %*% W[genes_t, ]
  cors <- cor(H_orig[, 1], H_new, method = "spearman")
  best_f <- which.max(cors)
  icaf_fn <- paste0("Factor", best_f)

  train_df <- build_train_df(W, X_train, surv_train, all_train_sub)
  val_df <- build_val_df(W, val_data)

  adj_cols <- c("time", "event", "PurIST", "DeCAF", paste0("Factor", seq_len(k)))
  train_adj <- train_df[complete.cases(train_df[, adj_cols]), ]
  val_adj <- val_df[complete.cases(val_df[, c("time", "event", "PurIST", "DeCAF")]), ]

  fml_t <- as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", icaf_fn))
  cx_t <- coxph(fml_t, data = train_adj)
  s_t <- summary(cx_t)
  row_t <- grep(paste0("^", icaf_fn, "$"), rownames(s_t$coefficients))
  hr_t <- exp(s_t$coefficients[row_t, "coef"])
  p_t_adj <- s_t$coefficients[row_t, "Pr(>|z|)"]

  base_t <- coxph(Surv(time, event) ~ PurIST + DeCAF, data = train_adj)
  lrt_t <- anova(base_t, cx_t)
  lrt_p_t <- lrt_t[["Pr(>|Chi|)"]][2]

  fml_v <- as.formula(paste0("Surv(time, event) ~ PurIST + DeCAF + ", icaf_fn, " + strata(dataset)"))
  cx_v <- coxph(fml_v, data = val_adj)
  s_v <- summary(cx_v)
  row_v <- grep(paste0("^", icaf_fn, "$"), rownames(s_v$coefficients))
  hr_v <- exp(s_v$coefficients[row_v, "coef"])
  p_v_adj <- s_v$coefficients[row_v, "Pr(>|z|)"]

  base_v <- coxph(Surv(time, event) ~ PurIST + DeCAF + strata(dataset), data = val_adj)
  lrt_v <- anova(base_v, cx_v)
  lrt_p_v <- lrt_v[["Pr(>|Chi|)"]][2]

  cat(sprintf("%-10s %3d %5.2f %6d | %8.3f %10.4e %10.4e | %8.3f %10.4e %10.4e\n",
    nm, k, alpha, best_f,
    hr_t, p_t_adj, lrt_p_t,
    hr_v, p_v_adj, lrt_p_v))
}
cat("\n")
