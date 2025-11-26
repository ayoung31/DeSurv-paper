# Utility functions for simulation metrics ---------------------------------

.null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x)) {
    y
  } else {
    x
  }
}

.shared_gene_matrix <- function(true_mat, est_mat) {
  genes <- intersect(rownames(true_mat), rownames(est_mat))
  if (!length(genes)) {
    stop("No overlapping genes between matrices.")
  }
  list(
    truth = true_mat[genes, , drop = FALSE],
    est = est_mat[genes, , drop = FALSE]
  )
}

.match_factors_by_correlation <- function(true_mat, est_mat) {
  mats <- .shared_gene_matrix(true_mat, est_mat)
  corr_mat <- stats::cor(mats$truth, mats$est, use = "pairwise.complete.obs")
  if (!is.matrix(corr_mat)) {
    corr_mat <- matrix(corr_mat, nrow = 1)
  }
  n_true <- nrow(corr_mat)
  n_est <- ncol(corr_mat)
  used_true <- rep(FALSE, n_true)
  used_est <- rep(FALSE, n_est)
  matches <- list()
  while (TRUE) {
    remaining <- corr_mat
    remaining[used_true, ] <- NA_real_
    remaining[, used_est] <- NA_real_
    idx <- which.max(abs(remaining))
    if (!length(idx) || is.na(remaining[idx])) break
    rc <- arrayInd(idx, dim(remaining))
    matches[[length(matches) + 1]] <- list(
      true_idx = rc[1],
      est_idx = rc[2],
      correlation = remaining[idx]
    )
    used_true[rc[1]] <- TRUE
    used_est[rc[2]] <- TRUE
    if (all(used_true) || all(used_est)) break
  }
  tibble::tibble(
    true_factor = vapply(matches, `[[`, integer(1), "true_idx"),
    est_factor = vapply(matches, `[[`, integer(1), "est_idx"),
    correlation = vapply(matches, `[[`, numeric(1), "correlation")
  )
}

.principal_angle_summary <- function(A, B) {
  if (ncol(A) == 0L || ncol(B) == 0L) {
    return(list(angles = numeric(), mean = NA_real_, max = NA_real_))
  }
  Qa <- qr.Q(qr(A))
  Qb <- qr.Q(qr(B))
  M <- t(Qa) %*% Qb
  sv <- La.svd(M)
  vals <- pmin(pmax(sv$d, -1), 1)
  angles <- acos(vals)
  list(
    angles = angles,
    mean = mean(angles),
    max = max(angles)
  )
}

.project_expression_onto_W <- function(X, W) {
  XtW <- crossprod(W)
  RHS <- crossprod(W, X)
  sol <- tryCatch(
    qr.solve(XtW, RHS),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Need MASS package to use generalized inverse.")
      }
      MASS::ginv(XtW) %*% RHS
    }
  )
  sol
}

.predict_simulation_risk <- function(fit, truth) {
  stopifnot(!is.null(truth$counts), !is.null(truth$surv))
  samp <- truth$surv
  samp$ID <- samp$patient
  samp$dataset <- .null_coalesce(truth$config_name, "simulation")
  rownames(samp) <- samp$ID
  dataset <- list(
    ex = truth$counts,
    sampInfo = samp,
    dataname = samp$dataset[1]
  )
  preds <- desurv_predict_dataset(fit, dataset)
  dplyr::left_join(
    preds,
    truth$surv,
    by = c("sample_id" = "patient")
  )
}

.marker_sets_from_W <- function(W, top_n = 200) {
  apply(W, 2, function(col) {
    genes <- names(sort(col, decreasing = TRUE))
    genes[seq_len(min(top_n, length(genes)))]
  })
}

# Column-wise correlation metric -------------------------------------------

metric_column_correlations <- function(truth, fit) {
  stopifnot(!is.null(truth$W_true), !is.null(fit$W))
  tibble::as_tibble(.match_factors_by_correlation(truth$W_true, fit$W))
}

# Survival subspace distance -----------------------------------------------

metric_survival_subspace_distance <- function(truth, fit, beta_threshold = 1e-3) {
  stopifnot(!is.null(truth$W_true), !is.null(truth$beta_true), !is.null(fit$W), !is.null(fit$beta))
  mats <- .shared_gene_matrix(truth$W_true, fit$W)
  true_idx <- which(abs(truth$beta_true) > beta_threshold)
  est_idx <- which(abs(fit$beta) > beta_threshold)
  summary <- .principal_angle_summary(
    mats$truth[, true_idx, drop = FALSE],
    mats$est[, est_idx, drop = FALSE]
  )
  tibble::tibble(
    mean_angle = summary$mean,
    max_angle = summary$max,
    n_true = length(true_idx),
    n_est = length(est_idx)
  )
}

# Reconstruction error -----------------------------------------------------

metric_reconstruction_error <- function(truth, fit) {
  stopifnot(!is.null(truth$counts), !is.null(fit$W))
  mats <- .shared_gene_matrix(truth$counts, fit$W)
  H_hat <- .project_expression_onto_W(mats$truth, mats$est)
  recon <- mats$est %*% H_hat
  num <- base::norm(mats$truth - recon, type = "F")
  den <- base::norm(mats$truth, type = "F")
  tibble::tibble(
    reconstruction_error = num / den
  )
}

# Factor detection metrics -------------------------------------------------

metric_factor_detection <- function(truth, fit, beta_threshold = 0.1) {
  stopifnot(!is.null(truth$beta_true), !is.null(fit$beta), !is.null(truth$W_true), !is.null(fit$W))
  matches <- metric_column_correlations(truth, fit)
  true_active <- which(abs(truth$beta_true) > beta_threshold)
  est_active <- which(abs(fit$beta) > beta_threshold)
  detected <- matches$est_factor[matches$true_factor %in% true_active]
  tp <- sum(est_active %in% detected)
  fn <- length(true_active) - tp
  fp <- sum(est_active %in% setdiff(matches$est_factor, detected))
  tn <- (ncol(fit$W) - length(est_active)) - fn
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  fdr <- if ((tp + fp) > 0) fp / (tp + fp) else NA_real_
  tibble::tibble(
    precision = precision,
    recall = recall,
    specificity = specificity,
    fdr = fdr,
    true_active = length(true_active),
    est_active = length(est_active)
  )
}

metric_effect_size_correlation <- function(truth, fit) {
  stopifnot(!is.null(truth$beta_true), !is.null(fit$beta), !is.null(truth$W_true), !is.null(fit$W))
  matches <- metric_column_correlations(truth, fit)
  est_beta <- fit$beta[matches$est_factor]
  true_beta <- truth$beta_true[matches$true_factor]
  corr <- if (length(est_beta) > 1) stats::cor(true_beta, est_beta) else NA_real_
  tibble::tibble(
    effect_size_correlation = corr,
    n_matched = length(est_beta)
  )
}

metric_factor_roc_curve <- function(truth, fit, beta_threshold = 0.1) {
  stopifnot(!is.null(truth$beta_true), !is.null(fit$beta), !is.null(truth$W_true), !is.null(fit$W))
  matches <- metric_column_correlations(truth, fit)
  scores <- abs(fit$beta[matches$est_factor])
  labels <- as.integer(abs(truth$beta_true[matches$true_factor]) > beta_threshold)
  ord <- order(scores, decreasing = TRUE)
  scores <- scores[ord]
  labels <- labels[ord]
  tp <- cumsum(labels)
  fp <- cumsum(1 - labels)
  total_pos <- sum(labels)
  total_neg <- length(labels) - total_pos
  tpr <- ifelse(total_pos > 0, tp / total_pos, NA_real_)
  fpr <- ifelse(total_neg > 0, fp / total_neg, NA_real_)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  recall <- tpr
  tibble::tibble(
    threshold = scores,
    tpr = tpr,
    fpr = fpr,
    precision = precision,
    recall = recall
  )
}

# Patient-level metrics ----------------------------------------------------

metric_test_cindex <- function(truth, fit) {
  preds <- .predict_simulation_risk(fit, truth)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for C-index.")
  }
  surv_obj <- survival::Surv(preds$time, preds$status)
  cc <- survival::survConcordance(surv_obj ~ preds$risk_score)
  tibble::tibble(c_index = cc$concordance)
}

metric_time_dependent_auc <- function(truth, fit, times) {
  preds <- .predict_simulation_risk(fit, truth)
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package `timeROC` is required for time-dependent AUC.")
  }
  roc <- timeROC::timeROC(
    T = preds$time,
    delta = preds$status,
    marker = preds$risk_score,
    cause = 1,
    times = times,
    iid = FALSE
  )
  tibble::tibble(
    time = roc$times,
    auc = roc$AUC
  )
}

metric_integrated_brier_score <- function(truth, fit, times) {
  preds <- .predict_simulation_risk(fit, truth)
  if (!requireNamespace("pec", quietly = TRUE)) {
    stop("Package `pec` is required for Brier score.")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for Brier score.")
  }
  data <- data.frame(time = preds$time, status = preds$status, risk = preds$risk_score)
  cox_fit <- survival::coxph(survival::Surv(time, status) ~ risk, data = data, x = TRUE, y = TRUE)
  pec_obj <- pec::pec(
    object = list(desurv = cox_fit),
    formula = survival::Surv(time, status) ~ 1,
    data = data,
    times = times
  )
  ibs <- pec::crps(pec_obj, times = times)[, "desurv"]
  tibble::tibble(time = times, brier = ibs)
}

metric_km_separation <- function(truth, fit, quantile = 0.5) {
  preds <- .predict_simulation_risk(fit, truth)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for KM curves.")
  }
  cutoff <- stats::quantile(preds$risk_score, probs = quantile)
  group <- ifelse(preds$risk_score >= cutoff, "high", "low")
  surv_obj <- survival::Surv(preds$time, preds$status)
  km <- survival::survfit(surv_obj ~ group)
  tidy <- broom::tidy(km)
  logrank <- survival::survdiff(surv_obj ~ group)
  pval <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  tibble::tibble(
    time = tidy$time,
    surv = tidy$estimate,
    std_err = tidy$std.error,
    conf_low = tidy$conf.low,
    conf_high = tidy$conf.high,
    strata = tidy$strata,
    logrank_p = pval
  )
}

metric_calibration_curve <- function(truth, fit, horizon, bins = 5) {
  preds <- .predict_simulation_risk(fit, truth)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for calibration.")
  }
  data <- data.frame(time = preds$time, status = preds$status, risk = preds$risk_score)
  cox_fit <- survival::coxph(survival::Surv(time, status) ~ risk, data = data, x = TRUE, y = TRUE)
  basehaz <- survival::basehaz(cox_fit, centered = FALSE)
  haz <- stats::approx(basehaz$time, basehaz$hazard, xout = horizon, rule = 2)$y
  surv_prob <- exp(-haz * exp(data$risk))
  cuts <- stats::quantile(surv_prob, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
  groups <- cut(surv_prob, breaks = cuts, include.lowest = TRUE)
  observed <- tapply(
    X = as.numeric(data$time > horizon),
    INDEX = groups,
    FUN = mean,
    na.rm = TRUE
  )
  tibble::tibble(
    group = levels(groups),
    predicted = tapply(surv_prob, groups, mean, na.rm = TRUE),
    observed = observed
  )
}

# Stability metrics --------------------------------------------------------

metric_subspace_stability <- function(fit_list, beta_threshold = 0.1) {
  stopifnot(length(fit_list) >= 2)
  combos <- utils::combn(seq_along(fit_list), 2)
  vals <- apply(combos, 2, function(idx) {
    fit_a <- fit_list[[idx[1]]]
    fit_b <- fit_list[[idx[2]]]
    mats <- .shared_gene_matrix(fit_a$W, fit_b$W)
    idx_a <- which(abs(fit_a$beta) > beta_threshold)
    idx_b <- which(abs(fit_b$beta) > beta_threshold)
    summary <- .principal_angle_summary(
      mats$truth[, idx_a, drop = FALSE],
      mats$est[, idx_b, drop = FALSE]
    )
    summary$mean
  })
  tibble::tibble(
    mean_subspace_angle = mean(vals),
    max_subspace_angle = max(vals),
    comparisons = ncol(combos)
  )
}

metric_risk_score_stability <- function(fit_list, truth) {
  stopifnot(length(fit_list) >= 2)
  preds <- lapply(fit_list, function(fit) .predict_simulation_risk(fit, truth)$risk_score)
  mat <- do.call(cbind, preds)
  corr <- stats::cor(mat, use = "pairwise.complete.obs")
  tibble::tibble(
    mean_correlation = mean(corr[upper.tri(corr)]),
    min_correlation = min(corr[upper.tri(corr)])
  )
}

# Gene-set recovery --------------------------------------------------------

metric_gene_set_recovery <- function(truth, fit, top_n = 200, beta_threshold = 0.1) {
  truth_sets <- .marker_sets_from_W(truth$W_true, top_n = top_n)
  est_sets <- .marker_sets_from_W(fit$W, top_n = top_n)
  matches <- metric_column_correlations(truth, fit)
  active <- which(abs(truth$beta_true) > beta_threshold)
  rows <- lapply(matches$true_factor, function(idx) {
    est_idx <- matches$est_factor[matches$true_factor == idx]
    true_set <- truth_sets[[idx]]
    est_set <- est_sets[[est_idx]]
    tp <- length(intersect(true_set, est_set))
    precision <- tp / length(est_set)
    recall <- tp / length(true_set)
    f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
    list(
      true_factor = idx,
      est_factor = est_idx,
      precision = precision,
      recall = recall,
      f1 = f1,
      survival_associated = idx %in% active
    )
  })
  dplyr::bind_rows(rows)
}

metric_global_gene_set_recovery <- function(truth, fit, top_n = 200, beta_threshold = 0.1) {
  truth_sets <- .marker_sets_from_W(truth$W_true, top_n = top_n)
  active <- which(abs(truth$beta_true) > beta_threshold)
  truth_union <- unique(unlist(truth_sets[active]))
  est_sets <- .marker_sets_from_W(fit$W, top_n = top_n)
  est_active <- which(abs(fit$beta) > beta_threshold)
  est_union <- unique(unlist(est_sets[est_active]))
  tp <- length(intersect(truth_union, est_union))
  precision <- tp / length(est_union)
  recall <- tp / length(truth_union)
  f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
  jaccard <- tp / length(union(truth_union, est_union))
  tibble::tibble(
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = jaccard
  )
}

metric_gene_set_overlap_matrix <- function(truth, fit, top_n = 200) {
  truth_sets <- .marker_sets_from_W(truth$W_true, top_n = top_n)
  est_sets <- .marker_sets_from_W(fit$W, top_n = top_n)
  mat <- outer(
    seq_along(truth_sets),
    seq_along(est_sets),
    Vectorize(function(i, j) {
      inter <- length(intersect(truth_sets[[i]], est_sets[[j]]))
      union_size <- length(union(truth_sets[[i]], est_sets[[j]]))
      if (union_size == 0) 0 else inter / union_size
    })
  )
  rownames(mat) <- names(truth_sets)
  colnames(mat) <- names(est_sets)
  mat
}

# BO diagnostics -----------------------------------------------------------

metric_bo_hyperparameter_distribution <- function(history_df) {
  stopifnot(all(c("k", "alpha", "lambdaH") %in% names(history_df)))
  history_df %>%
    dplyr::summarise(
      k_median = stats::median(k, na.rm = TRUE),
      alpha_median = stats::median(alpha, na.rm = TRUE),
      lambdaH_median = stats::median(lambdaH, na.rm = TRUE),
      k_iqr = stats::IQR(k, na.rm = TRUE),
      alpha_iqr = stats::IQR(alpha, na.rm = TRUE),
      lambdaH_iqr = stats::IQR(lambdaH, na.rm = TRUE)
    )
}

metric_performance_vs_oracle <- function(summary_df, k_true) {
  stopifnot(all(c("k", "mean_cindex") %in% names(summary_df)))
  best_idx <- which.max(summary_df$mean_cindex)
  best_row <- summary_df[best_idx, , drop = FALSE]
  tibble::tibble(
    k_true = k_true,
    k_selected = best_row$k,
    cindex_selected = best_row$mean_cindex
  )
}

metric_bo_convergence <- function(history_df) {
  stopifnot(all(c("iteration", "mean_cindex") %in% names(history_df)))
  history_df %>%
    dplyr::arrange(iteration) %>%
    dplyr::mutate(best_so_far = cummax(mean_cindex))
}
