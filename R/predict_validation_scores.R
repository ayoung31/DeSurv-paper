desurv_predict_dataset <- function(fit, dataset) {
  stopifnot(!is.null(fit$W), !is.null(dataset$ex))
  dataset_name <- dataset$dataname
  if (is.null(dataset_name) || !nzchar(dataset_name)) {
    dataset_name <- unique(dataset$sampInfo$dataset)
    dataset_name <- dataset_name[!is.na(dataset_name)]
    dataset_name <- if (length(dataset_name)) dataset_name[[1]] else "validation"
  }

  genes <- intersect(rownames(fit$W), rownames(dataset$ex))
  if (!length(genes)) {
    warning("No overlapping genes between fit and dataset: ", dataset_name)
    return(tibble::tibble())
  }

  W_sub <- fit$W[genes, , drop = FALSE]
  X_sub <- dataset$ex[genes, , drop = FALSE]

  Z <- t(X_sub) %*% W_sub
  ns <- which(fit$sdZ > 1e-12)
  if (!length(ns)) {
    warning("No non-zero loadings found for dataset: ", dataset_name)
    return(tibble::tibble())
  }

  Z_use <- Z[, ns, drop = FALSE]
  meanZ <- fit$meanZ[, ns, drop = FALSE]
  sdZ <- fit$sdZ[, ns, drop = FALSE]

  Z_centered <- sweep(Z_use, 2, meanZ, FUN = "-")
  Z_scaled <- sweep(Z_centered, 2, sdZ, FUN = "/")

  beta <- fit$beta[ns]
  risk <- as.numeric(Z_centered %*% (beta * sdZ))

  score_df <- as.data.frame(Z_scaled)
  score_df$sample_id <- rownames(Z_scaled)
  score_df$risk_score <- risk
  score_df$dataset <- dataset_name

  tibble::as_tibble(score_df, .name_repair = "minimal")
}

desurv_predict_validation <- function(fit, data_list) {
  preds <- lapply(data_list, function(dataset) desurv_predict_dataset(fit, dataset))
  dplyr::bind_rows(preds)
}
