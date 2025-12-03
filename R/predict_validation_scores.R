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
  if (ncol(Z) == 0L) {
    warning("No factors available for dataset: ", dataset_name)
    return(tibble::tibble())
  }

  beta <- as.numeric(fit$beta)
  risk <- drop(Z %*% beta)

  compute_train_stats <- function(fit) {
    if (is.null(fit$data) || is.null(fit$data$X)) {
      return(NULL)
    }
    train_genes <- intersect(rownames(fit$W), rownames(fit$data$X))
    if (!length(train_genes)) {
      return(NULL)
    }
    X_train <- fit$data$X[train_genes, , drop = FALSE]
    W_train <- fit$W[train_genes, , drop = FALSE]
    Z_train <- t(X_train) %*% W_train
    list(
      mean = colMeans(Z_train),
      sd = apply(Z_train, 2, stats::sd)
    )
  }

  stats <- compute_train_stats(fit)
  if (!is.null(stats)) {
    sd_vec <- stats$sd
    sd_vec[sd_vec < 1e-12] <- 1
    Z_scaled <- sweep(Z, 2, stats$mean, FUN = "-")
    Z_scaled <- sweep(Z_scaled, 2, sd_vec, FUN = "/")
  } else {
    Z_scaled <- scale(Z)
    center_attr <- attr(Z_scaled, "scaled:center")
    scale_attr <- attr(Z_scaled, "scaled:scale")
    attr(Z_scaled, "scaled:center") <- NULL
    attr(Z_scaled, "scaled:scale") <- NULL
    if (any(!is.finite(Z_scaled))) {
      Z_scaled <- Z
    }
  }

  score_df <- as.data.frame(Z_scaled)
  colnames(score_df) <- paste0("X", seq_len(ncol(Z_scaled)))
  score_df$sample_id <- rownames(Z_scaled)
  score_df$risk_score <- risk
  score_df$dataset <- dataset_name

  tibble::as_tibble(score_df, .name_repair = "minimal")
}

desurv_predict_validation <- function(fit, data_list) {
  preds <- lapply(data_list, function(dataset) desurv_predict_dataset(fit, dataset))
  dplyr::bind_rows(preds)
}
