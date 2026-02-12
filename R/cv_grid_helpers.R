# Helper functions for CV grid search over k x alpha x ntop
#
# This module supports _targets_cv_grid.R for exhaustive cross-validation
# over a grid of k (2-12), alpha (0 to 0.95 by 0.05), and ntop (NULL, 300)
# with fixed lambda=0.3, nu=0.05, lambdaW=0, lambdaH=0.

#' Create a grid of (k, alpha, ntop) parameter combinations
#'
#' @param k_values Integer vector of k values to test
#' @param alpha_values Numeric vector of alpha values to test
#' @param ntop_values List of ntop values to test (NULL means all genes)
#' @return List of lists, each with $k, $alpha, and $ntop elements
create_cv_grid <- function(k_values = 2:12,
                           alpha_values = seq(0, 0.95, by = 0.05),
                           ntop_values = list(NULL)) {
  # Build base grid of k x alpha
  base_grid <- expand.grid(
    k = as.integer(k_values),
    alpha = alpha_values,
    ntop_idx = seq_along(ntop_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Convert to list of lists for targets branching
  lapply(seq_len(nrow(base_grid)), function(i) {
    list(
      k = base_grid$k[i],
      alpha = base_grid$alpha[i],
      ntop = ntop_values[[base_grid$ntop_idx[i]]]
    )
  })
}

#' Run cross-validation for a single (k, alpha) point
#'
#' @param data List with $ex (expression matrix), $sampInfo (data frame with
#'   time, event, dataset columns)
#' @param k Integer number of latent factors
#' @param alpha Numeric supervision parameter in [0, 1]
#' @param fixed_params List of fixed parameters (lambda, nu, lambdaW, lambdaH, ntop)
#' @param nfolds Integer number of CV folds
#' @param n_starts Integer number of random initializations
#' @param seed Integer random seed
#' @param verbose Logical; print progress messages
#' @return List with k, alpha, mean_cindex, se_cindex, cv_results
run_cv_grid_point <- function(data,
                              k,
                              alpha,
                              fixed_params = list(
                                lambda = 0.3,
                                nu = 0.05,
                                lambdaW = 0,
                                lambdaH = 0,
                                ntop = NULL
                              ),
                              nfolds = 5,
                              n_starts = 30,
                              seed = 123,
                              verbose = TRUE) {

  if (verbose) {
    ntop_str <- if (is.null(fixed_params$ntop)) "ALL" else as.character(fixed_params$ntop)
    message(sprintf("Running CV for k=%d, alpha=%.2f, ntop=%s", k, alpha, ntop_str))
  }

  # Extract fixed parameters with defaults
  lambda <- fixed_params$lambda %||% 0.3
  nu <- fixed_params$nu %||% 0.05
  lambdaW <- fixed_params$lambdaW %||% 0
  lambdaH <- fixed_params$lambdaH %||% 0
  ntop <- fixed_params$ntop  # NULL means use all genes for cindex

  # Run cross-validation using DeSurv::desurv_cv with cv_only=TRUE
  # This avoids the final refit and just returns CV results
  cv_result <- tryCatch({
    DeSurv::desurv_cv(
      X = data$ex,
      y = data$sampInfo$time,
      d = data$sampInfo$event,
      dataset = data$sampInfo$dataset,
      k_grid = k,
      alpha_grid = alpha,
      lambda_grid = lambda,
      nu_grid = nu,
      lambdaW_grid = lambdaW,
      lambdaH_grid = lambdaH,
      n_starts = n_starts,
      nfolds = nfolds,
      seed = seed,
      tol = 1e-5,
      maxit = 3000,
      engine = "warmstart",
      rule = "max",
      parallel_grid = TRUE,
      ncores_grid = n_starts,
      ntop = ntop,
      verbose = verbose,
      preprocess = FALSE,
      cv_only = TRUE
    )
  }, error = function(e) {
    warning(sprintf(
      "CV failed for k=%d, alpha=%.2f: %s",
      k, alpha, conditionMessage(e)
    ))
    return(NULL)
  })

  # Extract summary statistics
  if (is.null(cv_result) || is.null(cv_result$summary) || nrow(cv_result$summary) == 0) {
    return(list(
      k = k,
      alpha = alpha,
      ntop = ntop,
      mean_cindex = NA_real_,
      se_cindex = NA_real_,
      cv_results = NULL
    ))
  }

  # The summary should have exactly one row for our single (k, alpha) point
  summ <- cv_result$summary

  list(
    k = k,
    alpha = alpha,
    ntop = ntop,
    mean_cindex = summ$mean_cindex[1],
    se_cindex = summ$se_cindex[1],
    cv_results = cv_result$cv_results
  )
}

#' Aggregate results from all CV grid points
#'
#' @param result_list List of results from run_cv_grid_point calls
#' @return Tibble with columns: k, alpha, mean_cindex, se_cindex
aggregate_cv_grid_results <- function(result_list) {
  # Filter out NULL results
  valid_results <- purrr::compact(result_list)

  if (length(valid_results) == 0) {
    warning("No valid CV results to aggregate")
    return(tibble::tibble(
      k = integer(),
      alpha = numeric(),
      ntop = numeric(),
      mean_cindex = numeric(),
      se_cindex = numeric()
    ))
  }

  tibble::tibble(
    k = purrr::map_int(valid_results, ~ as.integer(.x$k)),
    alpha = purrr::map_dbl(valid_results, ~ as.numeric(.x$alpha)),
    ntop = purrr::map_dbl(valid_results, ~ if (is.null(.x$ntop)) NA_real_ else as.numeric(.x$ntop)),
    mean_cindex = purrr::map_dbl(valid_results, ~ .x$mean_cindex %||% NA_real_),
    se_cindex = purrr::map_dbl(valid_results, ~ .x$se_cindex %||% NA_real_)
  )
}

#' Fit DeSurv on full training data for a single (k, alpha, ntop) grid point
#'
#' @param data List with $ex, $sampInfo (time, event, dataset columns)
#' @param k Integer number of latent factors
#' @param alpha Numeric supervision parameter
#' @param fixed_params List of fixed parameters (lambda, nu, lambdaW, lambdaH, ntop)
#' @param n_starts Integer number of random initializations (default 100)
#' @param seed Integer random seed
#' @param tol Numeric convergence tolerance
#' @param maxit Integer max iterations
#' @param verbose Logical
#' @return List with k, alpha, ntop, fit (desurv_fit object or NULL), cindex, convergence
fit_grid_point <- function(data, k, alpha,
                           fixed_params = list(lambda = 0.3, nu = 0.05,
                                               lambdaW = 0, lambdaH = 0,
                                               ntop = NULL),
                           n_starts = 100,
                           seed = 123,
                           tol = 1e-5,
                           maxit = 3000,
                           verbose = TRUE) {

  ntop <- fixed_params$ntop
  ntop_str <- if (is.null(ntop)) "ALL" else as.character(ntop)

  if (verbose) {
    message(sprintf("Fitting full model for k=%d, alpha=%.2f, ntop=%s, ninit=%d",
                    k, alpha, ntop_str, n_starts))
  }

  lambda <- fixed_params$lambda %||% 0.3
  nu <- fixed_params$nu %||% 0.05
  lambdaW <- fixed_params$lambdaW %||% 0
  lambdaH <- fixed_params$lambdaH %||% 0

  fit <- tryCatch({
    DeSurv::desurv_fit(
      X = data$ex,
      y = data$sampInfo$time,
      d = data$sampInfo$event,
      k = k,
      alpha = alpha,
      lambda = lambda,
      nu = nu,
      lambdaW = lambdaW,
      lambdaH = lambdaH,
      ninit = n_starts,
      parallel_init = TRUE,
      ncores_init = n_starts,
      seed = seed,
      tol = tol,
      maxit = maxit,
      verbose = verbose
    )
  }, error = function(e) {
    warning(sprintf(
      "Full fit failed for k=%d, alpha=%.2f, ntop=%s: %s",
      k, alpha, ntop_str, conditionMessage(e)
    ))
    return(NULL)
  })

  list(
    k = k,
    alpha = alpha,
    ntop = ntop,
    fit = fit,
    cindex = if (!is.null(fit)) fit$cindex else NA_real_,
    convergence = if (!is.null(fit)) fit$convergence else NA
  )
}

#' Aggregate results from all full-data fit grid points
#'
#' @param result_list List of results from fit_grid_point calls
#' @return Tibble with columns: k, alpha, ntop, cindex, convergence
aggregate_cv_grid_fit_results <- function(result_list) {
  valid_results <- purrr::compact(result_list)

  if (length(valid_results) == 0) {
    warning("No valid fit results to aggregate")
    return(tibble::tibble(
      k = integer(),
      alpha = numeric(),
      ntop = numeric(),
      cindex = numeric(),
      convergence = logical()
    ))
  }

  tibble::tibble(
    k = purrr::map_int(valid_results, ~ as.integer(.x$k)),
    alpha = purrr::map_dbl(valid_results, ~ as.numeric(.x$alpha)),
    ntop = purrr::map_dbl(valid_results, ~ if (is.null(.x$ntop)) NA_real_ else as.numeric(.x$ntop)),
    cindex = purrr::map_dbl(valid_results, ~ .x$cindex %||% NA_real_),
    convergence = purrr::map_lgl(valid_results, ~ .x$convergence %||% NA)
  )
}

#' Validate a single grid point against external datasets
#'
#' For each validation dataset, computes concordance using two methods:
#'   - transfer_beta: project via W and beta (or top-gene theta) from training
#'   - refit_cox: project via W, then refit Cox on validation Z scores
#'
#' Also computes pooled results across all datasets.
#'
#' @param grid_fit List from fit_grid_point() with $fit, $k, $alpha, $ntop
#' @param val_datasets Named list of preprocessed validation datasets, each with
#'   $ex (genes x samples), $sampInfo (data.frame with time, event columns)
#' @return Tibble with columns: k, alpha, ntop, dataset, method, cindex,
#'   cindex_se, n_samples, n_events
validate_grid_point <- function(grid_fit, val_datasets) {
  empty_result <- tibble::tibble(
    k = integer(), alpha = numeric(), ntop = numeric(),
    dataset = character(), method = character(),
    cindex = numeric(), cindex_se = numeric(),
    n_samples = integer(), n_events = integer()
  )

  if (is.null(grid_fit$fit)) {
    return(empty_result)
  }

  fit <- grid_fit$fit
  k_val <- as.integer(grid_fit$k)
  alpha_val <- as.numeric(grid_fit$alpha)
  ntop_val <- if (is.null(grid_fit$ntop)) NA_real_ else as.numeric(grid_fit$ntop)
  ntop_raw <- grid_fit$ntop  # keep NULL for logic below

  W <- fit$W
  beta <- fit$beta
  train_genes <- rownames(W)

  # Per-dataset results
  per_ds_rows <- list()
  # Collectors for pooling
  pooled_lp <- numeric(0)
  pooled_Z <- list()
  pooled_time <- numeric(0)
  pooled_event <- integer(0)
  pooled_ds_label <- character(0)

  for (ds_name in names(val_datasets)) {
    val_ds <- val_datasets[[ds_name]]
    val_genes <- rownames(val_ds$ex)
    common_genes <- intersect(train_genes, val_genes)

    if (length(common_genes) < 2) {
      per_ds_rows[[ds_name]] <- tibble::tibble(
        k = rep(k_val, 2), alpha = rep(alpha_val, 2),
        ntop = rep(ntop_val, 2), dataset = rep(ds_name, 2),
        method = c("transfer_beta", "refit_cox"),
        cindex = rep(NA_real_, 2), cindex_se = rep(NA_real_, 2),
        n_samples = rep(0L, 2), n_events = rep(0L, 2)
      )
      next
    }

    # Subset to common genes
    X_val <- val_ds$ex[common_genes, , drop = FALSE]
    W_common <- W[common_genes, , drop = FALSE]
    si <- val_ds$sampInfo
    time_val <- si$time
    event_val <- as.integer(si$event)

    # Filter valid samples
    valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
    if (length(valid_idx) < 2) {
      per_ds_rows[[ds_name]] <- tibble::tibble(
        k = rep(k_val, 2), alpha = rep(alpha_val, 2),
        ntop = rep(ntop_val, 2), dataset = rep(ds_name, 2),
        method = c("transfer_beta", "refit_cox"),
        cindex = rep(NA_real_, 2), cindex_se = rep(NA_real_, 2),
        n_samples = rep(length(valid_idx), 2),
        n_events = rep(sum(event_val[valid_idx]), 2)
      )
      next
    }
    X_val <- X_val[, valid_idx, drop = FALSE]
    time_val <- time_val[valid_idx]
    event_val <- event_val[valid_idx]
    n_samp <- length(valid_idx)
    n_evt <- sum(event_val)

    # Compute Z = t(X) %*% W (samples x k)
    Z_val <- t(X_val) %*% W_common

    # Compute LP using the same logic as .desurv_lp_with_top_genes
    lp_val <- tryCatch({
      if (is.null(ntop_raw)) {
        # No top-gene subsetting: lp = Z %*% beta
        drop(Z_val %*% beta)
      } else {
        # Top-gene subsetting with theta normalization
        top_info <- DeSurv::desurv_get_top_genes(W_common, as.integer(ntop_raw))
        idx_mat <- top_info$top_indices
        idx <- unique(as.integer(unlist(idx_mat, use.names = FALSE)))
        idx <- idx[!is.na(idx) & idx >= 1L & idx <= nrow(W_common)]
        if (!length(idx)) {
          drop(Z_val %*% beta)
        } else {
          theta <- W_common %*% beta
          theta_sub <- theta[idx, , drop = FALSE]
          theta_norm <- sqrt(sum(theta_sub^2))
          if (!is.finite(theta_norm)) {
            theta_sub[] <- 0
          } else if (theta_norm > 0) {
            theta_sub <- theta_sub / theta_norm
          }
          X_sub <- X_val[idx, , drop = FALSE]
          drop(t(X_sub) %*% theta_sub)
        }
      }
    }, error = function(e) rep(NA_real_, n_samp))

    # Method 1: transfer_beta
    m1_cindex <- NA_real_
    m1_se <- NA_real_
    tryCatch({
      if (!anyNA(lp_val) && !any(is.infinite(lp_val)) && n_evt > 0) {
        cc <- survival::concordance(
          survival::Surv(time_val, event_val) ~ lp_val,
          reverse = TRUE
        )
        m1_cindex <- cc$concordance
        m1_se <- sqrt(cc$var)
      }
    }, error = function(e) NULL)

    # Method 2: refit_cox
    m2_cindex <- NA_real_
    m2_se <- NA_real_
    tryCatch({
      if (n_evt > 0) {
        cox_df <- as.data.frame(Z_val)
        colnames(cox_df) <- paste0("Z", seq_len(ncol(Z_val)))
        cox_df$time <- time_val
        cox_df$event <- event_val
        z_terms <- paste(colnames(cox_df)[seq_len(ncol(Z_val))], collapse = " + ")
        cox_fit <- survival::coxph(
          as.formula(paste0("survival::Surv(time, event) ~ ", z_terms)),
          data = cox_df
        )
        m2_cindex <- cox_fit$concordance[6]
        m2_se <- cox_fit$concordance[7]
      }
    }, error = function(e) NULL)

    per_ds_rows[[ds_name]] <- tibble::tibble(
      k = rep(k_val, 2), alpha = rep(alpha_val, 2),
      ntop = rep(ntop_val, 2), dataset = rep(ds_name, 2),
      method = c("transfer_beta", "refit_cox"),
      cindex = c(m1_cindex, m2_cindex),
      cindex_se = c(m1_se, m2_se),
      n_samples = rep(n_samp, 2),
      n_events = rep(n_evt, 2)
    )

    # Collect for pooling
    pooled_lp <- c(pooled_lp, lp_val)
    pooled_Z[[ds_name]] <- Z_val
    pooled_time <- c(pooled_time, time_val)
    pooled_event <- c(pooled_event, event_val)
    pooled_ds_label <- c(pooled_ds_label, rep(ds_name, n_samp))
  }

  # Pooled results
  pooled_row <- tibble::tibble(
    k = rep(k_val, 2), alpha = rep(alpha_val, 2),
    ntop = rep(ntop_val, 2), dataset = rep("pooled", 2),
    method = c("transfer_beta", "refit_cox"),
    cindex = rep(NA_real_, 2), cindex_se = rep(NA_real_, 2),
    n_samples = rep(length(pooled_time), 2),
    n_events = rep(sum(pooled_event), 2)
  )

  if (length(pooled_time) >= 2 && sum(pooled_event) > 0) {
    # Pooled Method 1: transfer_beta
    tryCatch({
      if (!anyNA(pooled_lp) && !any(is.infinite(pooled_lp))) {
        cc <- survival::concordance(
          survival::Surv(pooled_time, pooled_event) ~ pooled_lp,
          reverse = TRUE
        )
        pooled_row$cindex[1] <- cc$concordance
        pooled_row$cindex_se[1] <- sqrt(cc$var)
      }
    }, error = function(e) NULL)

    # Pooled Method 2: refit_cox with strata(dataset)
    tryCatch({
      Z_pooled <- do.call(rbind, pooled_Z)
      cox_df <- as.data.frame(Z_pooled)
      colnames(cox_df) <- paste0("Z", seq_len(ncol(Z_pooled)))
      cox_df$time <- pooled_time
      cox_df$event <- pooled_event
      cox_df$dataset <- factor(pooled_ds_label)
      z_terms <- paste(colnames(cox_df)[seq_len(ncol(Z_pooled))], collapse = " + ")
      cox_fit <- survival::coxph(
        as.formula(paste0("survival::Surv(time, event) ~ ", z_terms,
                          " + strata(dataset)")),
        data = cox_df
      )
      pooled_row$cindex[2] <- cox_fit$concordance[6]
      pooled_row$cindex_se[2] <- cox_fit$concordance[7]
    }, error = function(e) NULL)
  }

  dplyr::bind_rows(c(per_ds_rows, list(pooled_row)))
}

#' Aggregate validation results from all grid points
#'
#' @param result_list List of tibbles from validate_grid_point calls
#' @return Combined tibble
aggregate_cv_grid_val_results <- function(result_list) {
  dplyr::bind_rows(purrr::compact(result_list))
}

# Null-coalescing operator if not already defined
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
