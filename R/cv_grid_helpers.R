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

# Null-coalescing operator if not already defined
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
