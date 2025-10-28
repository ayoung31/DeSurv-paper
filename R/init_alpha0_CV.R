init_alpha0_CV <- function(X, y, delta, param_grid, path, NINIT) {
  # ensure output directory exists
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  # read existing results if present
  existing_inits <- if (file.exists(path)) {
    readRDS(path)
  } else {
    NULL
  }
  
  # which seeds are already done?
  done_seeds <- if (!is.null(existing_inits) && "seed" %in% names(existing_inits)) {
    unique(existing_inits$seed)
  } else {
    integer(0)
  }
  
  # determine seeds still to run
  inits_new <- setdiff(seq_len(NINIT), done_seeds)
  if (length(inits_new) == 0L) {
    # nothing to do; keep file as-is
    return(invisible(path))
  }
  
  # unpack params (single combo)
  p <- param_grid
  k       <- p$k
  lambda  <- p$lambda
  eta     <- p$eta
  lambdaW <- p$lambdaW
  lambdaH <- p$lambdaH
  tol     <- p$tol
  maxit   <- p$maxit
  imaxit  <- p$imaxit
  
  # optional: stable combo id (useful for joins/debug)
  combo_id <- digest::digest(
    c(p$train_prefix, p$method_trans_train, p$ngene, k, lambda, eta, lambdaW, lambdaH, tol, maxit, imaxit),
    algo = "xxhash64"
  )
  
  res_list <- vector("list", length(inits_new))
  for (j in seq_along(inits_new)) {
    seed_j <- inits_new[j]
    
    fit <- tryCatch(
      run_coxNMF(
        X = X, y = y, delta = delta, k = k,
        alpha = 0, lambda = lambda, eta = eta, lambdaW = lambdaW, lambdaH = lambdaH,
        seed = seed_j, tol = tol, maxit = maxit, verbose = FALSE, ninit = 1, imaxit = imaxit
      ),
      error = function(e) NULL
    )
    
    # defensive extraction
    has_loss <- !is.null(fit) && !is.null(fit$loss)
    surv_loss <- if (has_loss) -1 * fit$loss$surv_loss else NA_real_
    nmf_loss  <- if (has_loss)        fit$loss$nmf_loss else NA_real_
    flag_nan  <- if (!is.null(fit) && !is.null(fit$`NaN flag`)) isTRUE(fit$`NaN flag`) else is.null(fit)
    
    res = tibble::tibble(
      combo_id = combo_id,
      seed = seed_j,
      surv_loss = surv_loss,
      nmf_loss  = nmf_loss,
      flag_nan  = flag_nan
    )
    
    res_list[[j]] <- dplyr::bind_cols(param_grid,res)
  }
  
  res_new <- dplyr::bind_rows(res_list)
  
  # append to existing, preserving column types
  out_tbl <- if (is.null(existing_inits)) res_new else dplyr::bind_rows(existing_inits, res_new)
  
  # atomic write
  tmp <- paste0(path, ".tmp")
  saveRDS(out_tbl, tmp, compress = "gzip")   # faster than "xz" for frequent writes
  ok <- file.rename(tmp, path)
  if (!ok) warning("file.rename() failed for ", path)
  
  invisible(path)
}