init_alpha0 <- function(X, y, delta, p, path, ninit, imaxit, maxit, tol) {
  # read existing results if present
  existing <- if (file.exists(path)) readRDS(path) else NULL
  existing_fits <- if (!is.null(existing) && !is.null(existing$fits)) existing$fits else list()
  
  # which seeds are already done? (only those with non-NULL fits)
  done_seeds <- {
    nm <- names(existing_fits)
    if (is.null(nm)) integer(0) else as.integer(nm[vapply(existing_fits, Negate(is.null), logical(1))])
  }
  
  # determine seeds still to run
  inits_new <- setdiff(seq_len(ninit), done_seeds)
  if (length(inits_new) == 0L) {
    return(existing %||% list(meta = list(), fits = existing_fits))
  }
  
  # unpack params (single combo)
  k       <- p$k
  lambda  <- p$lambda
  eta     <- p$eta
  lambdaW <- p$lambdaW
  lambdaH <- p$lambdaH
  
  # run new seeds
  fits_new <- vector("list", length(inits_new))
  names(fits_new) <- as.character(inits_new)
  
  for (j in seq_along(inits_new)) {
    seed_j <- inits_new[j]
    fits_new[[j]] <- tryCatch(
      run_coxNMF(
        X = X, y = y, delta = delta, k = k,
        alpha = 0, lambda = lambda, eta = eta,
        lambdaW = lambdaW, lambdaH = lambdaH,
        seed = seed_j, tol = tol, maxit = maxit,
        verbose = FALSE, ninit = 1, imaxit = imaxit
      ),
      error = function(e) NULL
    )
  }
  
  # drop NULLs (failed runs)
  keep <- !vapply(fits_new, is.null, logical(1))
  fits_new <- fits_new[keep]
  
  # merge with existing (overwrite by seed if re-run)
  fits_merged <- existing_fits
  if (length(fits_new)) {
    fits_merged[names(fits_new)] <- fits_new
  }
  
  out <- list(
    meta = list(
      k = k, lambda = lambda, eta = eta,
      lambdaW = lambdaW, lambdaH = lambdaH,
      imaxit = imaxit, maxit = maxit, tol = tol, ninit = ninit
    ),
    fits = fits_merged
  )
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_atomic_rds(out, path)
  
  .extract_losses <- function(fit) {
    if (is.null(fit) || is.null(fit$loss)) {
      return(c(nmf_loss = NA_real_, surv_loss = NA_real_, loss = NA_real_))
    }
    L <- fit$loss
    nmf  <- if (!is.null(L$nmf_loss)) L$nmf_loss else NA_real_
    surv <- if (!is.null(L$surv_loss)) L$surv_loss else NA_real_
    tot  <- if (!is.null(L$loss)) L$loss else NA_real_
    c(nmf_loss = nmf, surv_loss = surv, loss = tot)
  }
  
  # main: build a data.frame of losses + seeds
  if (!length(fits_merged)) {
    return(data.frame(seed=character(), nmf_loss=double(), surv_loss=double(), loss=double()))
  }
  
  seeds <- names(fits_merged)
  if (is.null(seeds)) seeds <- as.character(seq_along(fits_merged))
  
  mat <- t(sapply(fits_merged, .extract_losses))
  loss_df <- data.frame(seed = seeds, as.data.frame(mat, optional = TRUE), row.names = NULL)
  loss_df
  
}

`%||%` <- function(a, b) if (is.null(a)) b else a