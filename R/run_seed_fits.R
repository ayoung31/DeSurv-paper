run_seed_fits <- function(data_filtered, params_best, lambdaW, lambdaH) {
  seeds <- seq_len(NINIT_FULL)
  fits <- vector("list", length(seeds))
  cindex <- rep(NA_real_, length(seeds))
  for (i in seq_along(seeds)) {
    fit_i <- try(
      desurv_fit(
        X = data_filtered$ex,
        y = data_filtered$sampInfo$time,
        d = data_filtered$sampInfo$event,
        k = params_best$k,
        alpha = params_best$alpha,
        lambda = params_best$lambda,
        nu = params_best$nu,
        lambdaW = lambdaW,
        lambdaH = lambdaH,
        seed = seeds[i],
        tol = TOL / 100,
        tol_init = TOL,
        maxit = MAXIT,
        imaxit = MAXIT,
        ninit = 1,
        parallel_init = FALSE,
        verbose = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(fit_i, "try-error") && inherits(fit_i, "desurv_fit")) {
      fits[[i]] <- fit_i
      cindex[i] <- if (!is.null(fit_i$cindex)) fit_i$cindex else NA_real_
    }
  }
  keep <- which(!vapply(fits, is.null, logical(1)))
  if (!length(keep)) {
    stop("No successful seed fits were obtained.")
  }
  list(fits = fits[keep], seeds = seeds[keep], cindex = cindex[keep])
}
