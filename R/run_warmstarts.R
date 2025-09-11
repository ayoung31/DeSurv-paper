run_warmstarts <- function(
    X, y, delta,
    params,
    alpha_vec,
    verbose = TRUE,
    path
){

  maxit              = params$maxit
  tol                = params$tol
  imaxit             = params$imaxit
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  best_seed          = params$seed

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  fit_prev = NULL
  fits = list()
  for (a in alpha_vec) {
    if (verbose) message(sprintf("fit: k=%d λ=%.4g η=%.4g λW=%.4g λH=%.4g α=%.3f",
                                 k, lambda, eta, lambdaW, lambdaH, a))
    
    if (is.null(fit_prev)) {
      # --- NEW: if no warm start available at the start of a block,
      #     re-run a single seeded init (the best seed you selected upstream).
      if (!is.null(best_seed)) {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          seed=best_seed, tol=tol, maxit=maxit, verbose=verbose,
          ninit=1, imaxit=imaxit
        )
      } else {
        # fall back to your original multi-start behavior
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          seed=123, tol=tol, maxit=maxit, verbose=verbose,
          ninit=100, imaxit=imaxit
        )
        warning("best seed not found, starting at seed 123")
      }
    } else {
      fit <- run_coxNMF(
        X=X, y=y, delta=delta, k=k,
        alpha=a, lambda=lambda, eta=eta,
        lambdaW=lambdaW, lambdaH=lambdaH,
        tol=tol, maxit=maxit, verbose=verbose,
        ninit=1, imaxit=imaxit,
        W0=fit_prev$W, H0=fit_prev$H, beta0=fit_prev$beta
      )
    }
    
    fits[[as.character(a)]] <- fit
    fit_prev <- fit
  }

  
  # assemble / merge bundle (unchanged, plus best_seed in meta)
  meta <- list(
    k=k, lambda=lambda, eta=eta, lambdaW=lambdaW, lambdaH=lambdaH,
    alpha_vec = alpha_vec, imaxit=imaxit, maxit=maxit, tol=tol,
    ngene = nrow(X), nsamples = ncol(X),
    best_seed = best_seed
  )

  bundle <- list(meta = meta, fits = fits)

  
  tmp <- paste0(path, ".tmp")
  saveRDS(bundle, tmp, compress = "xz")
  file.rename(tmp, path)
  path
}