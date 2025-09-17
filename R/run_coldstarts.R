run_coldstarts <- function(
    X, y, delta,
    params,
    alpha_vec,
    verbose = TRUE,
    path
){

  maxit              = MAXIT
  tol                = TOL
  imaxit             = IMAXIT
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  alpha              = params$alpha
  best_seed          = params$seed

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  

  if (verbose) message(sprintf("fit: k=%d λ=%.4g η=%.4g λW=%.4g λH=%.4g α=%.3f",
                               k, lambda, eta, lambdaW, lambdaH, alpha))
  

    if (!is.null(best_seed)) {
      fit <- run_coxNMF(
        X=X, y=y, delta=delta, k=k,
        alpha=alpha, lambda=lambda, eta=eta,
        lambdaW=lambdaW, lambdaH=lambdaH,
        seed=best_seed, tol=tol, maxit=maxit, verbose=verbose,
        ninit=1, imaxit=imaxit
      )
    } else {
      # fall back to your original multi-start behavior
      fit <- run_coxNMF(
        X=X, y=y, delta=delta, k=k,
        alpha=alpha, lambda=lambda, eta=eta,
        lambdaW=lambdaW, lambdaH=lambdaH,
        seed=123, tol=tol, maxit=maxit, verbose=verbose,
        ninit=100, imaxit=imaxit
      )
      warning("best seed not found, starting at seed 123")
    }


  
  # assemble / merge bundle (unchanged, plus best_seed in meta)
  meta <- list(
    k=k, lambda=lambda, eta=eta, lambdaW=lambdaW, lambdaH=lambdaH,
    alpha = alpha, imaxit=imaxit, maxit=maxit, tol=tol,
    ngene = nrow(X), nsamples = ncol(X),
    best_seed = best_seed
  )

  bundle <- list(meta = meta, fit = fit)

  
  tmp <- paste0(path, ".tmp")
  saveRDS(bundle, tmp, compress = "xz")
  file.rename(tmp, path)
  path
}