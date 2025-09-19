run_warmstarts <- function(
    X, y, delta,
    params,
    alpha_vec,
    verbose = TRUE,
    path
){
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  seeds              = 1:NINIT
  flag_exist         = FALSE

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if(file.exists(path)){
    bundle_old = readRDS(path)
    exists = unlist(lapply(bundle_old,function(x) as.numeric(names(x$fits))==alpha_vec))
    if(all(exists)){
      if(length(bundle_old)==NINIT){
        return(path)
      }else{
        seeds = setdiff(1:NINIT,1:length(bundle_old))
        flag_exist = TRUE
      }
    }
    
  }
  bundle=list()
  for(i in seeds){
    fit_prev = NULL
    fits = list()
    for (a in alpha_vec) {
  
      if (is.null(fit_prev)) {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          seed=best_seed, tol=tol, maxit=maxit, verbose=verbose,
          ninit=1, imaxit=imaxit
        )
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
      
      ## compute metrics train
      
      ## compute metrics test
      
      fits[[as.character(a)]] <- fit
      fit_prev <- fit
    }
  
    
    # assemble / merge bundle (unchanged, plus best_seed in meta)
    meta <- list(
      k=k, lambda=lambda, eta=eta, lambdaW=lambdaW, lambdaH=lambdaH,
      alpha_vec = alpha_vec, imaxit=imaxit, maxit=maxit, tol=tol,
      ngene = nrow(X), nsamples = ncol(X),
      seed = seed
    )
  
    bundle[[i]] <- list(meta = meta, fits = fits)
  }
  
  if(flag_exists){
    bundle = c(bundle_old,bundle)
  }
  
  tmp <- paste0(path, ".tmp")
  saveRDS(bundle, tmp, compress = "xz")
  file.rename(tmp, path)
  path
}