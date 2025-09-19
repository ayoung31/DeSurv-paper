run_warmstarts <- function(
    X, y, delta,
    params,
    verbose = TRUE,
    path_fits,
    path_mets
){
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  seeds              = 1:NINIT
  flag_exist         = FALSE

  dir.create(dirname(path_fits), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(path_mets), recursive = TRUE, showWarnings = FALSE)
  if(file.exists(path_fits) & file.exists(path_mets)){
    bundle_old = readRDS(path_fits)
    load(path_mets)
    exists = unlist(lapply(bundle_old,function(x) names(x$fits)==as.character(ALPHA)))
    if(all(exists)){
      if(length(bundle_old)==NINIT){
        return(path_fits)
      }else{
        seeds = setdiff(1:NINIT,1:length(bundle_old))
        mets_train_old=mets_train
        rm(mets_train)
        flag_exist = TRUE
      }
    }
    
  }
  bundle=list()
  mets_train=list()
  j=1
  for(i in seeds){
    fit_prev = NULL
    fits = list()
    for (a in ALPHA) {
      print(sprintf("seed %d, alpha %.2f",i,a))
      if (is.null(fit_prev)) {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          seed=i, tol=TOL, maxit=MAXIT, verbose=verbose,
          ninit=1, imaxit=MAXIT
        )
      } else {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          tol=TOL, maxit=MAXIT, verbose=verbose,
          ninit=1, imaxit=MAXIT,
          W0=fit_prev$W, H0=fit_prev$H, beta0=fit_prev$beta
        )
      }
      
      ## compute metrics train
      mets_train[[j]] = compute_metrics(fit,X,y,delta,a)


      
      fits[[as.character(a)]] <- fit
      fit_prev <- fit
      j=j+1
    }
  
    
    # assemble / merge bundle (unchanged, plus best_seed in meta)
    meta <- list(
      k=k, lambda=lambda, eta=eta, lambdaW=lambdaW, lambdaH=lambdaH,
      alpha_vec = ALPHA, maxit=MAXIT, tol=TOL,
      ngene = nrow(X), nsamples = ncol(X),
      seed = i
    )
  
    bundle[[i]] <- list(meta = meta, fits = fits)
  }
  
  mets_train = dplyr::bind_rows(mets_train)
  
  if(flag_exist){
    bundle = c(bundle_old,bundle)
    mets_train = dplyr::bind_rows(mets_train_old,mets_train)
  }
  

  save(mets_train,file=path_mets)
  
  tmp <- paste0(path_fits, ".tmp")
  saveRDS(bundle, tmp, compress = "xz")
  file.rename(tmp, path_fits)
  
  path_fits
}
