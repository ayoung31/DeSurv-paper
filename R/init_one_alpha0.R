# One randomized init at alpha=0 (short run)
init_one_alpha0 <- function(X, y, delta, param_grid, path) {
  
  p=param_grid
  
  train_prefix = p$train_prefix
  method_trans_train = p$method_trans_train
  ngene = p$ngene
  maxit = p$maxit
  tol = p$tol
  imaxit = p$imaxit
  init = p$init
  k       = p$k
  lambda  = p$lambda
  eta     = p$eta
  lambdaW = p$lambdaW
  lambdaH = p$lambdaH
  
  fit = tryCatch({
    run_coxNMF(
      X = X, y = y, delta = delta, k = k,
      alpha = 0, lambda = lambda, eta = eta, lambdaW = lambdaW, lambdaH = lambdaH,
      seed = init, tol = tol, maxit = maxit, verbose = FALSE, ninit = 1, imaxit = imaxit
    )
  }, error = function(e) NULL
  )
    
  
  
  score = tibble(
    train_prefix = train_prefix,
    method_trans_train = method_trans_train,
    ngene     = ngene,
    maxit     = maxit,
    tol       = tol,
    imaxit    = imaxit,
    k         = k,
    lambda    = lambda,
    eta       = eta,
    lambdaW   = lambdaW,
    lambdaH   = lambdaH,
    seed      = init
  )
  
  if(!is.null(fit)){
    score$surv_loss = -1 * fit$loss$surv_loss
    score$nmf_loss  = fit$loss$nmf_loss
    score$flag_nan  = fit$`NaN flag`
  }else{
    score$surv_loss = NA
    score$nmf_loss = NA
    score$flag_nan = TRUE
  }
  
  tmp = paste0(path,".tmp")
  saveRDS(score,tmp,compress="xz")
  file.rename(tmp,path)
  path
}