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
  eta     = if (!is.null(p$eta)) p$eta else p$nu
  lambdaW = p$lambdaW
  lambdaH = p$lambdaH
  nu      = if (!is.null(p$nu)) p$nu else eta
  
  fit = tryCatch({
    DeSurv::desurv_fit(
      X = X, y = y, d = delta, k = k,
      alpha = 0, lambda = lambda, nu = nu, lambdaW = lambdaW, lambdaH = lambdaH,
      seed = init, tol = tol, maxit = maxit, verbose = FALSE, ninit = 1, imaxit = imaxit
    )
  }, error = function(e) NULL)
    
  
  
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
    surv_metric <- if (!is.null(fit$cindex)) fit$cindex else NA_real_
    score$surv_loss = surv_metric
    score$nmf_loss  = NA_real_
    score$sloss     = surv_metric
    score$nloss     = NA_real_
    score$loss      = -surv_metric
    score$flag_nan  = !isTRUE(fit$convergence)
  }else{
    score$surv_loss = NA_real_
    score$nmf_loss = NA_real_
    score$sloss = NA_real_
    score$nloss = NA_real_
    score$loss = NA_real_
    score$flag_nan = TRUE
  }
  
  tmp = paste0(path,".tmp")
  saveRDS(score,tmp,compress="xz")
  file.rename(tmp,path)
  path
}
