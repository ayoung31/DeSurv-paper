# Compute metrics from a bundle file (no refitting).
compute_metrics <- function(path, X, y, delta) {
  b <- readRDS(path)
  pars <- b$meta
  
  fit <- b$fit
  W <- fit$W; beta <- fit$beta
  lp <- t(X) %*% W %*% beta
  cix <- cvwrapr::getCindex(lp, survival::Surv(y, delta))
  
  loss <- fit$loss
  sl <- loss$surv_loss
  nl <- loss$nmf_loss
  met <- data.frame(
    k = pars$k, alpha = pars$alpha, lambda = pars$lambda, eta = pars$eta,
    lambdaW = pars$lambdaW, lambdaH = pars$lambdaH,
    c = cix,
    loss = loss$loss,
    sloss = sl, nloss = nl,
    pen_beta = loss$penalty_beta, penW = loss$penalty_W, penH = loss$penalty_H,
    bic = -2*sl + pars$k * log(ncol(X)),
    converged = fit$iter < (pars$maxit %||% Inf),
    niter = fit$iter,
    flag_nan = isTRUE(fit$nan_flag),
    stringsAsFactors = FALSE
  )
  return(met)
}
`%||%` <- function(a, b) if (!is.null(a)) a else b