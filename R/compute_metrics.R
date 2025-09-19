# Compute metrics from a bundle file (no refitting).
compute_metrics <- function(fit, X, y, delta, alpha, fold) {
  W <- fit$W
  beta <- fit$beta
  lp <- t(X) %*% W %*% beta
  cix <- cvwrapr::getCindex(lp, survival::Surv(y, delta))
  
  loss <- fit$loss
  sl <- loss$surv_loss
  nl <- loss$nmf_loss
  met <- data.frame(
    alpha = alpha,
    fold= fold,
    c = cix,
    loss = loss$loss,
    sloss = sl, 
    nloss = nl,
    bic = -2*sl + ncol(W) * log(ncol(X)),
    converged = fit$convergence,
    niter = fit$iter,
    flag_nan = fit$nan_flag
  )
  return(met)
}