# Compute metrics from a bundle file (no refitting).
compute_metrics <- function(bundle_path, X, y, delta) {
  b <- readRDS(bundle_path)
  stopifnot(is.list(b$fits), length(b$fits) > 0)
  pars <- b$meta
  out <- vector("list", length(b$fits))
  j <- 1
  for (nm in sort(names(b$fits))) {
    a <- as.numeric(nm)
    fit <- b$fits[[nm]]
    W <- fit$W; beta <- fit$beta
    lp <- t(X) %*% W %*% beta
    cix <- cvwrapr::getCindex(lp, survival::Surv(y, delta))
    
    loss <- fit$loss
    sl <- loss$surv_loss; nl <- loss$nmf_loss
    mrow <- data.frame(
      k = pars$k, alpha = a, lambda = pars$lambda, eta = pars$eta,
      lambdaW = pars$lambdaW, lambdaH = pars$lambdaH,
      c = cix,
      loss = loss$loss,
      sloss = sl, nloss = nl,
      sloss_std = sl / loss$std_surv,
      nloss_std = nl / loss$std_nmf,
      pen_beta = loss$penalty_beta, penW = loss$penalty_W, penH = loss$penalty_H,
      bic = -2*sl + pars$k * log(ncol(X)),
      converged = fit$iter < (pars$maxit %||% Inf),
      niter = fit$iter,
      flag_nan = isTRUE(fit$nan_flag),
      stringsAsFactors = FALSE
    )
    out[[j]] <- mrow; j <- j + 1
  }
  dplyr::bind_rows(out)
}
`%||%` <- function(a, b) if (!is.null(a)) a else b