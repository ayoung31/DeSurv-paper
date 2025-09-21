# Compute metrics from a bundle file (no refitting).
compute_metrics <- function(fit, X, y, delta, alpha, fold, test=FALSE) {
  W <- fit$W
  beta <- fit$beta
  XtW = t(X) %*% W
  ns = which(fit$sdZ > 1e-12)
  XtW = XtW[,ns]
  meanZ = fit$meanZ[,ns]
  sdZ = fit$meanZ[,ns]
  beta = beta[ns]
  XtW = sweep(XtW,2,meanZ,FUN="-")
  XtW = sweep(XtW,2,sdZ,FUN="/")
  lp <- XtW %*% (beta*sdZ)
  cix <- cvwrapr::getCindex(lp, survival::Surv(y, delta))
  dat = data.frame(matrix(nrow=length(lp),ncol=0))
  dat$lp = lp
  dat$y = y
  dat$delta = delta
  
  fit_null = survival::coxph(survival::Surv(y,delta)~1,data=dat,ties = "breslow")
  fit_new = survival::coxph(survival::Surv(y,delta)~offset(lp),data=dat,
                            control = coxph.control(iter.max = 0),ties="breslow")
  if(test){
    sl = 2*logLik(fit_new)/sum(delta)
    dev=-2*(fit_new$loglik - fit_null$loglik)/sum(delta)
    ol = NA
    nl = NA
  }else{
    loss <- fit$loss
    ol = loss$loss
    sl <- loss$surv_loss
    dev=-2*(fit_new$loglik - fit_null$loglik)/sum(delta)
    nl <- loss$nmf_loss
  }

  met <- data.frame(
    alpha = alpha,
    fold= fold,
    c = cix,
    loss = ol,
    sloss = sl,
    dev = dev,
    nloss = nl,
    bic = -2*sl + ncol(W) * log(ncol(X)),
    converged = fit$convergence,
    niter = fit$iter,
    flag_nan = fit$nan_flag
  )
  return(met)
}
