# Compute metrics from a bundle file (no refitting).
compute_metrics <- function(fit, X, y, delta, test=FALSE) {
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
  dat = as.data.frame(XtW)
  colnames(dat) = paste0("xtw",1:ncol(XtW))
  dat$lp = t(X)%*%W%*%beta
  dat$y = y
  dat$delta = delta
  dat$lp_bin = dat$lp > median(dat$lp)
  
  sd_eta = sd(lp)

  fit_null = survival::coxph(survival::Surv(y,delta)~1,data=dat,ties = "breslow")
  fit_new = survival::coxph(survival::Surv(y,delta)~offset(lp),data=dat,
                            control = coxph.control(iter.max = 0),ties="breslow")
  fit_bin = survival::coxph(survival::Surv(y,delta)~lp_bin,data=dat)
  refit = survival::coxph(as.formula(paste0("Surv(y,delta) ~ ",
                                            paste(paste0("xtw",1:ncol(XtW)),collapse = " + "))),data=dat)
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
    sd_eta = sd_eta,
    c = cix,
    c_refit = refit$concordance[6],
    pl_refit = 2*logLik(refit)/sum(delta),
    hr = summary(fit_bin)[[8]][,1],
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
