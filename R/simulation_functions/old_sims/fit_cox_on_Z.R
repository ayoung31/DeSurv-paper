library(survival)

fit_cox_on_Z <- function(Z, surv_df) {
  # Z: N x K matrix
  # surv_df: data.frame with columns time, status, rownames matching rows of Z (or we'll assume same order)
  df <- data.frame(
    time = surv_df$time,
    status = surv_df$status,
    Z
  )
  # build formula: Surv(time, status) ~ Z1 + Z2 + ... + ZK
  k <- ncol(Z)
  covariates <- paste0("X", seq_len(k))
  colnames(df)[-(1:2)] <- covariates
  form <- as.formula(paste("Surv(time, status) ~", paste(covariates, collapse = " + ")))
  fit <- coxph(form, data = df, ties = "breslow")
  list(
    fit = fit,
    lp  = predict(fit, type = "lp") # linear predictor for those N samples
  )
}