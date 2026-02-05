compute_cindex <- function(lp, surv_df) {
  # lp: numeric vector of risk scores (higher = riskier)
  # surv_df: same N rows, with time, status
  cc <- survConcordance(Surv(surv_df$time, surv_df$status) ~ lp)
  # C-index = concordant / comparable
  cc$concordance
}