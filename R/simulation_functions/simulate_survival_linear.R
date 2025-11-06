simulate_survival_linear <- function(
    H,
    beta,                # length K, log-hazard weights per program
    baseline_hazard = 0.01,
    censor_rate = 0.005
) {
  # linear predictor (risk score)
  linpred <- as.vector(t(beta) %*% H)  # length N
  
  # event times ~ exponential with rate = baseline_hazard * exp(linpred)
  event_time <- rexp(
    n = length(linpred),
    rate = baseline_hazard * exp(linpred)
  )
  
  # random censoring
  censor_time <- rexp(
    n = length(linpred),
    rate = censor_rate
  )
  
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)  # 1 = event, 0 = censored
  
  data.frame(
    patient = colnames(H),
    time = time,
    status = status,
    linpred_true = linpred
  )
}
