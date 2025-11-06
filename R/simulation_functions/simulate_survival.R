simulate_survival <- function(H, beta = c(1.0, -0.5, 0.8, -1.2),
                              baseline_hazard = 0.01,
                              censor_rate = 0.005) {
  
  K <- nrow(H)
  N <- ncol(H)
  stopifnot(length(beta) == K)
  
  linpred <- as.vector(t(beta) %*% H)  # length N
  
  # Exponential survival with rate = baseline_hazard * exp(linpred)
  # T ~ Exp(rate); so T = rexp(1, rate)
  event_time <- rexp(N, rate = baseline_hazard * exp(linpred))
  
  # Censoring times
  censor_time <- rexp(N, rate = censor_rate)
  
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)  # 1=event, 0=censored
  
  data.frame(
    patient = colnames(H),
    time = time,
    status = status,
    linpred_true = linpred
  )
}
