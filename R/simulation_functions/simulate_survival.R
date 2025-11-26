simulate_survival <- function(H, beta = c(1.0, -0.5, 0.8, -1.2),
                              baseline_hazard = 0.01,
                              censor_rate = 0.2) {
  
  K <- nrow(H)
  N <- ncol(H)
  stopifnot(length(beta) == K)
  
  linpred <- as.vector(t(beta) %*% H)  # length N
  
  # Exponential survival with rate = baseline_hazard * exp(linpred)
  # T ~ Exp(rate); so T = rexp(1, rate)
  event_time <- rexp(N, rate = baseline_hazard * exp(linpred))
  
  time <- event_time
  status <- rep(1L, N)
  if (censor_rate > 0) {
    ncensor <- round(censor_rate * N)
    ncensor <- max(0L, min(N, ncensor))
    if (ncensor > 0) {
      censor_idx <- sample.int(N, ncensor)
      time[censor_idx] <- runif(
        n = ncensor,
        min = 0,
        max = event_time[censor_idx]
      )
      status[censor_idx] <- 0L
    }
  }
  
  data.frame(
    patient = colnames(H),
    time = time,
    status = status,
    linpred_true = linpred
  )
}
