simulate_survival_interaction <- function(X,
                                          W,
                                          baseline_hazard = 0.01,
                                          censor_rate = 0.2,
                                          pair = c(1,3),
                                          gamma = 1.5) {
  # pair = which programs interact to kill you
  stopifnot(nrow(X) == nrow(W))
  stopifnot(all(pair %in% seq_len(ncol(W))))
  patient_scores <- crossprod(X, W)
  h1 <- patient_scores[, pair[1]]
  h2 <- patient_scores[, pair[2]]
  
  # nonlinear risk: high only if both high
  linpred <- gamma * (h1 * h2)
  
  n <- length(linpred)
  event_time <- rexp(n, rate = baseline_hazard * exp(linpred))
  time <- event_time
  status <- rep(1L, n)
  if (censor_rate > 0) {
    ncensor <- round(censor_rate * n)
    ncensor <- max(0L, min(n, ncensor))
    if (ncensor > 0) {
      censor_idx <- sample.int(n, ncensor)
      time[censor_idx] <- runif(
        n = ncensor,
        min = 0,
        max = event_time[censor_idx]
      )
      status[censor_idx] <- 0L
    }
  }
  
  data.frame(
    patient = colnames(X),
    time = time,
    status = status,
    linpred_true = linpred  # true risk driver is interaction
  )
}
