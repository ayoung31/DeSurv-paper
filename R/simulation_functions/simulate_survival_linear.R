simulate_survival_linear <- function(
    X,
    W,
    beta,                # length K, log-hazard weights per program
    baseline_hazard = 0.01,
    censor_rate = 0.2    # fraction of samples to censor
) {
  stopifnot(nrow(X) == nrow(W))
  stopifnot(length(beta) == ncol(W))
  patient_scores <- crossprod(X, W)        # N x K
  # linear predictor (risk score) built from X^T W
  linpred <- as.vector(patient_scores %*% beta)
  
  linpred = linpred / 1e6
  # event times ~ exponential with rate = baseline_hazard * exp(linpred)
  event_time <- rexp(
    n = length(linpred),
    rate = baseline_hazard * exp(linpred)
  )
  
  n <- length(linpred)
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
    linpred_true = linpred
  )
}
