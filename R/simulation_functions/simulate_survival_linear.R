simulate_survival_linear <- function(
    X,
    W,
    beta,                # length K, log-hazard weights per program
    baseline_hazard = 0.01,
    censor_rate = 0.2    # fraction of samples to censor
) {
  stopifnot(nrow(X) == nrow(W))
  stopifnot(length(beta) == ncol(W))

  # X: G x N, W: G x K
  patient_scores <- crossprod(X, W)   # N x K (S)
  
  # K <- ncol(patient_scores)
  # if (K < 2L) stop("Need K >= 2 for mean-of-others contrast.")
  # 
  # # Per-sample mean of *other* programs: N x K
  # other_means <- (rowSums(patient_scores) - patient_scores) / (K - 1)
  
  # Contrast: each program minus mean of the others in that sample
  patient_scores_contrast <- patient_scores #- other_means  # N x K
  
  # Optional: standardize each contrast column across patients (for stability)
  patient_scores_contrast <- scale(patient_scores_contrast)
  
  # Linear predictor (risk score) built from contrast(X^T W)
  linpred <- as.vector(patient_scores_contrast %*% beta)
  
  # Event times ~ exponential with rate = baseline_hazard * exp(linpred)
  event_time <- rexp(
    n    = length(linpred),
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
