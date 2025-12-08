## -----------------------------
## 4. Simulate survival from XtW with marker masking
## -----------------------------
simulate_survival_from_XtW <- function(
    X,            # G x N (genes x samples)
    W,            # G x K
    marker_sets,  # list of length K, each element = gene names or indices
    beta,         # numeric vector of length K
    baseline_hazard = 0.05,
    censor_rate = 0.02,
    scale_scores = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  G <- nrow(X)
  N <- ncol(X)
  K <- ncol(W)
  
  stopifnot(length(marker_sets) == K)
  stopifnot(length(beta) == K)
  
  # Initialize Wtilde with zeros, same dimnames as W
  Wtilde <- matrix(0, nrow = G, ncol = K,
                   dimnames = dimnames(W))
  
  # Factor-specific masking: only keep own markers in each column
  for (k in seq_len(K)) {
    mk <- marker_sets[[k]]
    
    # If marker sets are gene names, match to rownames(W)
    if (is.character(mk)) {
      idx <- match(mk, rownames(W))
    } else {
      idx <- mk
    }
    
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= G]
    if (length(idx) > 0) {
      Wtilde[idx, k] <- W[idx, k]
    }
  }
  
  # Scores: N x K (samples x factors)
  scores <- crossprod(X, Wtilde)  # N x K
  
  if (scale_scores) {
    m <- colMeans(scores)
    s <- apply(scores, 2, sd)
    s[s == 0] <- 1
    scores <- sweep(scores, 2, m, "-")
    scores <- sweep(scores, 2, s, "/")
  }
  
  # linear predictor from supplied beta vector
  linpred <- as.vector(scores %*% beta)
  
  # Exponential survival, independent exponential censoring
  event_time  <- rexp(N, rate = baseline_hazard * exp(linpred))
  censor_time <- rexp(N, rate = censor_rate)
  
  time   <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  
  list(
    time   = time,
    status = status,
    linpred = linpred,
    scores  = scores,
    Wtilde  = Wtilde
  )
}

