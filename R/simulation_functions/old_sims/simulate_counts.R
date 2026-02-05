simulate_counts <- function(X_mean,
                            lib_mu = 1e5, lib_sdlog = 0.5,
                            dispersion = 0.2) {
  # X_mean: G x N
  G <- nrow(X_mean)
  N <- ncol(X_mean)
  
  # library sizes (scaling factors)
  lib_size <- rlnorm(N, meanlog = log(lib_mu), sdlog = lib_sdlog)
  
  counts <- matrix(NA, nrow = G, ncol = N)
  for (i in seq_len(N)) {
    mu_i <- X_mean[, i] * lib_size[i]
    
    # Negative binomial parameterization:
    # Var = mu + mu^2 * dispersion
    # size = 1/dispersion, prob = size / (size + mu)
    size  <- 1 / dispersion
    prob  <- size / (size + mu_i)
    
    counts[, i] <- rnbinom(G, size = size, prob = prob)
  }
  
  rownames(counts) <- rownames(X_mean)
  colnames(counts) <- colnames(X_mean)
  
  list(counts = counts, lib_size = lib_size)
}
