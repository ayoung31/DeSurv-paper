add_nuisance <- function(X_mean,
                         batch_frac = 0.5,
                         batch_amp  = 3,
                         purity_amp = 2) {
  G <- nrow(X_mean)
  N <- ncol(X_mean)
  
  # Batch effect: pick ~50% of samples as "batch B"
  batch_group <- rbinom(N, size = 1, prob = batch_frac)
  # Gene-wise batch shift for that group
  batch_shift_gene <- rnorm(G, mean = batch_amp, sd = 0.5)
  batch_shift_mat <- batch_shift_gene %*% t(batch_group)
  
  # Purity-like scaling: some samples get globally higher tumor-like signal
  purity_factor <- 1 + rlnorm(N, meanlog = 0, sdlog = 0.5) * purity_amp
  
  X_dirty <- sweep(X_mean, 2, purity_factor, `*`) + batch_shift_mat
  pmax(X_dirty, 0)
}
