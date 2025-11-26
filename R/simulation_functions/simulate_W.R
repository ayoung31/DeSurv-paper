simulate_W <- function(G = 5000, K = 4,
                       n_marker = 200,
                       high_mean = 8,  # log-scale
                       low_mean = 2,
                       noise_sd = 0.5) {
  
  W <- matrix(low_mean, nrow = G, ncol = K)
  
  for (k in seq_len(K)) {
    marker_genes <- sample(seq_len(G), n_marker)
    W[marker_genes, k] <- high_mean
  }
  
  # Ensure nonnegativity and maybe add small noise
  W <- pmax(W + matrix(rnorm(G * K, sd = noise_sd), G, K), 0)
  
  rownames(W) <- paste0("gene", seq_len(G))
  colnames(W) <- paste0("prog", seq_len(K))
  W
}
