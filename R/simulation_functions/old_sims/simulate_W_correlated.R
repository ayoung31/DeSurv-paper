simulate_W_correlated <- function(
    G = 5000, K = 4,
    overlap_pair = c(1,3),
    n_marker_shared = 150,
    n_marker_unique = 50,
    high_mean = 8,
    low_mean = 2,
    noise_sd = 0.5
) {
  W <- matrix(low_mean, nrow = G, ncol = K)
  
  # shared CAF-like block for prog1 and prog3
  shared_genes <- sample(seq_len(G), n_marker_shared)
  
  # unique subsets
  unique_1 <- sample(setdiff(seq_len(G), shared_genes), n_marker_unique)
  unique_3 <- sample(setdiff(seq_len(G), c(shared_genes, unique_1)), n_marker_unique)
  
  k1 <- overlap_pair[1]
  k2 <- overlap_pair[2]
  
  # fill in high expression for shared genes
  W[shared_genes, k1] <- high_mean
  W[shared_genes, k2] <- high_mean * 0.9  # almost same profile
  
  # add some private markers so they're not 100% identical
  W[unique_1, k1] <- high_mean
  W[unique_3, k2] <- high_mean
  
  # other programs get their own markers (low correlation)
  other_k <- setdiff(seq_len(K), c(k1, k2))
  for (k in other_k) {
    mk <- sample(setdiff(seq_len(G),
                         c(shared_genes, unique_1, unique_3)),
                 n_marker_shared)
    W[mk, k] <- high_mean
  }
  
  # add noise and floor at 0
  W <- pmax(W + matrix(rnorm(G*K, sd = noise_sd), G, K), 0)
  
  rownames(W) <- paste0("gene", seq_len(G))
  colnames(W) <- paste0("prog", seq_len(K))
  
  W
}