simulate_W_subtle <- function(G = 5000, K = 4, subtle_prog = 4) {
  W <- simulate_W(G, K)
  
  # Make prog4 tiny (rare/low contribution to X)
  W[, subtle_prog] <- W[, subtle_prog] * 0.2
  
  # Make other programs big/noisy so they dominate variance
  W[, setdiff(seq_len(K), subtle_prog)] <-
    W[, setdiff(seq_len(K), subtle_prog)] * 2.5
  
  W
}