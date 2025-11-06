simulate_W_easy <- function(G = 5000, K = 4, big_prog = 1) {
  W <- simulate_W(G, K)
  
  # Make prog1 SUPER strong everywhere (bigger amplitudes)
  W[, big_prog] <- W[, big_prog] * 3  # crank up variance contribution
  
  W
}