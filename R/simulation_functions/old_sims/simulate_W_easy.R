simulate_W_easy <- function(G = 5000, K = 4, big_prog = 1,
                            high_mean = 8, low_mean=2,
                            big_prog_multiplier = 3, noise_sd=1) {
  W <- simulate_W(
    G, K, 
    high_mean = high_mean,
    low_mean = low_mean,
    noise_sd = noise_sd)
  marker_info <- attr(W, "marker_info")
  
  # Make prog1 SUPER strong everywhere (bigger amplitudes)
  W[, big_prog] <- W[, big_prog] * big_prog_multiplier  # crank up variance contribution
  
  attr(W, "marker_info") <- marker_info
  
  W
}
