simulate_H <- function(
    N = 200,
    K = 4,
    correlated_pairs = list(c(1,3), c(2,4)),  # which programs tend to co-occur
    shape_bg = 2,
    scale_bg = 2,
    shape_noise = 1,
    scale_noise = 1
) {
  H <- matrix(NA, nrow = K, ncol = N)
  
  # For each correlated pair, generate a shared base signal
  for (pair in correlated_pairs) {
    base_signal <- rgamma(N, shape = shape_bg, scale = scale_bg)
    for (k in pair) {
      H[k, ] <- base_signal + rgamma(N, shape = shape_noise, scale = scale_noise)
    }
  }
  
  # If K has programs not covered in correlated_pairs, fill them independently
  uncovered <- setdiff(seq_len(K), unlist(correlated_pairs))
  for (k in uncovered) {
    if (all(is.na(H[k, ]))) {
      H[k, ] <- rgamma(N, shape = shape_bg, scale = scale_bg)
    }
  }
  
  rownames(H) <- paste0("prog", seq_len(K))
  colnames(H) <- paste0("patient", seq_len(N))
  H
}
