simulate_X_mean <- function(W, H) {
  X_mean <- W %*% H  # G x N
  rownames(X_mean) <- rownames(W)
  colnames(X_mean) <- colnames(H)
  X_mean
}
