project_scores <- function(X, W) {
  # X: G x N
  # W: G x K
  # return: N x K
  t(X) %*% W
}