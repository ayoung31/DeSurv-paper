simulate_desurv_easy <- function(G=5000, N=200, K=4) {
  W_true <- simulate_W_easy(G, K, big_prog = 1)
  H_true <- simulate_H(N, K, correlation = TRUE)
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(X_mean)
  
  beta <- c(2.0, 0.0, 0.0, 0.0)
  surv_df <- simulate_survival_linear(H_true, beta = beta)
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = beta,
    scenario  = "easy"
  )
}
