simulate_desurv_subtle <- function(G=5000, N=200, K=4) {
  W_true <- simulate_W_subtle(G, K, subtle_prog = 4)
  H_true <- simulate_H(N, K, correlation = TRUE)
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(X_mean)
  
  beta <- c(0.0, 0.0, 0.0, 2.0)
  surv_df <- simulate_survival_linear(H_true, beta = beta)
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = beta,
    scenario  = "subtle"
  )
}
