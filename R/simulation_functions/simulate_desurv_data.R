simulate_desurv_data <- function(
    G = 5000, N = 200, K = 4,
    beta = c(1.0, -0.5, 0.8, -1.2)
) {
  W_true <- simulate_W(G = G, K = K)
  H_true <- simulate_H(N = N, K = K, correlation = TRUE)
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(X_mean)
  
  surv_df <- simulate_survival(H_true, beta = beta)
  
  list(
    counts = sim_counts$counts,   # raw-ish RNA-seq counts
    lib_size = sim_counts$lib_size,
    surv = surv_df,               # time, status
    W_true = W_true,
    H_true = H_true,
    beta_true = beta
  )
}