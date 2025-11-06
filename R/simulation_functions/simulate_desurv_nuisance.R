simulate_desurv_nuisance <- function(G=5000, N=200, K=4) {
  # biological programs
  W_true <- simulate_W_subtle(G, K, subtle_prog = 4) 
  # subtle_prog=4 still low-ish magnitude biologically relevant prog
  H_true <- simulate_H(N, K, correlation=TRUE)
  
  X_mean_clean <- simulate_X_mean(W_true, H_true)
  X_mean_dirty <- add_nuisance(X_mean_clean,
                               batch_frac = 0.5,
                               batch_amp  = 3,
                               purity_amp = 2)
  
  sim_counts <- simulate_counts(X_mean_dirty)
  
  # survival driven by subtle program (prog4)
  beta <- c(0, 0, 0, 2.0)
  surv_df <- simulate_survival_linear(H_true, beta = beta)
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = beta,
    scenario  = "nuisance"
  )
}