simulate_desurv_interaction <- function(G=5000, N=200, K=4) {
  W_true <- simulate_W(G, K)                 # normal
  H_true <- simulate_H(N, K, correlation=TRUE) # prog1+prog3 co-occur in some pts
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(X_mean)
  
  surv_df <- simulate_survival_interaction(
    H_true,
    pair = c(1,3),
    gamma = 1.5
  )
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = NA,          # not linear beta anymore
    scenario  = "interaction"
  )
}
