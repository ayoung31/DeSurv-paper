#' Simulate a synergistic hazard scenario
#'
#' Two programs (e.g. proliferation plus hypoxia) only become lethal when they
#' co-occur. The marginal effect of each program is weak; the interaction term
#' encoded by `gamma` dictates risk.
#'
#' @inheritParams simulate_desurv_data
#' @param interaction_pair Indices of the two programs that synergize.
#' @param interaction_gamma Strength of the interaction (1-3 works well).
#'
#' @return List with counts, survival data, true W/H, `beta_true = NA`, and the
#'   `"interaction"` scenario label.
simulate_desurv_interaction <- function(
    G = 5000,
    N = 200,
    K = 4,
    interaction_pair = c(1, 3),
    interaction_gamma = 1.5,
    n_marker = 200,
    high_mean = 8,
    low_mean = 2,
    w_noise_sd = 0.5,
    correlated_pairs = list(c(1, 3), c(2, 4)),
    h_shape_bg = 2,
    h_scale_bg = 2,
    h_shape_noise = 1,
    h_scale_noise = 1,
    lib_mu = 1e5,
    lib_sdlog = 0.5,
    dispersion = 0.2,
    baseline_hazard = 0.01,
    censor_rate = 0.2
) {
  W_true <- simulate_W(
    G = G,
    K = K,
    n_marker = n_marker,
    high_mean = high_mean,
    low_mean = low_mean,
    noise_sd = w_noise_sd
  )                 # normal
  marker_info <- attr(W_true, "marker_info")
  H_true <- simulate_H(
    N = N,
    K = K,
    correlated_pairs = correlated_pairs,
    shape_bg = h_shape_bg,
    scale_bg = h_scale_bg,
    shape_noise = h_shape_noise,
    scale_noise = h_scale_noise
  )
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(
    X_mean,
    lib_mu = lib_mu,
    lib_sdlog = lib_sdlog,
    dispersion = dispersion
  )
  
  surv_df <- simulate_survival_interaction(
    X = X_mean,
    W = W_true,
    pair = interaction_pair,
    gamma = interaction_gamma,
    baseline_hazard = baseline_hazard,
    censor_rate = censor_rate
  )
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = NA,          # not linear beta anymore
    marker_info = marker_info,
    scenario  = "interaction"
  )
}
