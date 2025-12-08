#' Simulate a nuisance-dominated scenario
#'
#' Extends `simulate_desurv_subtle()` by injecting strong batch and purity
#' effects that obscure the lethal program, mimicking real-world confounding.
#'
#' @inheritParams simulate_desurv_subtle
#' @param batch_frac Fraction of samples affected by the batch (0.2-0.8).
#' @param batch_amp Mean log-expression shift for the batch (1-5).
#' @param purity_amp Multiplicative scaling factor for purity shifts (0.5-3).
#' @param batch_shift_sd,purity_sdlog Spread of the nuisance effects.
#'
#' @return Same list structure with `scenario = "nuisance"`.
simulate_desurv_nuisance <- function(
    G = 5000,
    N = 200,
    K = 4,
    subtle_prog = 4,
    subtle_scale = 0.2,
    dominant_scale = 2.5,
    batch_frac = 0.5,
    batch_amp  = 3,
    purity_amp = 2,
    batch_shift_sd = 0.5,
    purity_sdlog = 0.5,
    lethal_effect = 2.0,
    background_effect = 0.0,
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
  # biological programs
  W_true <- simulate_W_subtle(
    G = G,
    K = K,
    subtle_prog = subtle_prog,
    subtle_scale = subtle_scale,
    other_scale = dominant_scale
  )
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
  
  X_mean_clean <- simulate_X_mean(W_true, H_true)
  X_mean_dirty <- add_nuisance(
    X_mean_clean,
    batch_frac = batch_frac,
    batch_amp  = batch_amp,
    purity_amp = purity_amp,
    batch_shift_sd = batch_shift_sd,
    purity_sdlog = purity_sdlog
  )
  
  sim_counts <- simulate_counts(
    X_mean_dirty,
    lib_mu = lib_mu,
    lib_sdlog = lib_sdlog,
    dispersion = dispersion
  )
  
  # survival driven by subtle program
  beta <- rep(background_effect, K)
  beta[subtle_prog] <- lethal_effect
  surv_df <- simulate_survival_linear(
    X = X_mean_clean,
    W = W_true,
    beta = beta,
    baseline_hazard = baseline_hazard,
    censor_rate = censor_rate
  )
  
  list(
    counts = sim_counts$counts,
    surv   = surv_df,
    W_true = W_true,
    H_true = H_true,
    beta_true = beta,
    marker_info = marker_info,
    scenario  = "nuisance"
  )
}
