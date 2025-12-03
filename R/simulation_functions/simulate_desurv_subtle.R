#' Simulate a "subtle" lethal-program scenario
#'
#' Here the biologically lethal program is low amplitude and often drowned out
#' by larger nuisance programs, mimicking an immune-suppressed or stromal
#' compartment whose signal is faint in bulk RNA-seq. Discovering the signal
#' requires methods that look for survival associations rather than pure
#' variance.
#'
#' @inheritParams simulate_desurv_easy
#' @param subtle_prog Index of the lethal but low-amplitude program.
#' @param subtle_scale Scale factor applied to the subtle program within
#'   `simulate_W_subtle` (0.05-0.5 keeps it tiny).
#' @param dominant_scale Scale factor for the other programs (1.5-4.0 keeps
#'   them dominant).
#'
#' @return Same structure as `simulate_desurv_easy` with `scenario = "subtle"`.
simulate_desurv_subtle <- function(
    G = 5000,
    N = 200,
    K = 4,
    subtle_prog = 4,
    subtle_scale = 0.2,
    dominant_scale = 2.5,
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
  X_mean <- simulate_X_mean(W_true, H_true)
  sim_counts <- simulate_counts(
    X_mean,
    lib_mu = lib_mu,
    lib_sdlog = lib_sdlog,
    dispersion = dispersion
  )
  
  beta <- rep(background_effect, K)
  beta[subtle_prog] <- lethal_effect
  surv_df <- simulate_survival_linear(
    X = X_mean,
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
    scenario  = "subtle"
  )
}
