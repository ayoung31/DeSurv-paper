#' Simulate an "easy" single-program hazard scenario
#'
#' This simulator reflects a cohort where one canonical malignant program
#' (think high proliferation) drives both most of the expression variation and
#' essentially all of the mortality, while other programs behave as background
#' heterogeneity. Methods that cannot recover the dominant program should fail.
#'
#' @param G Number of genes (500-10000 recommended).
#' @param N Number of patients (100-500 keeps the task light-weight).
#' @param K Number of latent programs.
#' @param big_prog Index of the lethal/dominant program (1-based).
#' @param big_prog_multiplier Factor used inside `simulate_W_easy` to inflate
#'   the dominant program's amplitude (2-5 keeps it obvious but not absurd).
#' @param lethal_effect Log-hazard coefficient for the dominant program
#'   (1-3 gives meaningful separation).
#' @param background_effect Log-hazard coefficient for the non-dominant programs
#'   (keep at -0.5 to 0.5 to mimic neutral biology).
#' @param correlated_pairs,h_shape_bg,h_scale_bg,h_shape_noise,h_scale_noise
#'   Parameters passed to `simulate_H`; keep gamma shapes/scales between 1 and 5.
#' @param lib_mu,lib_sdlog,dispersion RNA-seq count parameters (same ranges as
#'   `simulate_desurv_data`).
#' @param baseline_hazard Base event rate.
#' @param censor_rate Fraction of subjects censored (0-1).
#'
#' @return A list with counts, survival data, true W/H, beta vector, and the
#'   scenario label `"easy"`.
simulate_desurv_controlled_cor <- function(
    G = 5000,
    N = 200,
    K = 4,
    big_prog = 1,
    lethal_effect = 2.0,
    background_effect = 0.0,
    high_mean = 8,
    low_mean = 2,
    noise_sd = 1,
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

  W_true = simulate_W_controlled_cor(
    G = G, K = K,
    n_marker = 100,
    base_sd = 1,
    mix_strength = 100,    # 0 = more orthogonal, >0 = more correlated
    shift = 0,           # baseline shift to make entries positive
    marker_boost = 2.0,    # multiply marker genes in their main program
    other_shrink = 0.5,    # multiply those marker genes in other programs
    noise_sd = 0.1
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
  beta[big_prog] <- lethal_effect
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
    X_mean = X_mean,
    marker_info = marker_info,
    scenario  = "easy"
  )
}
