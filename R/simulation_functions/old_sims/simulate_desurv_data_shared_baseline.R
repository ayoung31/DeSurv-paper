#' Simulate data with a shared housekeeping baseline
#'
#' Models tumors where most genes follow a ubiquitous housekeeping program and
#' each biological program adds small deviations (e.g. lineage-specific bumps).
#' Survival is driven by a single lethal program layered on top of that baseline.
#'
#' @inheritParams simulate_desurv_data
#' @param lethal_prog Index of the lethal program (1-based).
#' @param lethal_effect Log-hazard weight for the lethal program (1-3 typical).
#' @param nonlethal_effect Log-hazard weight for other programs (-0.5-0.5).
#' @param baseline_mean,baseline_sd Mean/sd of the baseline expression backdrop.
#' @param bump_mean,bump_sd Magnitude of program-specific bumps (1-4 keeps them
#'   subtle).
#' @param bump_genes_per_prog Number of genes that receive a bump per program
#'   (100-400 captures sparse pathway deviations).
#' @param noise_sd_baseline Column-wise noise injected into W (0.1-0.5).
simulate_desurv_data_shared_baseline <- function(
    G = 5000,
    N = 200,
    K = 4,
    lethal_prog = 3,          # which program actually drives survival
    lethal_effect = 1.5,      # log-hazard weight for that program
    nonlethal_effect = 0.0,   # weight for the others
    baseline_mean = 5,
    baseline_sd   = 1,
    bump_mean     = 3,
    bump_sd       = 0.5,
    bump_genes_per_prog = 150,
    noise_sd_baseline = 0.2,
    lib_mu = 1e5,
    lib_sdlog = 0.5,
    dispersion = 0.2,
    baseline_hazard = 0.01,
    censor_rate = 0.2,
    correlated_pairs = list(c(1, 3), c(2, 4)),
    h_shape_bg = 2,
    h_scale_bg = 2,
    h_shape_noise = 1,
    h_scale_noise = 1
) {
  # 1. Program matrix W with shared baseline + program-specific bumps
  W_true <- simulate_W_shared_baseline(
    G = G,
    K = K,
    baseline_mean = baseline_mean,
    baseline_sd   = baseline_sd,
    bump_mean     = bump_mean,
    bump_sd       = bump_sd,
    bump_genes_per_prog = bump_genes_per_prog,
    noise_sd_baseline = noise_sd_baseline
  )
  marker_info <- attr(W_true, "marker_info")
  
  # 2. Patient loadings H (program usage per patient)
  H_true <- simulate_H(
    N = N,
    K = K,
    correlated_pairs = correlated_pairs,
    shape_bg = h_shape_bg,
    scale_bg = h_scale_bg,
    shape_noise = h_shape_noise,
    scale_noise = h_scale_noise
  )
  
  # 3. Expected expression X_mean = W %*% H
  X_mean <- simulate_X_mean(W_true, H_true)
  
  # 4. Noisy RNA-seq counts from X_mean
  sim_counts <- simulate_counts(
    X_mean,
    lib_mu = lib_mu,
    lib_sdlog = lib_sdlog,
    dispersion = dispersion
  )
  
  X_counts <- sim_counts$counts  # G x N
  
  # 5. Survival driven by ONE lethal program
  beta_vec <- rep(nonlethal_effect, K)
  beta_vec[lethal_prog] <- lethal_effect
  
  surv_df <- simulate_survival_linear(
    X = X_mean,
    W = W_true,
    beta = beta_vec,
    baseline_hazard = baseline_hazard,
    censor_rate = censor_rate
  )
  
  # Make sure survival rows align with columns of X
  # (they should, but let's be explicit)
  stopifnot(identical(colnames(X_counts), surv_df$patient))
  
  list(
    counts      = X_counts,     # G x N raw-ish counts
    surv        = surv_df,      # data.frame(patient,time,status,linpred_true)
    W_true      = W_true,       # G x K (shared baseline + bumps)
    H_true      = H_true,       # K x N
    beta_true   = beta_vec,     # length K, which program is lethal
    marker_info = marker_info,
    lethal_prog = lethal_prog,
    params      = list(
      G = G,
      N = N,
      K = K,
      lib_mu = lib_mu,
      dispersion = dispersion,
      lethal_prog = lethal_prog,
      lethal_effect = lethal_effect
    )
  )
}
