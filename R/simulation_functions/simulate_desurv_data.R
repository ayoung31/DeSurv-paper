#' Simulate a baseline DeSurv cohort
#'
#' This simulator generates a "vanilla" bulk RNA-seq cohort where several
#' transcriptional programs contribute to expression and a user-specified
#' linear combination of those programs modulates hazard. The default resembles
#' a typical multi-lineage tumor where one proliferative program is lethal while
#' others are partially protective.
#'
#' @param G Number of genes (recommended 1000-20000 to mimic standard panels).
#' @param N Number of patients/samples (recommended 100-1000 for survival power).
#' @param K Number of latent programs (4-8 usually captures broad biology).
#' @param beta Length-K log-hazard weights; positive values increase risk.
#' @param n_marker Number of marker genes per program (use 50-300 for sparse
#'   programs, larger to make them diffuse).
#' @param high_mean Mean log-expression for marker genes (6-10 keeps clear
#'   signal while staying in bulk RNA-seq range).
#' @param low_mean Mean log-expression for non-marker genes (1-4 keeps them low).
#' @param w_noise_sd Gene-level noise injected into W to avoid identical columns
#'   (0.2-1 widens heterogeneity).
#' @param correlated_pairs List of program indices that co-occur in patients;
#'   keep length small (1-3 pairs) to mimic biologically coupled pathways.
#' @param h_shape_bg,h_scale_bg Shape/scale of the shared gamma background
#'   program intensity (values between 1 and 5 work well).
#' @param h_shape_noise,h_scale_noise Gamma parameters for per-program noise
#'   (0.5-3 gives realistic scatter).
#' @param lib_mu,lib_sdlog Mean and log-sd of library sizes (1e4-1e6 and
#'   0.2-0.8 cover most bulk assays).
#' @param dispersion Negative-binomial dispersion of the counts (0.05-0.5 for
#'   bulk RNA-seq).
#' @param baseline_hazard Base event rate for the exponential survival model
#'   (0.005-0.05 for multi-year horizons).
#' @param censor_rate Fraction of samples that should be censored (0-1).
#'
#' @return A list with raw counts, library sizes, survival data frame, true W/H,
#'   and `beta_true`.
simulate_desurv_data <- function(
    G = 5000, N = 200, K = 4,
    beta = c(1.0, -0.5, 0.8, -1.2),
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
  )
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
  
  surv_df <- simulate_survival(
    H_true,
    beta = beta,
    baseline_hazard = baseline_hazard,
    censor_rate = censor_rate
  )
  
  list(
    counts = sim_counts$counts,   # raw-ish RNA-seq counts
    lib_size = sim_counts$lib_size,
    surv = surv_df,               # time, status
    W_true = W_true,
    H_true = H_true,
    beta_true = beta
  )
}
