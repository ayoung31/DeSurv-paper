# targets_bo_configs.R - Local desktop FULL mode
# Full HPC-quality iterations, but capped at 19 CPUs
#
# Parameters matching HPC quality:
#   - ninit = 30: Full random initializations during BO CV
#   - bo_n_init = 20: Full initial BO sample points
#   - bo_n_iter = 50: Full BO optimization iterations
#   - bo_max_refinements = 0: No early stopping (matches original)
#   - desurv_ncores_grid = 4  # Reduced: 4 crew workers × 4 cores = 16 total: Capped to local CPU count
#
# Note: This will take significantly longer than quick mode

targets_bo_configs <- function() {
  list(
    # TCGA+CPTAC combined training (main analysis)
    tcgacptac = list(
      data_mode = "external",
      data_loader = "load_data",
      train_datasets = c("TCGA_PAAD", "CPTAC"),
      method_trans_train = "rank",
      desurv_bo_bounds = list(
        k_grid = list(lower = 2L, upper = 15L, type = "integer"),
        alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
        lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
        nu_grid = list(lower = 0, upper = 1, type = "continuous")
      ),
      ngene_config = c(2000),
      ntop_config = c(100, 200),
      lambdaw_config = c(0),
      lambdah_config = c(0),
      ninit = 30,  # Full HPC value
      bo_n_init = 20,  # Full HPC value
      bo_n_iter = 50,  # Full HPC value
      bo_candidate_pool = 4000,  # Full HPC value
      bo_max_refinements = 0,  # No early stopping (matches original)
      bo_tol_gain = 0.002,
      bo_plateau = 1,
      bo_top_k = 10,
      bo_shrink_base = 0.3,
      bo_importance_gain = 0.1,
      bo_coarse_control = NULL,
      bo_refine_control = NULL,
      bo_tol = 1e-5,
      bo_maxit = 4000,
      nfold = 5,
      desurv_parallel_grid = TRUE,  # Safe with crew_controller_local (mirai-based, not forked)
      desurv_ncores_grid = 5  # One core per CV fold (nfold=5); 2 workers × 5 = 10 CPUs max
    ),

    # Bladder cancer analysis
    bladder = list(
      data_mode = "split",
      data_loader = "load_data_bladder_vig",
      train_datasets = c("Bladder"),
      split_raw_files = "data/original/IMVigor210_sampInfo_ex_clinical.rds",
      split_train_frac = 0.7,
      split_seed = 123,
      split_strata_vars = c("event", "dataset"),
      method_trans_train = "none",
      desurv_bo_bounds = list(
        k_grid = list(lower = 2L, upper = 10L, type = "integer"),
        alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
        lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
        nu_grid = list(lower = 0, upper = 1, type = "continuous")
      ),
      ngene_config = 2000,
      ntop_config = c(100, 200),
      lambdaw_config = c(0),
      lambdah_config = c(0),
      ninit = 30,  # Full HPC value
      bo_n_init = 20,  # Full HPC value
      bo_n_iter = 50,  # Full HPC value
      bo_candidate_pool = 4000,  # Full HPC value
      bo_max_refinements = 1,  # Standard refinement
      bo_tol_gain = 0.002,
      bo_plateau = 1,
      bo_top_k = 10,
      bo_shrink_base = 0.3,
      bo_importance_gain = 0.1,
      bo_coarse_control = NULL,
      bo_refine_control = NULL,
      bo_tol = 1e-5,
      bo_maxit = 4000,
      nfold = 5,
      desurv_parallel_grid = TRUE,  # Safe with crew_controller_local (mirai-based, not forked)
      desurv_ncores_grid = 5  # One core per CV fold (nfold=5); 2 workers × 5 = 10 CPUs max
    )
  )
}
