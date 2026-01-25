# targets_bo_configs.R - Local desktop version
# Reduced iterations for quick local testing with 20 CPUs
#
# Key differences from HPC version:
#   - ninit = 4 (vs 30): Only 4 random initializations during BO CV
#   - bo_n_init = 4 (vs 20): Only 4 initial BO sample points
#   - bo_n_iter = 4 (vs 50): Only 4 BO optimization iterations
#   - bo_max_refinements = 0 (vs 1-2): No refinement rounds
#   - desurv_ncores_grid = 19: Limited to local CPU count

targets_bo_configs <- function() {
  list(
    # TCGA+CPTAC combined training (main analysis)
    tcgacptac = list(
      data_mode = "external",
      data_loader = "load_data",
      train_datasets = c("TCGA_PAAD", "CPTAC"),
      method_trans_train = "none",
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
      ninit = 4,  # Reduced from 30
      bo_n_init = 4,  # Reduced from 20
      bo_n_iter = 4,  # Reduced from 50
      bo_candidate_pool = 200,  # Reduced from 4000
      bo_max_refinements = 0,  # Reduced from 1-2
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
      desurv_parallel_grid = TRUE,
      desurv_ncores_grid = 19  # Local desktop limit (leave 1 CPU for OS)
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
      ninit = 4,
      bo_n_init = 4,
      bo_n_iter = 4,
      bo_candidate_pool = 200,
      bo_max_refinements = 0,
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
      desurv_parallel_grid = TRUE,
      desurv_ncores_grid = 19
    )
  )
}
