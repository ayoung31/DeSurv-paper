# targets_bo_configs.R - UNC Longleaf HPC mode
# Matched to student's master branch for exact reproduction.
#
# WARNING: ninit and desurv_ncores_grid are inside the hashed bo_config.
# Changing them invalidates all BO targets. This is expected since Longleaf
# uses a fresh store.

targets_bo_configs <- function() {
  list(
    # TCGA+CPTAC combined training (main analysis)
    tcgacptac = list(
      data_mode = "external",
      data_loader = "load_data",
      train_datasets = c("TCGA_PAAD", "CPTAC"),
      method_trans_train = "rank",
      desurv_bo_bounds = list(
        k_grid = list(lower = 2L, upper = 10L, type = "integer"),
        alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
        lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
        nu_grid = list(lower = 0, upper = 1, type = "continuous")
      ),
      ngene_config = c(500, 5000),
      ntop_config = c(50, 250),
      lambdaw_config = c(0),
      lambdah_config = c(1e-7, 1e2),
      ninit = 50,
      bo_n_init = 20,
      bo_n_iter = 100,
      bo_candidate_pool = 4000,
      bo_max_refinements = 2,
      bo_tol_gain = 0.002,
      bo_plateau = 1,
      bo_top_k = 10,
      bo_shrink_base = 0.5,
      bo_importance_gain = 0.3,
      bo_coarse_control = NULL,
      bo_refine_control = NULL,
      bo_tol = 1e-5,
      bo_maxit = 4000,
      nfold = 5,
      desurv_parallel_grid = TRUE,
      desurv_ncores_grid = 50
    )
  )
}
