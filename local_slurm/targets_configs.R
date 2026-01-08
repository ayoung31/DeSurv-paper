# targets_configs.R - Local desktop version with easy/quick configurations
# These configs use very few initializations and should run ok locally on a desktop.
#
# Key differences from full configs:
#   - ninit = 2 (vs 30): Only 2 random initializations during BO CV
#   - bo_n_init = 4 (vs 20): Only 4 initial BO sample points
#   - bo_n_iter = 4 (vs 50): Only 4 BO optimization iterations
#   - bo_max_refinements = 0 (vs 1-2): No refinement rounds
#   - ninit_full = 10 (vs 100): Only 10 inits for final model
#   - desurv_ncores_grid = 19: Limited to local CPU count
#
# To use: source this file instead of the original targets_configs.R
# or copy these configs into your local targets_configs.R

targets_bo_configs <- function() {
  list(
    # Easy TCGA config - quick functionality test
    easy = list(
        data_mode = "external",
        data_loader = "load_data",
        train_datasets = c("TCGA_PAAD"),
        method_trans_train = "none",
        desurv_bo_bounds = list(
          k_grid = list(lower = 2L, upper = 8L, type = "integer"),
          alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
          lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
          nu_grid = list(lower = 0, upper = 1, type = "continuous")
        ),
        ngene_config = c(1000),
        ntop_config = c(100),
        lambdaw_config = c(0),
        lambdah_config = c(0),
        ninit = 4,  # 4 inits × 5 folds = 20 jobs, good for 19 CPUs
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
        desurv_ncores_grid = 19  # Local desktop limit
    ),

    # Easy bladder config - quick functionality test
    bladdereasy = list(
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
      ngene_config = 1000,
      lambdaw_config = c(0),
      lambdah_config = c(0),
      ninit = 4,  # 4 inits × 5 folds = 20 jobs, good for 19 CPUs
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
      desurv_ncores_grid = 19  # Local desktop limit
    ),

    # Easy PDAC split config - quick functionality test
    pdacspliteasy = list(
      data_mode = "split",
      data_loader = "load_data_pdac_split",
      train_datasets = c(
        "TCGA_PAAD",
        "CPTAC",
        "Dijk",
        "Moffitt_GEO_array",
        "PACA_AU_array",
        "PACA_AU_seq",
        "Puleo_array"
      ),
      split_raw_files = c(
        "data/original/TCGA_PAAD.rds",
        "data/original/TCGA_PAAD.survival_data.rds",
        "data/original/TCGA_PAAD_subtype.csv",
        "data/original/CPTAC.rds",
        "data/original/CPTAC.survival_data.rds",
        "data/original/CPTAC_subtype.csv",
        "data/original/Dijk.rds",
        "data/original/Dijk.survival_data.rds",
        "data/original/Dijk_subtype.csv",
        "data/original/Moffitt_GEO_array.rds",
        "data/original/Moffitt_GEO_array.survival_data.rds",
        "data/original/Moffitt_GEO_array_subtype.csv",
        "data/original/PACA_AU_array.rds",
        "data/original/PACA_AU_array.survival_data.rds",
        "data/original/PACA_AU_array_subtype.csv",
        "data/original/PACA_AU_seq.rds",
        "data/original/PACA_AU_seq.survival_data.rds",
        "data/original/PACA_AU_seq_subtype.csv",
        "data/original/Puleo_array.rds",
        "data/original/Puleo_array.survival_data.rds",
        "data/original/Puleo_array_subtype.csv"
      ),
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
      ngene_config = 1000,
      ntop_config = c(100),
      lambdaw_config = c(0),
      lambdah_config = c(0),
      ninit = 4,  # 4 inits × 5 folds = 20 jobs, good for 19 CPUs
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
      desurv_ncores_grid = 19  # Local desktop limit
    )
  )
}

targets_run_configs <- function() {
  list(
    # Easy run config - references all easy BO configs
    easy = list(
      bo_key = c("easy", "bladdereasy", "pdacspliteasy"),
      ninit_full = 19,  # Match available CPUs for final model fitting
      run_tol = 1e-5,
      run_maxit = 4000,
      std_nmf_k_grid = 2:12,
      std_nmf_nrun = NULL,
      coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
      coxnet_alpha_grid = seq(0, 1, by = 0.1)
    )
  )
}

targets_val_configs <- function() {
  list(
    # Easy validation config
    easy = list(
      run_key = "easy",
      mode = "external",
      val_datasets = c(
        "CPTAC",
        "Dijk",
        "Moffitt_GEO_array",
        "PACA_AU_array",
        "PACA_AU_seq",
        "Puleo_array"
      )
    )
  )
}
