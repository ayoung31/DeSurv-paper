# Add new entries to run multiple configurations side-by-side.
targets_bo_configs <- function() {
  list(
    default = list(
      data_mode = "external",
      data_loader = "load_data",
      train_datasets = c("TCGA_PAAD"),
      method_trans_train = "rank",
      desurv_bo_bounds = list(
        k_grid = list(lower = 2L, upper = 10L, type = "integer"),
        alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
        lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
        nu_grid = list(lower = 0, upper = 1, type = "continuous")
      ),
      ngene_config = c(500, 5000),
      lambdaw_config = c(0),
      lambdah_config = c(1e-7, 1e2),
      ninit = 50,
      bo_n_init = 15,
      bo_n_iter = 60,
      bo_candidate_pool = 4000,
      bo_max_refinements = 2,
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
      desurv_ncores_grid = NULL
    ),
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
      ngene_config = c(500, 5000),
      lambdaw_config = c(0),
      lambdah_config = c(0),
      ninit = 50,
      bo_n_init = 15,
      bo_n_iter = 60,
      bo_candidate_pool = 4000,
      bo_max_refinements = 2,
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
      desurv_ncores_grid = NULL
    )
  )
}

targets_run_configs <- function() {
  list(
    default = list(
      bo_key = "default",
      ninit_full = 100,
      ntop_config = c(50, 250),
      run_tol = 1e-5,
      run_maxit = 4000,
      std_nmf_k_grid = 2:12,
      std_nmf_nrun = NULL,
      coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
      coxnet_alpha_grid = seq(0, 1, by = 0.1)
    ),
    bladder_default = list(
      bo_key = "bladder",
      ninit_full = 100,
      ntop_config = c(50, 250),
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
    default = list(
      run_key = "default",
      mode = "external",
      val_datasets = c(
        "CPTAC",
        "Dijk",
        "Moffitt_GEO_array",
        "PACA_AU_array",
        "PACA_AU_seq",
        "Puleo_array"
      ),
      use_train_genes_for_val = FALSE,
      val_cluster_maxk = 6L,
      val_cluster_reps = 1000L,
      val_cluster_pitem = 0.8,
      val_cluster_pfeature = 1,
      val_cluster_seed = 9999L
    ),
    bladder_holdout = list(
      run_key = "bladder_default",
      mode = "train_split",
      val_datasets = character(0),
      use_train_genes_for_val = TRUE,
      val_cluster_maxk = 6L,
      val_cluster_reps = 1000L,
      val_cluster_pitem = 0.8,
      val_cluster_pfeature = 1,
      val_cluster_seed = 9999L
    )
  )
}
