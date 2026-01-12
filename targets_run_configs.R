targets_run_configs <- function() {
  list(
    easy = list(
      bo_key = c("easy"),#,"bladdereasy","pdacspliteasy"
      ninit_full = 10,
      run_tol = 1e-5,
      run_maxit = 4000,
      std_nmf_k_grid = 2:12,
      std_nmf_nrun = NULL,
      coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
      coxnet_alpha_grid = seq(0, 1, by = 0.1)
    )
    # full = list(
    #   bo_key = names(targets_bo_configs()),
    #   ninit_full = 100,
    #   run_tol = 1e-5,
    #   run_maxit = 4000,
    #   std_nmf_k_grid = 2:12,
    #   std_nmf_nrun = NULL,
    #   coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    #   coxnet_alpha_grid = seq(0, 1, by = 0.1)
    # )
    # bladdereasy = list(
    #   bo_key = "bladdereasy",
    #   ninit_full = 10,
    #   run_tol = 1e-5,
    #   run_maxit = 4000,
    #   std_nmf_k_grid = 2:12,
    #   std_nmf_nrun = NULL,
    #   coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    #   coxnet_alpha_grid = seq(0, 1, by = 0.1)
    # )
    # default = list(
    #   bo_key = "default",
    #   ninit_full = 100,
    #   run_tol = 1e-5,
    #   run_maxit = 4000,
    #   std_nmf_k_grid = 2:12,
    #   std_nmf_nrun = NULL,
    #   coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    #   coxnet_alpha_grid = seq(0, 1, by = 0.1)
    # ),
    # bladder_default = list(
    #   bo_key = "bladder",
    #   ninit_full = 100,
    #   run_tol = 1e-5,
    #   run_maxit = 4000,
    #   std_nmf_k_grid = 2:12,
    #   std_nmf_nrun = NULL,
    #   coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    #   coxnet_alpha_grid = seq(0, 1, by = 0.1)
    # )
  )
}
