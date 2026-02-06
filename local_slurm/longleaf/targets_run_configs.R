# targets_run_configs.R - UNC Longleaf HPC mode
# Same as local full mode — ninit_full = 100 is not capped by LOCAL_CPU_LIMIT
# (it affects crew cpus_per_task via DEFAULT_NINIT_FULL, not the config value)

targets_run_config <- function(label) {
  base_full <- list(
    ninit_full = 100,  # Full HPC value
    run_tol = 1e-5,
    run_maxit = 4000,
    std_nmf_k_grid = 2:12,
    std_nmf_nrun = NULL,
    coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    coxnet_alpha_grid = seq(0, 1, by = 0.1)
  )
  base_full
}

targets_run_configs <- function(bo_labels = names(targets_bo_configs())) {
  if (is.null(bo_labels)) {
    bo_labels <- character(0)
  }
  configs <- lapply(bo_labels, targets_run_config)
  names(configs) <- bo_labels
  configs
}
