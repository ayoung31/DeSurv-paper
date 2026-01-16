targets_run_config <- function(label) {
  base_full <- list(
    ninit_full = 100,
    run_tol = 1e-5,
    run_maxit = 4000,
    std_nmf_k_grid = 2:12,
    std_nmf_nrun = NULL,
    coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    coxnet_alpha_grid = seq(0, 1, by = 0.1)
  )
  if (label %in% c("easy", "bladdereasy")) {
    base_full$ninit_full <- 10
  }
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
