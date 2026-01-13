targets_run_configs <- function() {
  bo_labels <- names(targets_bo_configs())
  base_full <- list(
    ninit_full = 100,
    run_tol = 1e-5,
    run_maxit = 4000,
    std_nmf_k_grid = 2:12,
    std_nmf_nrun = NULL,
    coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    coxnet_alpha_grid = seq(0, 1, by = 0.1)
  )
  light_full <- modifyList(base_full, list(ninit_full = 10))
  configs <- lapply(bo_labels, function(label) base_full)
  names(configs) <- bo_labels
  if ("easy" %in% bo_labels) {
    configs$easy <- light_full
  }
  if ("bladdereasy" %in% bo_labels) {
    configs$bladdereasy <- light_full
  }
  configs
}
