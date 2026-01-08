skip_if_not_installed("DeSurv")
skip_if_not_installed("survival")

run_simulation_pipeline <- function(scenario = "R0",
                                    sim_overrides = list(),
                                    nmf_args = list()) {
  suppressPackageStartupMessages(library(DeSurv))
  default_sim <- list(
    G = 80,
    N = 20,
    markers_per_factor = 6,
    B_size = 12,
    noise_sd = 0.2
  )
  args <- utils::modifyList(default_sim, sim_overrides)
  args <- c(list(scenario = scenario, seed = 123), args)
  sim <- do.call(simulate_desurv_scenario, args)
  expr <- sim$X
  if (is.null(colnames(expr))) {
    colnames(expr) <- paste0("Sample", seq_len(ncol(expr)))
  }
  if (is.null(rownames(expr))) {
    rownames(expr) <- paste0("Gene", seq_len(nrow(expr)))
  }
  time <- pmax(sim$time, 1e-3)
  samp_info <- data.frame(
    ID = colnames(expr),
    dataset = scenario,
    time = time,
    event = sim$status,
    stringsAsFactors = FALSE
  )
  default_nmf <- list(
    alpha = 0.2,
    lambda = 0.1,
    nu = 0.9,
    lambdaW = 1e-3,
    lambdaH = 1e-3,
    tol = 1e-4,
    imaxit = 15,
    ninit = 1,
    maxit = 80,
    verbose = FALSE
  )
  nmf_params <- utils::modifyList(default_nmf, nmf_args)
  fit <- suppressWarnings(do.call(
    DeSurv::desurv_fit,
    c(list(
      X = as.matrix(expr),
      y = samp_info$time,
      d = samp_info$event,
      k = ncol(sim$W)
    ), nmf_params)
  ))
  list(sim = sim, fit = fit)
}

test_that("pipeline fits DeSurv on current simulation output", {
  pipeline <- run_simulation_pipeline(
    scenario = "R0",
    sim_overrides = list(G = 60, N = 18, markers_per_factor = 5, B_size = 10)
  )
  expect_equal(dim(pipeline$fit$W), dim(pipeline$sim$W))
  expect_equal(length(pipeline$fit$beta), ncol(pipeline$sim$W))
  expect_true(all(is.finite(pipeline$fit$beta)))
})
