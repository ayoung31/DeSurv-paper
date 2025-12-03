skip_if_not_installed("DeSurv")
skip_if_not_installed("survival")

run_simulation_pipeline <- function(sim_fun,
                                    sim_args = list(),
                                    nmf_args = list()) {
  suppressPackageStartupMessages(library(DeSurv))
  default_sim <- list(
    G = 220,
    N = 18,
    K = 3,
    correlated_pairs = list(c(1, 2))
  )
  args <- utils::modifyList(default_sim, sim_args)
  set.seed(123)
  sim <- do.call(sim_fun, args)
  sim$surv$time <- pmax(sim$surv$time, 1e-3)
  expr <- transform_simulation_expression(sim$counts, "log2")
  split <- split_train_test(ncol(expr), train_frac = 0.7)
  train <- build_simulation_dataset(expr, sim$surv, split$train, "train")
  default_nmf <- list(
    alpha = 0.2,
    lambda = 0.1,
    nu = 0.9,
    lambdaW = 1e-3,
    lambdaH = 1e-3,
    tol = 1e-4,
    imaxit = 20,
    ninit = 2,
    maxit = 120,
    verbose = FALSE
  )
  nmf_params <- utils::modifyList(default_nmf, nmf_args)
  train$sampInfo$time <- pmax(train$sampInfo$time, 1e-3)
  fit <- suppressWarnings(do.call(
    DeSurv::desurv_fit,
    c(list(
      X = as.matrix(train$ex),
      y = train$sampInfo$time,
      d = train$sampInfo$status,
      k = ncol(sim$W_true)
    ), nmf_params)
  ))
  list(sim = sim, fit = fit)
}

test_that("pipeline recovers survival metrics for easy scenario", {
  pipeline <- run_simulation_pipeline(
    simulate_desurv_easy,
    sim_args = list(correlated_pairs = list(c(1, 2)))
  )
  expect_equal(dim(pipeline$fit$W), dim(pipeline$sim$W_true))
  corr <- metric_column_correlations(pipeline$sim, pipeline$fit)
  expect_s3_class(corr, "tbl_df")
  expect_equal(nrow(corr), ncol(pipeline$sim$W_true))
  expect_true(all(is.finite(corr$correlation)))
  recon <- metric_reconstruction_error(pipeline$sim, pipeline$fit)
  expect_true(is.numeric(recon$reconstruction_error))
  cidx <- metric_test_cindex(pipeline$sim, pipeline$fit)
  expect_true(is.numeric(cidx$c_index))
  expect_true(cidx$c_index >= 0)
})

test_that("pipeline handles nuisance scenario with confounders", {
  pipeline <- run_simulation_pipeline(
    simulate_desurv_nuisance,
    sim_args = list(
      subtle_prog = 2,
      correlated_pairs = list(c(1, 3))
    )
  )
  expect_equal(dim(pipeline$fit$W), dim(pipeline$sim$W_true))
  detection <- metric_factor_detection(pipeline$sim, pipeline$fit, beta_threshold = 0.2)
  expect_s3_class(detection, "tbl_df")
  expect_true(all(names(detection) %in% c("precision", "recall", "specificity", "fdr", "true_active", "est_active")))
  cidx <- metric_test_cindex(pipeline$sim, pipeline$fit)
  expect_true(is.numeric(cidx$c_index))
})
