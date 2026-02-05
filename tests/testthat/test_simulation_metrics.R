build_simulation_inputs <- function() {
  set.seed(1)
  W_out <- simulate_W_marker_background(
    G = 40,
    K = 2,
    markers_per_factor = 5,
    B_size = 6,
    marker_overlap = 0
  )
  H <- simulate_H(N = 6, K = 2, correlated = FALSE, seed = 2)
  X <- simulate_X(W_out$W, H, noise_sd = 0, seed = 3)
  list(W_out = W_out, H = H, X = X)
}

test_that("simulate_survival_from_XtW respects marker-only selection", {
  sim <- build_simulation_inputs()
  set.seed(4)
  surv <- simulate_survival_from_XtW(
    X = sim$X,
    W = sim$W_out$W,
    marker_sets = sim$W_out$marker_sets,
    beta = c(1, 0.4),
    baseline_hazard = 0.05,
    censor_rate = 0.02,
    survival_gene_n = 2,
    survival_marker_frac = 1,
    background_genes = sim$W_out$background
  )
  expect_length(surv$survival_gene_sets, ncol(sim$W_out$W))
  expect_equal(length(surv$survival_gene_sets[[1]]), 2)
  expect_true(all(surv$survival_gene_sets[[1]] %in% sim$W_out$marker_sets[[1]]))
  expect_equal(dim(surv$Wtilde), dim(sim$W_out$W))
  expect_equal(dim(surv$scores), c(ncol(sim$X), ncol(sim$W_out$W)))
  expect_true(all(abs(colMeans(surv$scores)) < 1e-6))
  expect_equal(length(surv$time), ncol(sim$X))
  expect_true(all(surv$status %in% c(0, 1)))
})

test_that("simulate_survival_from_XtW can draw background-only genes", {
  sim <- build_simulation_inputs()
  set.seed(5)
  surv <- simulate_survival_from_XtW(
    X = sim$X,
    W = sim$W_out$W,
    marker_sets = sim$W_out$marker_sets,
    beta = c(0.6, -0.2),
    baseline_hazard = 0.05,
    censor_rate = 0.02,
    survival_gene_n = 3,
    survival_marker_frac = 0,
    background_genes = sim$W_out$background
  )
  expect_equal(length(surv$survival_gene_sets[[1]]), 3)
  expect_true(all(surv$survival_gene_sets[[1]] %in% sim$W_out$background))
  marker_only <- setdiff(sim$W_out$marker_sets[[1]], surv$survival_gene_sets[[1]])
  if (length(marker_only) > 0) {
    expect_true(all(surv$Wtilde[marker_only, 1] == 0))
  }
})

test_that("precision_recall summarizes per-factor accuracy", {
  truth <- list(c("a", "b", "c"), c("d", "e"))
  est <- list(c("a", "c", "f"), c("e"))
  pr <- precision_recall(est, truth)
  expect_equal(dim(pr), c(2, 2))
  expect_equal(as.numeric(pr["precision", 1]), 2 / 3)
  expect_equal(as.numeric(pr["recall", 1]), 2 / 3)
  expect_equal(as.numeric(pr["precision", 2]), 1)
  expect_equal(as.numeric(pr["recall", 2]), 1 / 2)
})
