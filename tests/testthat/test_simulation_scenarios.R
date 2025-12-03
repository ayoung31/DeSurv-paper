check_marker_info <- function(info, W) {
  expect_true(is.list(info))
  expect_length(info, ncol(W))
  for (entry in info) {
    expect_true(all(c("index", "genes") %in% names(entry)))
    expect_true(length(entry$index) >= 1)
    expect_equal(entry$genes, rownames(W)[entry$index])
  }
}

check_common_components <- function(sim, G, N, K, expect_lib_size = FALSE) {
  expect_equal(dim(sim$counts), c(G, N))
  if (expect_lib_size) {
    expect_length(sim$lib_size, N)
  }
  expect_equal(dim(sim$W_true), c(G, K))
  expect_equal(dim(sim$H_true), c(K, N))
  expect_equal(nrow(sim$surv), N)
  expect_true(all(c("patient", "time", "status", "linpred_true") %in% names(sim$surv)))
  check_marker_info(sim$marker_info, sim$W_true)
}

test_that("simulate_desurv_data basic structure", {
  set.seed(123)
  G <- 30; N <- 8; K <- 3
  sim <- simulate_desurv_data(
    G = G, N = N, K = K,
    n_marker = 5,
    baseline_hazard = 0.02,
    censor_rate = 0.1,
    correlated_pairs = list(c(1, 2)),
    beta = c(1, -0.4, 0.2)
  )
  check_common_components(sim, G, N, K, expect_lib_size = TRUE)
  expect_type(sim$beta_true, "double")
  expect_length(sim$beta_true, K)
})

test_that("simulate_desurv_data_shared_baseline structure", {
  set.seed(456)
  G <- 28; N <- 7; K <- 4
  sim <- simulate_desurv_data_shared_baseline(
    G = G, N = N, K = K,
    bump_genes_per_prog = 5,
    baseline_hazard = 0.015
  )
  check_common_components(sim, G, N, K, expect_lib_size = FALSE)
  expect_true(is.numeric(sim$beta_true))
  expect_equal(length(sim$beta_true), K)
  expect_true(sim$lethal_prog %in% seq_len(K))
  expect_true(is.list(sim$params))
})

test_that("simulate_desurv_easy structure", {
  set.seed(789)
  G <- 220; N <- 6; K <- 4
  sim <- simulate_desurv_easy(
    G = G, N = N, K = K,
    big_prog = 2,
    lethal_effect = 1.2,
    background_effect = -0.1
  )
  check_common_components(sim, G, N, K, expect_lib_size = FALSE)
  expect_equal(sim$scenario, "easy")
  expect_length(sim$beta_true, K)
  expect_gt(sim$beta_true[2], sim$beta_true[1])
})

test_that("simulate_desurv_interaction structure", {
  set.seed(321)
  G <- 40; N <- 7; K <- 4
  sim <- simulate_desurv_interaction(
    G = G, N = N, K = K,
    interaction_pair = c(1, 4),
    interaction_gamma = 1.1,
    n_marker = 10
  )
  check_common_components(sim, G, N, K, expect_lib_size = FALSE)
  expect_equal(sim$scenario, "interaction")
  expect_true(is.na(sim$beta_true))
})

test_that("simulate_desurv_nuisance structure", {
  set.seed(654)
  G <- 220; N <- 6; K <- 4
  sim <- simulate_desurv_nuisance(
    G = G, N = N, K = K,
    batch_frac = 0.4,
    lethal_effect = 1.0
  )
  check_common_components(sim, G, N, K, expect_lib_size = FALSE)
  expect_equal(sim$scenario, "nuisance")
  expect_length(sim$beta_true, K)
})

test_that("simulate_desurv_subtle structure", {
  set.seed(987)
  G <- 220; N <- 6; K <- 4
  sim <- simulate_desurv_subtle(
    G = G, N = N, K = K,
    subtle_prog = 3,
    subtle_scale = 0.3
  )
  check_common_components(sim, G, N, K, expect_lib_size = FALSE)
  expect_equal(sim$scenario, "subtle")
  expect_length(sim$beta_true, K)
})
