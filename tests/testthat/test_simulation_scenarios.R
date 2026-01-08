check_simulation_structure <- function(sim, G, N, K) {
  expect_equal(dim(sim$X), c(G, N))
  expect_equal(dim(sim$W), c(G, K))
  expect_equal(dim(sim$H), c(N, K))
  expect_equal(dim(sim$Wtilde), c(G, K))
  expect_equal(dim(sim$scores_XtWtilde), c(N, K))
  expect_length(sim$marker_sets, K)
  expect_length(sim$survival_gene_sets, K)
  expect_equal(names(sim$marker_sets), colnames(sim$W))
  expect_equal(names(sim$survival_gene_sets), colnames(sim$W))
  expect_equal(rownames(sim$W), rownames(sim$X))
  expect_equal(length(sim$beta), K)
  expect_equal(length(sim$time), N)
  expect_equal(length(sim$status), N)
  expect_true(all(sim$status %in% c(0, 1)))
}

test_that("simulate_desurv_scenario returns expected structure", {
  set.seed(123)
  G <- 60
  N <- 12
  K <- 3
  sim <- simulate_desurv_scenario(
    scenario = "R0",
    G = G,
    N = N,
    K = K,
    markers_per_factor = 5,
    B_size = 10,
    noise_sd = 0.2,
    seed = 123
  )
  check_simulation_structure(sim, G, N, K)
  expect_equal(sim$scenario, "R0")
  expect_equal(sim$params$G, G)
  expect_equal(sim$params$markers_per_factor, 5)
  unique_markers <- length(unique(unlist(sim$marker_sets)))
  expect_equal(unique_markers, K * 5)
  expect_length(sim$noise_genes, G - 10 - unique_markers)
  expect_equal(sim$survival_gene_sets[[1]], sim$marker_sets[[1]])
})

test_that("simulate_desurv_scenario honors survival gene overrides", {
  set.seed(456)
  G <- 50
  N <- 10
  K <- 2
  sim <- simulate_desurv_scenario(
    scenario = "R0",
    G = G,
    N = N,
    K = K,
    beta = c(1, 0),
    markers_per_factor = 6,
    B_size = 12,
    survival_gene_n = 2,
    survival_marker_frac = 0.5,
    seed = 456
  )
  check_simulation_structure(sim, G, N, K)
  expect_equal(length(sim$survival_gene_sets[[1]]), 2)
  expect_equal(
    length(intersect(sim$survival_gene_sets[[1]], sim$marker_sets[[1]])),
    1
  )
  expect_equal(
    length(intersect(sim$survival_gene_sets[[1]], sim$background)),
    1
  )
})

test_that("get_desurv_defaults exposes scenario configuration", {
  defaults <- get_desurv_defaults("R_mixed")
  expect_true(is.list(defaults))
  expect_true(all(c("G", "N", "K", "beta") %in% names(defaults)))
  expect_true(!is.null(defaults$survival_gene_n))
  expect_true(!is.null(defaults$survival_marker_frac))
})
