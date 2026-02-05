test_that("simulate_W_marker_background honors overlap and partitions", {
  set.seed(123)
  G <- 60
  K <- 3
  markers_per_factor <- 6
  B_size <- 10
  overlap_frac <- 0.5
  out <- simulate_W_marker_background(
    G = G,
    K = K,
    markers_per_factor = markers_per_factor,
    B_size = B_size,
    marker_overlap = overlap_frac
  )
  expect_equal(dim(out$W), c(G, K))
  expect_length(out$marker_sets, K)
  expect_equal(names(out$marker_sets), colnames(out$W))
  expect_true(all(vapply(out$marker_sets, length, integer(1)) == markers_per_factor))
  expected_overlap <- floor(markers_per_factor * overlap_frac)
  expect_equal(
    length(intersect(out$marker_sets[[1]], out$marker_sets[[2]])),
    expected_overlap
  )
  expect_equal(
    length(intersect(out$marker_sets[[1]], out$marker_sets[[3]])),
    expected_overlap
  )
  unique_markers <- length(unique(unlist(out$marker_sets)))
  expect_length(out$background, B_size)
  expect_length(out$noise_genes, G - B_size - unique_markers)
  norms <- colSums(out$W^2)
  expect_true(all(abs(norms - 1) < 1e-6))
})

test_that("simulate_W_marker_background validates inputs", {
  expect_error(
    simulate_W_marker_background(G = 10, K = 2, markers_per_factor = 0),
    "markers_per_factor"
  )
  expect_error(
    simulate_W_marker_background(G = 10, K = 2, markers_per_factor = 5, B_size = 10),
    "Not enough genes"
  )
  expect_error(
    simulate_W_marker_background(
      G = 20,
      K = 2,
      markers_per_factor = 3,
      B_size = 4,
      marker_overlap = 1.5
    ),
    "marker_overlap"
  )
})

test_that("select_k_by_lcb chooses smallest k within LCB of best observed", {
  k_df <- data.frame(
    k = c(2, 3, 4),
    mean_cindex = c(0.70, 0.75, 0.74),
    pred_mean = c(0.72, 0.73, 0.715),
    pred_sd = c(0.005, 0.01, 0.005),
    stringsAsFactors = FALSE
  )
  res <- select_k_by_lcb(k_df, lcb_level = 0.90)
  expect_equal(res$k_selected, 2L)
  expect_equal(res$k_best, 3L)
  expect_true(is.finite(res$lcb_threshold))
})

test_that("select_k_by_lcb falls back when top prediction is missing", {
  k_df <- data.frame(
    k = c(2, 3),
    mean_cindex = c(0.70, 0.75),
    pred_mean = c(0.71, NA_real_),
    pred_sd = c(0.01, NA_real_),
    stringsAsFactors = FALSE
  )
  res <- select_k_by_lcb(k_df, lcb_level = 0.90)
  expect_equal(res$k_selected, 3L)
  expect_equal(res$reason, "missing_prediction")
})

test_that("select_bo_k_lcb uses best observed when GP is unavailable", {
  history <- data.frame(
    eval_id = 1:4,
    k_grid = c(2, 2, 3, 3),
    mean_cindex = c(0.65, 0.68, 0.70, 0.69),
    status = "ok",
    stringsAsFactors = FALSE
  )
  bo_results <- list(history = history, km_fit = NULL, bounds = data.frame())
  res <- select_bo_k_lcb(bo_results, lcb_level = 0.90)
  expect_equal(res$k_selected, 3L)
  expect_equal(res$reason, "no_gp")
})
