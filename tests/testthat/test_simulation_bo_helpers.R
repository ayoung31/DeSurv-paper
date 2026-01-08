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
