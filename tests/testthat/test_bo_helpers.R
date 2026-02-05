test_that("select_bo_k_by_cv_se picks smallest k within best_mean - best_se", {
  history <- data.frame(
    eval_id = 1:5,
    k_grid = c(2, 2, 3, 4, 4),
    alpha_grid = c(0.1, 0.2, 0.1, 0.1, 0.2),
    mean_cindex = c(0.75, 0.79, 0.78, 0.80, 0.77),
    status = "ok",
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    eval_id = rep(1:5, each = 2),
    fold = rep(1:2, times = 5),
    val_cindex = c(
      0.74, 0.76,
      0.78, 0.80,
      0.77, 0.79,
      0.79, 0.81,
      0.76, 0.78
    ),
    stringsAsFactors = FALSE
  )
  bo_results <- list(history = history, diagnostics = diagnostics)
  res <- select_bo_k_by_cv_se(bo_results)
  expect_equal(res$k_selected, 2L)
  expect_equal(res$k_best, 4L)
  expect_true(is.finite(res$lcb_threshold))
  expect_equal(res$reason, "lcb")
})

test_that("select_bo_k_by_cv_se falls back to best observed when SE is missing", {
  history <- data.frame(
    eval_id = 1:3,
    k_grid = c(2, 3, 4),
    mean_cindex = c(0.70, 0.72, 0.74),
    status = "ok",
    stringsAsFactors = FALSE
  )
  bo_results <- list(history = history)
  res <- select_bo_k_by_cv_se(bo_results)
  expect_equal(res$k_selected, 4L)
  expect_equal(res$reason, "missing_se")
})
