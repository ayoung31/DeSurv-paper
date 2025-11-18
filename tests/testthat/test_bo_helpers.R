test_that("standardize_bo_params normalizes *_grid suffixes", {
  params <- list(
    k_grid = 4,
    alpha_grid = 0.35,
    lambda_grid = 0.01,
    nu_grid = 0.6,
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    n_starts = 20
  )
  cleaned <- standardize_bo_params(params)
  expect_equal(names(cleaned), c("k", "alpha", "lambda", "nu", "lambdaW", "lambdaH", "n_starts"))
  expect_equal(cleaned$k, 4)
  expect_equal(cleaned$alpha, 0.35)
  expect_equal(cleaned$n_starts, 20)
})

test_that("standardize_bo_params handles empty or NULL inputs", {
  expect_equal(standardize_bo_params(NULL), list())
  expect_equal(standardize_bo_params(numeric(0)), list())
})

test_that("maybe_add_numeric_bound adds bounds only when needed", {
  base <- list(existing = list(lower = 0, upper = 1))
  # single value: no change
  out <- maybe_add_numeric_bound(base, c(5), "ngene", type = "integer")
  expect_false("ngene" %in% names(out))
  # range with positive values triggers log scale
  out2 <- maybe_add_numeric_bound(base, c(1e-4, 1e-2), "lambdaW_grid", log_scale = TRUE)
  expect_true("lambdaW_grid" %in% names(out2))
  expect_equal(out2$lambdaW_grid$scale, "log10")
  # range including zero should omit log scale
  out3 <- maybe_add_numeric_bound(base, c(0, 1e-2), "lambdaH_grid", log_scale = TRUE)
  expect_false(isTRUE(out3$lambdaH_grid$scale == "log10"))
})
