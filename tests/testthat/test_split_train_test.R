test_that("split_train_test partitions indices without overlap", {
  set.seed(42)
  res <- split_train_test(10, train_frac = 0.5)
  expect_length(res$train, 5)
  expect_length(res$test, 5)
  expect_setequal(c(res$train, res$test), seq_len(10))
  expect_equal(intersect(res$train, res$test), integer(0))
})

test_that("split_train_test handles zero or full training fractions", {
  set.seed(1)
  none <- split_train_test(7, train_frac = 0)
  expect_length(none$train, 0)
  expect_setequal(none$test, seq_len(7))

  set.seed(2)
  all <- split_train_test(4, train_frac = 1)
  expect_setequal(all$train, seq_len(4))
  expect_length(all$test, 0)
})
