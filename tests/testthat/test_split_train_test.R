test_that("split_train_validation partitions data without overlap", {
  skip_if_not_installed("caret")
  ex <- matrix(rnorm(40), nrow = 5)
  colnames(ex) <- paste0("S", seq_len(ncol(ex)))
  samp_info <- data.frame(
    ID = colnames(ex),
    dataset = rep(c("a", "b"), each = 4),
    time = seq_len(ncol(ex)),
    event = rep(c(0, 1), each = 2, length.out = ncol(ex)),
    stringsAsFactors = FALSE
  )
  rownames(samp_info) <- samp_info$ID
  data <- list(
    ex = ex,
    sampInfo = samp_info,
    samp_keeps = seq_len(ncol(ex)),
    dataname = "unit_test"
  )
  res <- split_train_validation(data, train_frac = 0.5, seed = 42)
  train_ids <- colnames(res$train$ex)
  test_ids <- colnames(res$test$ex)
  expect_equal(length(intersect(train_ids, test_ids)), 0)
  expect_setequal(c(train_ids, test_ids), colnames(ex))
  expect_equal(
    nrow(res$train$sampInfo) + nrow(res$test$sampInfo),
    nrow(samp_info)
  )
  expect_equal(res$train$samp_keeps, seq_len(nrow(res$train$sampInfo)))
  expect_equal(res$test$samp_keeps, seq_len(nrow(res$test$sampInfo)))
})
