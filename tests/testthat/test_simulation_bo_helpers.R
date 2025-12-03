test_that("summarize_simulation_config flattens replicate results", {
  skip_if_not_installed("dplyr")
  config_result <- list(
    config_name = "scenario_a",
    config_signature = "sig__a",
    results = list(
      list(
        simulator = "simulate_fn",
        expression_transform = "log2",
        train_fraction = 0.7,
        replicate = 1L,
        run_id = "run_a",
        bo_summary = list(
          best_score = 0.83,
          training_dir = "results/run_a",
          history_path = "results/run_a/history.csv",
          params_best = list(k = 4L, alpha = 0.3)
        )
      ),
      list(
        simulator = "simulate_fn",
        expression_transform = "log2",
        train_fraction = 0.7,
        replicate = 2L,
        run_id = "run_b",
        bo_summary = list(
          best_score = 0.79,
          training_dir = "results/run_b",
          history_path = "results/run_b/history.csv",
          params_best = list(k = 5L, alpha = 0.4)
        )
      )
    )
  )
  tbl <- summarize_simulation_config(config_result)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 2)
  expect_equal(tbl$config_name, rep("scenario_a", 2))
  expect_equal(tbl$replicate, c(1L, 2L))
  expect_equal(tbl$k, c(4L, 5L))
  expect_equal(tbl$alpha, c(0.3, 0.4))
  expect_equal(tbl$best_score, c(0.83, 0.79))
})

test_that("summarize_simulation_config handles empty replicate lists", {
  empty <- list(config_name = "scenario_empty", config_signature = "sig", results = list())
  tbl <- summarize_simulation_config(empty)
  expect_equal(nrow(tbl), 0)
  expect_true(all(c("config_name", "config_signature") %in% names(tbl)))
})
