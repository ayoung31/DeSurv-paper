make_truth_object <- function() {
  genes <- paste0("gene", 1:4)
  factors <- paste0("prog", 1:2)
  samples <- paste0("sample", 1:3)
  
  W_true <- matrix(
    c(
      1.0, 0.3,
      0.8, 0.5,
      0.2, 1.1,
      0.1, 0.9
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(genes, factors)
  )
  H_true <- matrix(
    c(
      1.0, 0.2, 0.4,
      0.5, 1.0, 0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(factors, samples)
  )
  counts <- W_true %*% H_true
  surv <- data.frame(
    patient = samples,
    time = c(5, 6, 4),
    status = c(1, 0, 1),
    linpred_true = c(0.3, 0.2, 0.4),
    stringsAsFactors = FALSE
  )
  marker_info <- list(
    prog1 = list(index = c(1, 2), genes = genes[1:2]),
    prog2 = list(index = c(3, 4), genes = genes[3:4])
  )
  list(
    W_true = W_true,
    H_true = H_true,
    counts = counts,
    beta_true = c(0.8, -0.4),
    surv = surv,
    marker_info = marker_info,
    config_name = "unit_test"
  )
}

make_fit_object <- function(noise = 0.05, beta = c(0.75, -0.35)) {
  truth <- make_truth_object()
  W <- truth$W_true + noise
  dimnames(W) <- dimnames(truth$W_true)
  list(
    W = W,
    beta = beta
  )
}

test_that("column correlation metric matches factors", {
  truth <- make_truth_object()
  fit <- make_fit_object(noise = 0.02)
  res <- metric_column_correlations(truth, fit)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), ncol(truth$W_true))
  expect_setequal(res$true_factor, seq_len(ncol(truth$W_true)))
  expect_setequal(res$est_factor, seq_len(ncol(fit$W)))
  expect_true(all(abs(res$correlation) > 0.9))
})

test_that("survival subspace distance reports finite angles", {
  truth <- make_truth_object()
  fit <- make_fit_object(noise = 0.05)
  res <- metric_survival_subspace_distance(truth, fit, beta_threshold = 0.1)
  expect_s3_class(res, "tbl_df")
  expect_equal(res$n_true, 2)
  expect_equal(res$n_est, 2)
  expect_true(is.finite(res$mean_angle))
})

test_that("reconstruction error is near zero for perfect fit", {
  truth <- make_truth_object()
  fit <- list(W = truth$W_true)
  res <- metric_reconstruction_error(truth, fit)
  expect_true(res$reconstruction_error < 1e-10)
})

test_that("factor detection identifies true actives", {
  truth <- make_truth_object()
  fit <- list(
    W = truth$W_true,
    beta = c(0.9, -0.2)
  )
  res <- metric_factor_detection(truth, fit, beta_threshold = 0.5)
  expect_equal(res$precision, 1)
  expect_equal(res$recall, 1/1)  # only prog1 active above threshold
  expect_equal(res$true_active, 1)
})

test_that("effect size correlation equals one for identical betas", {
  truth <- make_truth_object()
  fit <- list(W = truth$W_true, beta = truth$beta_true)
  res <- metric_effect_size_correlation(truth, fit)
  expect_equal(res$effect_size_correlation, 1)
  expect_equal(res$n_matched, 2)
})

test_that("factor ROC curve contains monotonically increasing TPR", {
  truth <- make_truth_object()
  fit <- make_fit_object(beta = c(0.9, 0.1))
  res <- metric_factor_roc_curve(truth, fit, beta_threshold = 0.2)
  expect_true(all(diff(res$tpr) >= -1e-8, na.rm = TRUE))
  expect_true(all(res$precision <= 1, na.rm = TRUE))
})

test_that("gene set metrics report perfect recovery for identical loadings", {
  skip_if_not_installed("dplyr")
  truth <- make_truth_object()
  fit <- list(W = truth$W_true, beta = truth$beta_true)
  rec <- metric_gene_set_recovery(truth, fit, top_n = 2, beta_threshold = 0.2)
  expect_true(all(rec$precision == 1))
  expect_true(all(rec$recall == 1))
  global <- metric_global_gene_set_recovery(truth, fit, top_n = 2, beta_threshold = 0.2)
  expect_equal(global$precision, 1)
  expect_equal(global$recall, 1)
  overlap <- metric_gene_set_overlap_matrix(truth, fit, top_n = 2)
  expect_equal(dim(overlap), c(ncol(truth$W_true), ncol(fit$W)))
  expect_equal(unname(diag(overlap)), rep(1, 2))
})

test_that("BO summary metrics compute medians and convergence", {
  skip_if_not_installed("dplyr")
  history <- tibble::tibble(
    iteration = 1:5,
    k = c(4, 5, 4, 6, 5),
    alpha = seq(0.1, 0.5, length.out = 5),
    lambdaH = seq(0.01, 0.05, length.out = 5),
    mean_cindex = c(0.6, 0.62, 0.65, 0.63, 0.66)
  )
  summary_df <- metric_bo_hyperparameter_distribution(history)
  expect_equal(summary_df$k_median, stats::median(history$k))
  expect_equal(summary_df$alpha_iqr, stats::IQR(history$alpha))
  
  perf <- metric_performance_vs_oracle(summary_df = tibble::tibble(
    k = history$k,
    mean_cindex = history$mean_cindex
  ), k_true = 5)
  expect_equal(perf$k_true, 5)
  expect_equal(perf$cindex_selected, max(history$mean_cindex))
  
  conv <- metric_bo_convergence(history)
  expect_true(all(diff(conv$best_so_far) >= 0))
})

test_that("subspace and gene-set stability metrics behave as expected", {
  skip_if_not_installed("dplyr")
  fit_a <- make_fit_object(noise = 0.0)
  fit_b <- make_fit_object(noise = 0.01)
  stability <- metric_subspace_stability(list(fit_a, fit_b), beta_threshold = 0.2)
  expect_equal(stability$comparisons, 1)
  expect_true(stability$mean_subspace_angle >= 0)
  
  truth <- make_truth_object()
  original_predict <- get0("desurv_predict_dataset", envir = globalenv())
  assign(
    "desurv_predict_dataset",
    function(fit, dataset) {
      tibble::tibble(
        sample_id = colnames(dataset$ex),
        risk_score = seq_along(sample_id) + .null_coalesce(fit$id, 0)
      )
    },
    envir = globalenv()
  )
  on.exit({
    if (is.null(original_predict)) {
      rm("desurv_predict_dataset", envir = globalenv())
    } else {
      assign("desurv_predict_dataset", original_predict, envir = globalenv())
    }
  }, add = TRUE)
  fit_a$id <- 0
  fit_b$id <- 1
  fits <- list(fit_a, fit_b)
  risk <- metric_risk_score_stability(fits, truth)
  expect_equal(risk$mean_correlation, 1)
})

test_that("global gene metrics handle overlap matrices", {
  truth <- make_truth_object()
  fit <- make_fit_object(noise = 0.1)
  overlap <- metric_gene_set_overlap_matrix(truth, fit, top_n = 2)
  expect_equal(nrow(overlap), ncol(truth$W_true))
  expect_equal(ncol(overlap), ncol(fit$W))
})
