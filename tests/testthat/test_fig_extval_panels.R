testthat::test_that("build_fig_extval_panels handles named validation lists", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("pheatmap")
  testthat::skip_if_not_installed("cowplot")
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("survminer")

  genes <- c("gene1", "gene2")
  ex <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(genes, c("S1", "S2", "S3"))
  )
  sampInfo <- data.frame(
    time = c(1, 2, 3),
    event = c(1, 0, 1),
    PurIST = c("Basal-like", "Classical", "Basal-like"),
    stringsAsFactors = FALSE
  )
  rownames(sampInfo) <- c("S1", "S2", "S3")

  data_val_filtered <- list(
    DatasetA = list(
      ex = ex,
      sampInfo = sampInfo
    )
  )

  fit_desurv <- list(
    beta = c(1, -1),
    W = matrix(
      c(0.5, 0.2, 0.1, 0.3),
      nrow = 2,
      dimnames = list(genes, c("F1", "F2"))
    )
  )
  tops_desurv <- list(
    top_genes = data.frame(
      F1 = genes,
      F2 = rev(genes),
      stringsAsFactors = FALSE
    )
  )

  panels <- build_fig_extval_panels(
    data_val_filtered = data_val_filtered,
    fit_desurv = fit_desurv,
    tops_desurv = tops_desurv
  )
  testthat::expect_true(all(c("A", "B", "C") %in% names(panels)))
})
