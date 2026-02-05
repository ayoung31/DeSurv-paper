compute_test_metrics <- function(sim_truth, test_dataset, fit) {
  stopifnot(!is.null(test_dataset$sampInfo))
  split_ids <- colnames(test_dataset$ex)
  if (is.null(split_ids) || !length(split_ids)) {
    split_ids <- test_dataset$sampInfo$ID
  }
  if (is.null(split_ids) || !length(split_ids)) {
    split_ids <- rownames(test_dataset$sampInfo)
  }
  truth_test <- subset_truth_for_split(sim_truth, split_ids)
  list(
    c_index = metric_test_cindex(truth_test, fit, dataset = test_dataset)$c_index,
    reconstruction = metric_reconstruction_error(sim_truth, fit)$reconstruction_error,
    column_cor = metric_column_correlations(sim_truth, fit),
    effect_correlation = metric_effect_size_correlation(sim_truth, fit)$effect_size_correlation,
    factor_detection = metric_factor_detection(sim_truth, fit, beta_threshold = 0.2),
    gene_set_recovery = metric_gene_set_recovery(sim_truth, fit, beta_threshold = 0.2),
    global_gene_set_recovery = metric_global_gene_set_recovery(sim_truth, fit, beta_threshold = 0.2)
  )
}
