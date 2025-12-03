subset_truth_for_split <- function(sim_truth, sample_ids) {
  truth <- sim_truth
  idx <- match(sample_ids, sim_truth$surv$patient)
  truth$surv <- sim_truth$surv[idx, , drop = FALSE]
  truth$counts <- sim_truth$counts[, sample_ids, drop = FALSE]
  if (!is.null(sim_truth$expression)) {
    truth$expression <- sim_truth$expression[, sample_ids, drop = FALSE]
  }
  truth$H_true <- sim_truth$H_true[, idx, drop = FALSE]
  truth
}
