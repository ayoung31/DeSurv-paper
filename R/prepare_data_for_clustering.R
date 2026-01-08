# Prepare expression matrix based on selected factors
prepare_data_for_clustering <- function(tops, data, facs, weight) {

  facs <- clean_factor_ids(facs)
  if (!length(facs)) {
    stop("No valid factors supplied after filtering invalid selections.")
  }

  ## pull out key variables
  ex <- data$ex

  ## filter to genes of interest
  genes <- unlist(tops[, facs, drop = FALSE])
  keep_genes <- which(rownames(ex) %in% genes)
  Xtemp <- ex[keep_genes, ]
  if (nrow(Xtemp) == 0) stop("No matching genes found after filtering.")

  ## filter samples with zero variance
  sds <- apply(Xtemp, 2, sd)
  nonzero_samps <- which(sds > 0)
  n_dropped <- sum(sds == 0)

  # Finding 4 fix: Track which samples are kept for proper label alignment downstream
  original_sample_ids <- colnames(Xtemp)
  Xtemp <- Xtemp[, nonzero_samps, drop = FALSE]
  kept_sample_ids <- colnames(Xtemp)

  if (n_dropped > 0) {
    warning(sprintf(
      "%d sample(s) had zero expression variance for all genes in this factor and were excluded from clustering.",
      n_dropped
    ))
  }

  ## calculate sample weights for consensus clustering
  # Finding 4 fix: Only compute weights for kept samples
  if (weight) {
    temp <- data$sampInfo[rownames(data$sampInfo) %in% kept_sample_ids, , drop = FALSE]
    ndataset <- length(unique(temp$dataset))

    weights <- temp %>%
      dplyr::group_by(dataset) %>%
      dplyr::summarise(n = dplyr::n(), weight = 1 / (dplyr::n() * ndataset), .groups = "drop")
    temp <- temp %>% dplyr::left_join(weights, by = "dataset")
    # Ensure weights are in the same order as Xtemp columns
    weightsItem <- temp$weight[match(kept_sample_ids, rownames(temp))]
  } else {
    weightsItem <- NULL
  }

  # Finding 4 fix: Return sample indices for proper alignment
  return(list(
    Xtemp = Xtemp,
    weightsItem = weightsItem,
    kept_samples = nonzero_samps,
    kept_sample_ids = kept_sample_ids,
    n_dropped = n_dropped
  ))
}
