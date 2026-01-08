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

  ## filter samples with zero SD across all selected genes
  sds <- apply(Xtemp, 2, sd)
  nonzero_samps <- which(sds > 0)
  if (any(sds == 0)) {
    warning("Some subjects had zero variance across all genes in this factor and were excluded.")
  }
  Xtemp <- Xtemp[, nonzero_samps, drop = FALSE]

  ## track which sample IDs were kept (for proper label alignment downstream)
  kept_sample_ids <- colnames(Xtemp)

  ## calculate sample weights for consensus clustering (AFTER filtering)
  if (weight) {
    ## subset sampInfo to only the kept samples
    temp <- data$sampInfo[kept_sample_ids, , drop = FALSE]
    ndataset <- length(unique(temp$dataset))

    weights <- temp %>%
      dplyr::group_by(dataset) %>%
      dplyr::summarise(n = dplyr::n(), weight = 1 / (dplyr::n() * ndataset)) %>%
      dplyr::ungroup()
    temp <- temp %>% dplyr::left_join(weights, by = "dataset")
    weightsItem <- temp$weight
  } else {
    weightsItem <- NULL
  }

  return(list(
    Xtemp = Xtemp,
    weightsItem = weightsItem,
    kept_sample_ids = kept_sample_ids
  ))
}
