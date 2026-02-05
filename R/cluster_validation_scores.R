# Helpers for validation latent outputs and clustering

write_validation_latent_outputs <- function(latent_list, base_dir) {
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  lapply(latent_list, function(entry) {
    dataset_dir <- file.path(base_dir, entry$dataset)
    dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(entry$Z) && length(entry$Z)) {
      raw_df <- as.data.frame(entry$Z)
      raw_df <- cbind(sample_id = rownames(entry$Z), raw_df)
      utils::write.csv(
        raw_df,
        file = file.path(dataset_dir, "latent_scores_raw.csv"),
        row.names = FALSE
      )
    }

    if (!is.null(entry$Z_scaled) && length(entry$Z_scaled)) {
      scaled_df <- as.data.frame(entry$Z_scaled)
      scaled_df <- cbind(sample_id = rownames(entry$Z_scaled), scaled_df)
      utils::write.csv(
        scaled_df,
        file = file.path(dataset_dir, "latent_scores_scaled.csv"),
        row.names = FALSE
      )
    }

    if (!is.null(entry$risk_score) && length(entry$risk_score)) {
      risk_df <- data.frame(
        sample_id = rownames(entry$Z),
        risk_score = entry$risk_score,
        stringsAsFactors = FALSE
      )
      utils::write.csv(
        risk_df,
        file = file.path(dataset_dir, "risk_scores.csv"),
        row.names = FALSE
      )
    }
    invisible(entry)
  })
}

cluster_validation_latent <- function(latent_entry,
                                      maxK = VAL_CLUSTER_MAXK,
                                      reps = VAL_CLUSTER_REPS,
                                      pItem = VAL_CLUSTER_PITEM,
                                      pFeature = VAL_CLUSTER_PFEATURE,
                                      seed = VAL_CLUSTER_SEED,
                                      clusterAlg = "km",
                                      distance = "euclidean",
                                      dir = NULL) {
  Z <- latent_entry$Z_scaled
  if (is.null(Z) || !length(Z) || nrow(Z) < 2) {
    warning("Skipping clustering for dataset ", latent_entry$dataset,
            ": insufficient samples or latent factors.")
    return(list(dataset = latent_entry$dataset, assignments = tibble::tibble(), clus_res = NULL))
  }
  mat <- t(Z)
  sample_count <- ncol(mat)
  if (sample_count < 3) {
    warning("Skipping clustering for dataset ", latent_entry$dataset,
            ": need at least 3 samples, found ", sample_count, ".")
    return(list(dataset = latent_entry$dataset, assignments = tibble::tibble(), clus_res = NULL))
  }
  maxK_use <- min(maxK, sample_count - 1L)
  if (maxK_use < 2) {
    warning("Skipping clustering for dataset ", latent_entry$dataset,
            ": cannot form clusters with maxK < 2.")
    return(list(dataset = latent_entry$dataset, assignments = tibble::tibble(), clus_res = NULL))
  }
  pItem_use <- min(max(pItem, 0), 1)
  pFeature_use <- min(max(pFeature, 0), 1)

  clus_res <- run_consensus_clustering_v2(
    mat = mat,
    maxK = maxK_use,
    reps = reps,
    pItem = pItem_use,
    pFeature = pFeature_use,
    seed = seed,
    clusterAlg = clusterAlg,
    distance = distance,
    dir = dir
  )

  assignments <- lapply(seq_along(clus_res), function(k) {
    clus_k <- clus_res[[k]]$consensusClass
    if (is.null(clus_k)) {
      return(NULL)
    }
    tibble::tibble(
      dataset = latent_entry$dataset,
      sample_id = names(clus_k),
      k = k,
      cluster = as.integer(clus_k)
    )
  })
  assignments <- assignments[!vapply(assignments, is.null, logical(1))]
  assignments_tbl <- if (length(assignments)) {
    dplyr::bind_rows(assignments)
  } else {
    tibble::tibble()
  }
  list(dataset = latent_entry$dataset, assignments = assignments_tbl, clus_res = clus_res)
}
