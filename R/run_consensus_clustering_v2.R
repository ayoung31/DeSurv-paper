# Run consensus clustering on samples or genes
run_consensus_clustering_v2 <- function(mat, maxK, reps, pItem, pFeature, seed,
                                        clusterAlg, distance, cluster_by = "col",
                                        weightsItem = NULL, dir = NULL) {
  if (is.null(dir) || length(dir) != 1 || is.na(dir) || !nzchar(dir)) {
    dir <- file.path(tempdir(), "desurv_clustering")
  }
  dir <- normalizePath(path.expand(trimws(as.character(dir))),
                       mustWork = FALSE, winslash = "/")
  if (!nzchar(dir)) {
    stop("Failed to determine an output directory for consensus clustering.")
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir, mustWork = TRUE, winslash = "/")
  
  cluster_label <- if (identical(cluster_by, "row")) "genes" else "subjects"
  pdf_file <- file.path(dir, paste0(cluster_label, ".pdf"))
  ccp_folder_name <- paste0(cluster_label, "_ccp")
  ccp_output_dir <- file.path(dir, ccp_folder_name)
  if (dir.exists(ccp_output_dir)) {
    unlink(ccp_output_dir, recursive = TRUE)
  }
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(dir)
  
  mat <- as.matrix(mat)
  if (cluster_by == "row") {
    mat <- t(mat)
  }
  dmat <- as.dist(1 - cor(mat, method = "pearson"))
  
  grDevices::pdf(file = pdf_file)
  close_pdf <- TRUE
  on.exit({
    if (close_pdf && grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  }, add = TRUE)
  
  if (cluster_by == "row") {
    result <- ConsensusClusterPlus::ConsensusClusterPlus(
      dmat, maxK,
      reps = reps, pItem = pItem, pFeature = pFeature,
      seed = seed, clusterAlg = clusterAlg, distance = distance,
      weightsFeature = weightsItem, plot = "pdf", title = ccp_folder_name
    )
  } else {
    result <- ConsensusClusterPlus::ConsensusClusterPlus(
      dmat, maxK,
      reps = reps, pItem = pItem, pFeature = pFeature,
      seed = seed, clusterAlg = clusterAlg, distance = distance,
      weightsItem = weightsItem, plot = "pdf", title = ccp_folder_name
    )
  }
  
  if (close_pdf) {
    grDevices::dev.off()
    close_pdf <- FALSE
  }
  
  if (dir.exists(ccp_output_dir)) {
    files <- list.files(ccp_output_dir, full.names = TRUE)
    if (length(files)) {
      pref <- paste0(cluster_label, "_")
      for (src in files) {
        dst <- file.path(dir, paste0(pref, basename(src)))
        file.rename(src, dst)
      }
    }
    unlink(ccp_output_dir, recursive = TRUE)
  }
  
  return(result)
}
