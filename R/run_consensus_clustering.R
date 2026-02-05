# Run consensus clustering on samples or genes
run_consensus_clustering <- function(mat, maxK, reps, pItem, pFeature, seed,
                                     clusterAlg, distance,
                                     dir = NULL, file_name) {

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  pdf_file <- file.path(dir, file_name)

  mat <- as.matrix(mat)

  # Compute distance matrix based on the distance argument (previously hardcoded to pearson)
  dmat <- switch(
    distance,
    pearson = as.dist(1 - cor(mat, method = "pearson")),
    spearman = as.dist(1 - cor(mat, method = "spearman")),
    euclidean = dist(t(mat), method = "euclidean"),
    maximum = dist(t(mat), method = "maximum"),
    manhattan = dist(t(mat), method = "manhattan"),
    canberra = dist(t(mat), method = "canberra"),
    binary = dist(t(mat), method = "binary"),
    minkowski = dist(t(mat), method = "minkowski"),
    # Default: correlation-based distance
    as.dist(1 - cor(mat, method = "pearson"))
  )

  grDevices::pdf(file = pdf_file)
  close_pdf <- TRUE
  on.exit({
    if (close_pdf && grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  }, add = TRUE)


  result <- ConsensusClusterPlus::ConsensusClusterPlus(
    dmat, maxK,
    reps = reps, pItem = pItem, pFeature = pFeature,
    seed = seed, clusterAlg = clusterAlg, distance = distance
  )

  
  if (close_pdf) {
    grDevices::dev.off()
    close_pdf <- FALSE
  }
  

  
  return(result)
}
