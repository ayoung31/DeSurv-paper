# Run consensus clustering on samples or genes
run_consensus_clustering <- function(mat, maxK, reps, pItem, pFeature, seed,
                                     clusterAlg, distance, cluster_by = "col",
                                     weightsItem = NULL) {
  if (cluster_by == "row") {
    mat <- t(as.matrix(mat))
  } else {
    mat <- as.matrix(mat)
  }
  
  dmat <- as.dist(1 - cor(mat, method = "pearson"))
  
  pdf(file = NULL)
  if(cluster_by=="row"){
    result <- ConsensusClusterPlus::ConsensusClusterPlus(
      dmat, maxK,
      reps = reps, pItem = pItem, pFeature = pFeature,
      seed = seed, clusterAlg = clusterAlg, distance = distance, weightsFeature = weightsItem
    )
  }else{
    result <- ConsensusClusterPlus::ConsensusClusterPlus(
      dmat, maxK,
      reps = reps, pItem = pItem, pFeature = pFeature,
      seed = seed, clusterAlg = clusterAlg, distance = distance, weightsItem = weightsItem
    )
  }
  
  dev.off()
  
  return(result)
}
