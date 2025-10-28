
run_clustering_internal_v2 <- function(tops, data, gene_lists, color.lists = NULL, type = "bas/clas",
                                    save = TRUE, plot = TRUE, facs = NULL,
                                    dir = NULL,
                                    maxKcol = NULL, maxKrow = NULL,
                                    reps = 1000, pFeature = 0.8, pItem = 1, seed = 9999,
                                    clusterAlg = "km", distance = "euclidean", 
                                    weight=FALSE, replace=TRUE){
  
  
  # --- Determine genes and clustering name ---
  name = paste0("facs_",paste(facs,collapse="_"))
  

  if(length(facs)>0 ){
    # --- Prepare expression data ---
    clus_data <- prepare_data_for_clustering(tops=tops, 
                                         data=data, 
                                         facs=facs, 
                                         weight=weight)
    Xtemp = clus_data$Xtemp
    weightsItem = clus_data$weightsItem
    
    
    # --- Set defaults if needed ---
    if (is.null(maxKrow)) maxKrow <- length(facs)
    if (is.null(maxKcol)) maxKcol <- length(facs) + 2
    
    # --- Consensus clustering on samples ---
    clusCol <- run_consensus_clustering_v2(
      mat = Xtemp, maxK = maxKcol, reps = reps, pItem = pItem,
      pFeature = pFeature, seed = seed, clusterAlg = clusterAlg,
      distance = distance, cluster_by = "col", weightsItem = weightsItem,
      dir = dir
    )
    
    # --- Optional gene clustering ---
    clusRow <- NULL
    if (maxKrow >= 2) {
      clusRow <- run_consensus_clustering_v2(
        mat = Xtemp, maxK = maxKrow + 1, reps = reps, pItem = pItem,
        pFeature = pFeature, seed = seed, clusterAlg = clusterAlg,
        distance = distance, cluster_by = "row", weightsItem=weightsItem,
        dir = dir
      )
    }
    
    # --- Store clustering results ---
    clus_res <- list(clusRow = clusRow, clusCol = clusCol,factors=facs)
    
    
    # --- Annotate and plot ---
    for (k in 2:length(clus_res$clusCol)) {
      class <- clus_res$clusCol[[k]]$consensusClass
      clus_name <- paste0(name, "_with_", k, "clusters")
      
      data$sampInfo[[clus_name]] <- NA
      data$sampInfo[[clus_name]] <- class
      
      for(k2 in 2:length(clus_res$clusRow)){
        if (plot) {
          plot_heatmap(Xtemp=Xtemp, tops=tops, data=data, 
                       clusCol=clus_res$clusCol[[k]], 
                       clusRow=clus_res$clusRow[[k2]], 
                       cluster_name=clus_name, factors=facs,
                       save=FALSE)
          
        }
      }
      if(plot){
        plot_survival(data=data, factors=facs, cluster_name=clus_name)
      }
    }
  }
  
  return(list(data=data,clus_res=clus_res))
}

