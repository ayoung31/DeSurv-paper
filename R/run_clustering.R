
run_clustering <- function(tops,
                           data, 
                           gene_lists, 
                           color.lists = NULL,
                           facs = NULL,
                           base_dir,
                           WtX = FALSE,
                           reps = 1000, pFeature = 1, pItem = .8, seed = 9999,
                           clusterAlg = "km", distance = "euclidean"){
  
  
  clusCol=NULL
  if(length(facs)>0 ){
    if(WtX){
      dataname = data$dataset
      suffix = "WtX"
    }else{
      dataname = data$dataname
      suffix = "X"
    }
    file_name = paste0("clustering","_",suffix,".pdf")
    dir = file.path(base_dir,dataname)
    # --- Prepare expression data ---
    clus_data <- prepare_data_for_clustering(tops=tops,
                                         data=data,
                                         facs=facs,
                                         WtX = WtX)
    Xtemp <- clus_data$Xtemp
    kept_sample_ids <- clus_data$kept_sample_ids
    
    
    # --- Set defaults if needed ---
    maxK <- length(facs) * 3
    
    # --- Consensus clustering on samples ---
    clusCol <- run_consensus_clustering(
      mat = Xtemp, maxK = maxK, reps = reps, pItem = pItem,
      pFeature = pFeature, seed = seed, clusterAlg = clusterAlg,
      distance = distance, dir = dir, file_name
    )
  }

  
  return(clusCol)
}
