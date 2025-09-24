# Prepare expression matrix based on selected factors
prepare_data_for_clustering <- function(tops, data, facs, weight) {
  
  ## pull out key variables
  ex <- data$ex
  
  ## filter to genes of interest
  genes <- unlist(tops[, facs])
  keep_genes <- which(rownames(ex) %in% genes)
  Xtemp <- ex[keep_genes, ]
  if (nrow(Xtemp) == 0) stop("No matching genes found after filtering.")
  
  ## filter samples
  sds = apply(Xtemp,2,sd)
  nonzero_samps = which(sds>0)
  Xtemp = Xtemp[,nonzero_samps]
  if(any(sds==0)){
    warning("some subjects had all zero expression for all genes in this factor")
  }
  
  ## calculate sample weights for consensus clustering
  if(weight){
    temp=data$sampInfo
    ndataset = length(unique(temp$dataset))
    
    weights = temp %>% group_by(dataset) %>%
      summarise(n=n(),weight = 1/(n()*ndataset)) %>% ungroup()
    temp = temp %>% left_join(weights)
    weightsItem = temp$weight
  }else{
    weightsItem=NULL
  }
  
  
  return(list(Xtemp=Xtemp,weightsItem=weightsItem))
}
