preprocess_data_val = function(data,ngene=1000,genes=NULL,method_trans_train="rank"){
  ### restrict to samp_keeps
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  
  
  ### restrict genes
  if(is.null(genes)){
    # if gene list not supplied, restrict to top ngene highly expressed & variable
    D = list()
    datasets=unique(data$sampInfo$dataset)
    for(d in 1:length(datasets)){
      X = data$ex[,data$sampInfo$dataset==datasets[d]]
      if(ngene=="all"){
        D[[d]] = X
      }else{
        D[[d]] = gene_filter(X=X,ngene=ngene*2)
      }
      
    }
    keep_genes <- Reduce(intersect, lapply(D, function(x) rownames(x)))
    X = data$ex[keep_genes,]
    data$ex = X
    data$featInfo =rownames(X)
  }else{
    # otherwise restrict to supplied gene list
    genes = genes[genes %in% rownames(data$ex)]
    data$ex = data$ex[genes,,drop=FALSE]
    data$ex[is.na(data$ex)] = 0
    data$featInfo = rownames(data$ex)
  }
  
  
  if(method_trans_train=='rank'){
    data$ex = apply(data$ex, 2, rank, ties.method = "average")
  }else if(method_trans_train=='quant'){
    data$ex = preprocessCore::normalize.quantiles(data$ex,keep.names = TRUE)
  }else if(method_trans_train=="none"){
    return(data)
  }else{
    stop("This transformation method is not supported")
  }
  
  return(data)
}
