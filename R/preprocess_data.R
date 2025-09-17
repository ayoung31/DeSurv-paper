preprocess_data = function(data,ngene=1000,genes=NULL,method_trans_train="rank"){
  ### restrict to samp_keeps
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  
  ### restrict genes
  if(is.null(genes)){
    # if gene list not supplied, restrict to top ngene highly expressed & variable
    X = gene_filter(X=data$ex,ngene=ngene)
    data$ex = X
    data$featInfo =rownames(X)
  }else{
    # otherwise restrict to supplied gene list
    data$ex = data$ex[genes,,drop=FALSE]
    data$ex[is.na(data$ex)] = 0
    data$featInfo = genes
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
