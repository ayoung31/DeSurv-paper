subset_with_zeros <- function(df, rows) {
  # figure out which rows exist
  existing <- intersect(rows, rownames(df))
  missing  <- setdiff(rows, rownames(df))
  
  # subset existing rows
  out <- df[existing, , drop = FALSE]
  
  # create zero rows for missing
  if (length(missing) > 0) {
    zero_mat <- matrix(0, nrow = length(missing), ncol = ncol(df))
    colnames(zero_mat) <- colnames(df)
    rownames(zero_mat) <- missing
    out <- rbind(out, zero_mat)
  }
  
  # put back in requested order
  out[rows, , drop = FALSE]
}

preprocess_data = function(data,ngene=1000,genes=NULL,method_trans_train="rank",
                           quantile_targets = NULL){
  
  ### restrict genes
  if(is.null(genes)){
    # if gene list not supplied, restrict to top ngene highly expressed & variable
    X = gene_filter(X=data$ex,ngene=ngene)
    data$ex = X
    data$featInfo =rownames(X)
  }else{
    # otherwise restrict to supplied gene list
    data$ex = subset_with_zeros(data$ex,genes)
    data$featInfo = rownames(data$ex)
  }
  
  ### restrict to samp_keeps
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  
  if(method_trans_train=='rank'){
    data$ex = apply(data$ex, 2, rank, ties.method = "average")
  }else if(method_trans_train=='quant'){
    data$ex = preprocessCore::normalize.quantiles(data$ex,keep.names = TRUE)
  }else{
    stop("This transformation method is not supported")
  }
  

  
  return(data)
}
