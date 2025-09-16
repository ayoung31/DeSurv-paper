subset_with_zeros <- function(df, rows, method_trans_train) {
  # figure out which rows exist
  existing <- intersect(rows, rownames(df))
  missing  <- setdiff(rows, rownames(df))
  
  # subset existing rows
  out <- df[existing, , drop = FALSE]
  
  if(method_trans_train=='rank'){
    out = apply(out, 2, rank, ties.method = "average")
  }else if(method_trans_train=='quant'){
    out = preprocessCore::normalize.quantiles(out,keep.names = TRUE)
  }
  
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
  
  ### restrict to samp_keeps
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  
  
  ### restrict genes
  if(is.null(genes)){
    # if gene list not supplied, restrict to top ngene highly expressed & variable
    
    if(method_trans_train=='rank'){
      data$ex = apply(data$ex, 2, rank, ties.method = "average")
    }else if(method_trans_train=='quant'){
      data$ex = preprocessCore::normalize.quantiles(data$ex,keep.names = TRUE)
    }
 
    data$ex = gene_filter(X=data$ex,ngene=ngene)
    data$featInfo =rownames(data$ex)
    
  }else{
    # otherwise restrict to supplied gene list
    data$ex = subset_with_zeros(data$ex,genes,method_trans_train)
    data$featInfo = rownames(data$ex)
  }
  

  

  
  return(data)
}
