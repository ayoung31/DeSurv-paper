get_top_genes <- function(W,ntop) {
  maxes = apply(W,2,max)
  nc = sum(maxes>0)
  wc = which(maxes>0)
  # W = W[,maxes>0,drop=FALSE]
  
  if(nc){
    submat = W[,wc,drop=FALSE]
    if(nc>1){
      W[,wc] = submat%*%diag(1/apply(submat,2,max))
    }else{
      W[,wc] = submat/apply(submat,2,max)
    }
    
    if (is.null(W) || !is.matrix(W)) stop("W must be a non-null matrix.")
    if (ntop > nrow(W)) stop("ntop exceeds number of genes available.")
    
    top_genes <- vector("list",ncol(W))
    names(top_genes) = paste0("factor", 1:ncol(W))
    top_diffs <- vector("list",ncol(W))
    flag_empty <- FALSE
    
    for (i in 1:ncol(W)) {
      current_col <- W[, i]
      
      if (sum(current_col) > 0) {
        
        other_cols <- W[, -i, drop = FALSE]
        max_other <- if (ncol(W) > 1) apply(other_cols, 1, max) else rep(0, nrow(W))
        
        diff_vector <- current_col - max_other
        top_indices <- order(diff_vector, decreasing = TRUE)[1:ntop]
        
        top_genes[[paste0("factor", i)]] <- rownames(W)[top_indices]
        top_diffs[[paste0("factor", i)]] <- diff_vector[top_indices]
      } 
    }
    
    if (flag_empty) {
      warning("Some factors had zero weights for all genes.")
    }
    
    
    return(top_genes)
  }
  
  return(NULL)
}
