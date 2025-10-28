get_top_genes_abs <- function(W,ntop) {
  maxes = apply(W,2,max)
  W = W[,maxes>0,drop=FALSE]
  
  if(ncol(W)>0){

    if (is.null(W) || !is.matrix(W)) stop("W must be a non-null matrix.")
    if (ntop > nrow(W)) stop("ntop exceeds number of genes available.")
    
    top_genes <- list()
    flag_empty <- FALSE
    
    for (i in seq_len(ncol(W))) {
      current_col <- W[, i]
      
      if (sum(current_col) > 0) {

        top_indices <- order(current_col, decreasing = TRUE)[1:ntop]
        
        top_genes[[paste0("factor", i)]] <- rownames(W)[top_indices]
      } else {
        flag_empty <- TRUE
        top_genes[[paste0("factor", i)]] = NULL
      }
    }
    
    if (flag_empty) {
      warning("Some factors had zero weights for all genes.")
    }
    return(as.data.frame(top_genes))
  }
  
}
