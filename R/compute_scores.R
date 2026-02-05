compute_scores = function(tops,W,X,y,delta,score_bin=FALSE){
  g_common <- intersect(rownames(W), rownames(X))
  
  Wc <- W[g_common, , drop = FALSE]
  Xc <- X[g_common, , drop = FALSE]
  
  k <- ncol(tops)
  idx_list <- lapply(seq_len(k), function(i) {
    m <- match(tops[, i], g_common)
    m[!is.na(m)]
  })
  rows <- unlist(idx_list, use.names = FALSE)
  cols <- rep.int(seq_len(k), vapply(idx_list, length, integer(1)))
  
  M <- Matrix::sparseMatrix(
    i = rows, j = cols, x = 1,
    dims = dim(Wc), dimnames = dimnames(Wc)
  )
  
  if (score_bin) {
    W_eff <- Matrix::drop0(Matrix::Matrix((Wc > 0) * (M != 0), sparse = TRUE))
  } else {
    W_eff <- Matrix::drop0(Matrix::Matrix(Wc, sparse = TRUE) * M)
  }
  
  XtW = as.matrix(Matrix::crossprod(as.matrix(Xc), W_eff))
  
  score_data <- data.frame(scale(XtW),
                           time = y,
                           event = delta,
                           check.names = FALSE)
  return(score_data)
}


