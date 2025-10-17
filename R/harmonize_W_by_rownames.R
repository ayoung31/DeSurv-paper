harmonize_W_by_rownames <- function(W_list, rows = c("intersection","union"),
                                    fill = 0, dedup = c("sum","mean")) {
  rows  <- match.arg(rows)
  dedup <- match.arg(dedup)
  
  # Deduplicate rows within each W by rowname
  W_list_dd <- lapply(W_list, function(W) {
    rn <- rownames(W)
    if (is.null(rn)) stop("All W must have rownames to align by rowname.")
    spl <- split(seq_len(nrow(W)), rn)
    # aggregate duplicates
    ag <- lapply(spl, function(ix) {
      if (length(ix) == 1) W[ix,,drop=FALSE]
      else {
        if (dedup == "sum") colSums(W[ix,,drop=FALSE,])
        else colMeans(W[ix,,drop=FALSE,])
      }
    })
    M <- do.call(rbind, ag)
    M <- as.matrix(M)
    rownames(M) <- names(spl)
    M
  })
  
  # Determine target row set & order
  if (rows == "intersection") {
    target_rows <- Reduce(intersect, lapply(W_list_dd, rownames))
    if (length(target_rows) == 0)
      stop("No overlapping rownames across folds.")
  } else {
    target_rows <- sort(unique(unlist(lapply(W_list_dd, rownames))))
  }
  
  # Reindex each W to target_rows, fill missing with 'fill'
  W_h <- lapply(W_list_dd, function(W) {
    miss <- setdiff(target_rows, rownames(W))
    if (length(miss)) {
      add <- matrix(fill, nrow = length(miss), ncol = ncol(W),
                    dimnames = list(miss, colnames(W)))
      W <- rbind(W, add)
    }
    W[target_rows, , drop=FALSE]
  })
  
  list(W_list = W_h, rows = target_rows)
}