select_nclusters <- function(ccp,
                            k_min = 2,
                            k_max = NULL,
                            min_cluster_frac = 0.03,   # penalize clusters smaller than this fraction
                            size_penalty = 2.0,        # stronger = avoid tiny clusters more
                            verbose = TRUE) {
  
  ks_all <- seq_along(ccp)
  ks_all <- ks_all[!vapply(ccp, is.null, logical(1))]
  
  # ConsensusClusterPlus usually stores results starting at k=2 (so index != k)
  # We'll infer k from names if possible, else assume k = index.
  infer_k <- function(i) {
    nm <- names(ccp)[i]
    if (!is.null(nm) && grepl("^\\d+$", nm)) return(as.integer(nm))
    return(i)
  }
  
  ks <- vapply(ks_all, infer_k, integer(1))
  if (is.null(k_max)) k_max <- max(ks)
  keep <- ks >= k_min & ks <= k_max
  ks_all <- ks_all[keep]
  ks <- ks[keep]
  
  score_one_k <- function(res) {
    C <- res$consensusMatrix
    cl <- res$consensusClass
    if (is.null(C) || is.null(cl)) return(NA_real_)
    cl <- as.integer(cl)
    n <- length(cl)
    if (nrow(C) != n || ncol(C) != n) return(NA_real_)
    
    # Use upper triangle only (exclude diagonal)
    ut <- upper.tri(C, diag = FALSE)
    Ci <- C[ut]
    
    same <- outer(cl, cl, "==")[ut]
    if (!any(same) || all(same)) return(NA_real_)
    
    within <- mean(Ci[same])
    between <- mean(Ci[!same])
    base <- within - between
    
    # Penalize fragmentation: clusters smaller than min_cluster_frac
    tab <- table(cl)
    frac <- as.numeric(tab) / n
    small <- pmax(0, (min_cluster_frac - frac) / min_cluster_frac)  # 0 if >= threshold
    penalty <- sum(small^2)
    
    base - size_penalty * penalty
  }
  
  scores <- vapply(ks_all, function(i) score_one_k(ccp[[i]]), numeric(1))
  best_idx <- which.max(scores)
  k_star <- ks[best_idx]
  
  out <- data.frame(k = ks, score = scores)
  out <- out[order(out$k), , drop = FALSE]
  
  if (verbose) {
    message("Selected K = ", k_star)
  }
  
  list(k = k_star, diagnostics = out)
}
