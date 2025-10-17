
# Align W columns across folds and compute a consensus W
# W_list: list of P x k matrices (one per fold)
# beta_list: optional list of length-k numeric vectors (one per fold) to weight matches
# ref: reference fold index (integer) or "auto" to pick the best (highest average pairwise cosine)
# tau: similarity threshold; matches below this are dropped for consensus
# consensus: "mean" or "median" to aggregate W across folds
align_W_folds <- function(
    W_list,
    beta_list = NULL,
    ref = 1,
    tau = 0.2,
    consensus = c("median","mean"),
    rows = c("intersection","union"),
    fill = 0,
    dedup = c("sum","mean")
) {
  consensus <- match.arg(consensus)
  rows      <- match.arg(rows)
  dedup     <- match.arg(dedup)
  
  # 0) Harmonize rows by rowname
  harm <- harmonize_W_by_rownames(W_list, rows = rows, fill = fill, dedup = dedup)
  W_list <- harm$W_list
  P <- length(harm$rows)
  k <- ncol(W_list[[1]])
  F <- length(W_list)
  stopifnot(all(vapply(W_list, ncol, 1L) == k))
  
  nrm_cols <- function(A) sweep(A, 2, sqrt(colSums(A^2)) + 1e-12, "/")
  
  # If ref == "auto", pick fold most similar to others by average greedy cosine
  if (identical(ref, "auto")) {
    Wn_all <- lapply(W_list, nrm_cols)
    avg_cos <- sapply(seq_len(F), function(r) {
      Wr <- Wn_all[[r]]
      sims <- c()
      for (f in setdiff(seq_len(F), r)) {
        S <- t(Wr) %*% Wn_all[[f]]
        sims <- c(sims, apply(S, 1, max))
      }
      mean(sims)
    })
    ref <- which.max(avg_cos)
  }
  
  W_ref    <- W_list[[ref]]
  W_ref_n  <- nrm_cols(W_ref)
  beta_ref <- if (!is.null(beta_list)) beta_list[[ref]] else NULL
  
  perms     <- vector("list", F)
  sims      <- vector("list", F)
  W_aligned <- vector("list", F)
  
  sim_matrix <- function(Wr_n, Wf_n, br = NULL, bf = NULL) {
    S <- t(Wr_n) %*% Wf_n  # k x k cosine
    if (!is.null(br) && !is.null(bf)) {
      wr <- pmin(1, abs(br) / (mad(abs(br)) + 1e-12))
      wf <- pmin(1, abs(bf) / (mad(abs(bf)) + 1e-12))
      S <- (wr %o% wf) * S
    }
    S
  }
  
  for (f in seq_len(F)) {
    Wf <- W_list[[f]]
    if (f == ref) {
      perms[[f]] <- seq_len(k)
      sims[[f]]  <- rep(1, k)
      W_aligned[[f]] <- Wf
      next
    }
    S <- sim_matrix(W_ref_n, nrm_cols(Wf),
                    beta_ref, if (!is.null(beta_list)) beta_list[[f]] else NULL)
    perm <- as.integer(clue::solve_LSAP(1 - S, maximum = FALSE))
    matched_sim <- S[cbind(seq_len(k), perm)]
    
    # Permute
    Wp <- Wf[, perm, drop=FALSE]
    
    # Rescale to match reference columns (LS positive scale)
    num <- colSums(Wp * W_ref)
    den <- colSums(Wp * Wp) + 1e-12
    c   <- pmax(num / den, 0)
    Wm  <- sweep(Wp, 2, c, "*")
    
    perms[[f]]     <- perm
    sims[[f]]      <- matched_sim
    W_aligned[[f]] <- Wm
  }
  
  # Build consensus over rows in the harmonized order
  S_mat    <- do.call(cbind, lapply(sims, function(s) s))  # k x F
  keep_msk <- S_mat >= tau
  
  consensus_W <- matrix(0, nrow = P, ncol = k,
                        dimnames = list(harm$rows, colnames(W_ref)))
  for (j in seq_len(k)) {
    cols_j <- lapply(seq_len(F), function(f) if (keep_msk[j,f]) W_aligned[[f]][, j] else NULL)
    cols_j <- cols_j[!vapply(cols_j, is.null, TRUE)]
    if (length(cols_j) == 0L) {
      consensus_W[, j] <- W_ref[, j]
    } else {
      Wstack <- do.call(cbind, cols_j)
      consensus_W[, j] <- if (consensus == "median") apply(Wstack, 1, median) else rowMeans(Wstack)
    }
  }
  
  list(
    ref_index   = ref,
    rownames    = harm$rows,      # common (or union) row set used
    perms       = perms,
    sims        = sims,
    W_aligned   = W_aligned,
    consensus_W = consensus_W
  )
}
