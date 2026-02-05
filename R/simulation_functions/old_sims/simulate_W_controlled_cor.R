simulate_W_controlled_cor <- function(
    G = 5000, K = 4,
    n_marker = 150,
    base_sd = 1,
    mix_strength = 0.3,    # 0 = more orthogonal, >0 = more correlated
    shift = 1.0,           # baseline shift to make entries positive
    marker_boost = 2.0,    # multiply marker genes in their main program
    other_shrink = 0.5,    # multiply those marker genes in other programs
    noise_sd = 0.1
) {
  # Step 1: base structured programs (can be mildly correlated)
  B <- matrix(rnorm(G * K, sd = base_sd), nrow = G, ncol = K)
  qrB <- qr(B)
  Q  <- qr.Q(qrB)[, seq_len(K), drop = FALSE]  # G x K, approx orthonormal
  
  # Mixing to induce some correlation across programs
  A <- diag(K)
  if (mix_strength > 0) {
    off <- matrix(rnorm(K * K, sd = mix_strength), nrow = K, ncol = K)
    diag(off) <- 0
    A <- A + off
  }
  
  W <- Q %*% A
  
  # Shift to nonnegativity and truncate
  W <- W + shift
  W <- pmax(W, 0)
  
  # Step 2: impose marker structure
  marker_indices <- vector("list", K)
  
  for (k in seq_len(K)) {
    # top n_marker genes for program k
    ord <- order(W[, k], decreasing = TRUE)
    mk  <- ord[seq_len(min(n_marker, G))]
    
    marker_indices[[k]] <- mk
    
    # boost in their own program
    W[mk, k] <- W[mk, k] * marker_boost
    
    # shrink in other programs so they look more specific
    others <- setdiff(seq_len(K), k)
    if (length(others) > 0 && other_shrink != 1) {
      W[mk, others] <- W[mk, others] * other_shrink
    }
  }
  
  # Optional small noise + nonnegativity
  if (noise_sd > 0) {
    W <- pmax(W + matrix(rnorm(G * K, sd = noise_sd), nrow = G, ncol = K), 0)
  }
  
  rownames(W) <- paste0("gene", seq_len(G))
  colnames(W) <- paste0("prog", seq_len(K))
  
  # Attach marker info
  marker_info <- lapply(seq_len(K), function(k) {
    idx <- marker_indices[[k]]
    list(
      index = idx,
      genes = rownames(W)[idx]
    )
  })
  names(marker_info) <- colnames(W)
  attr(W, "marker_info") <- marker_info
  
  W
}
