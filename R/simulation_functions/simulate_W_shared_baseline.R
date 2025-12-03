simulate_W_shared_baseline <- function(
    G = 5000,
    K = 4,
    baseline_mean = 5,          # baseline expression for most genes
    baseline_sd   = 1,
    bump_mean     = 3,          # extra expression for program-specific genes
    bump_sd       = 0.5,
    bump_genes_per_prog = 150,
    noise_sd_baseline = 0.2
) {
  # 1. Baseline expression shared across all programs
  baseline <- rnorm(G, mean = baseline_mean, sd = baseline_sd)
  baseline <- pmax(baseline, 0)
  
  # 2. Start W as K copies of the baseline
  W <- matrix(rep(baseline, K), nrow = G, ncol = K)
  marker_indices <- vector("list", K)
  
  # 3. Add small program-specific "bumps"
  for (k in seq_len(K)) {
    bump_genes <- sample(seq_len(G), bump_genes_per_prog)
    marker_indices[[k]] <- bump_genes
    W[bump_genes, k] <- W[bump_genes, k] +
      pmax(rnorm(bump_genes_per_prog, mean = bump_mean, sd = bump_sd), 0)
  }
  
  # 4. Add tiny gene-level noise so columns aren't literally identical
  W <- W + matrix(rnorm(G * K, mean = 0, sd = noise_sd_baseline), nrow = G, ncol = K)
  W <- pmax(W, 0)
  
  rownames(W) <- paste0("gene", seq_len(G))
  colnames(W) <- paste0("prog", seq_len(K))
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
