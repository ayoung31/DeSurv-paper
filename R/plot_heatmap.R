plot_heatmap <- function(Xtemp, tops, data, clusCol, clusRow = NULL, cluster_name, factors,
                         save = FALSE, legend = NULL, fill = NULL, row.colors = NULL) {
  
  # Set up colors
  nclass <- length(unique(clusCol$consensusClass))
  colors <- determine_colors(unique(clusRow$consensusClass), nclass)
  row.colors <- colors$row.colors
  col.colors <- colors$col.colors
  
  ntop=nrow(tops)
  
  # Build expression dataset
  dataSet <- list(
    ex = as.matrix(Xtemp),
    featInfo = rownames(Xtemp),
    sampID = colnames(Xtemp)
  )
  
  # Build gene color table
  gene_vec <- unlist(tops[, factors])
  if (length(gene_vec) != ntop * length(factors)) {
    warning("Mismatch in number of genes expected from tops.")
  }
  
  gene_colors <- data.frame(
    geneSymbol = names(clusRow$consensusClass),
    Color = row.colors[clusRow$consensusClass],
    stringsAsFactors = FALSE
  )
  
  match_idx <- match(dataSet$featInfo, gene_colors$geneSymbol)
  valid_idx <- which(!is.na(match_idx))
  
  if (length(valid_idx) == 0) {
    stop("No matching genes found between expression matrix and top genes.")
  }
  
  dataSet$ex <- dataSet$ex[valid_idx, , drop = FALSE]
  dataSet$featInfo <- dataSet$featInfo[valid_idx]
  dataSet$genesTemp <- gene_colors[match_idx[valid_idx], ]
  
  # Ensure clustering exists
  if (!cluster_name %in% colnames(data$sampInfo)) {
    stop(paste("Cluster label", cluster_name, "not found in data$sampInfo"))
  }
  
  # Dendrograms
  dataSet$Colv <- as.dendrogram(clusCol$consensusTree)
  dataSet$Rowv <- if (!is.null(clusRow)) as.dendrogram(clusRow$consensusTree) else NULL
  
  # Sample subtype
  dataSet$subtype <- data$sampInfo[[cluster_name]]
  
  # Assign color to each sample
  ColSideColors <- col.colors[dataSet$subtype]

  Plot_heatmap_CC(dataSet, NULL, ColSideColors, cluster_name)
  
  invisible(NULL)
}

