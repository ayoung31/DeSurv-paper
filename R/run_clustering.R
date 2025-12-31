
#' Perform consensus clustering on selected gene sets
#' @param tops A dataframe of top genes per factor
#' @param data A list with expression matrix and sample metadata
#' @param gene_lists A named list of gene sets
#' @param color.lists Optional color annotations
#' @param type "each", "bas/clas", or integer factor index
#' @param save Logical, save plots
#' @param plot Logical, generate plots
#' @param facs Vector of factor indices to include (optional)
#' @param maxKcol Max K for sample clustering
#' @param maxKrow Max K for gene clustering
#' @param reps Number of clustering iterations
#' @param pFeature Proportion of features (genes) per iteration
#' @param pItem Proportion of items (samples) per iteration
#' @param seed Random seed
#' @param clusterAlg Clustering algorithm
#' @param distance Distance metric
#' @return data with clustering assignments added
run_clustering <- function(tops, data, gene_lists, color.lists = NULL, type = "bas/clas",
                            save = FALSE, plot = TRUE, facs = NULL, dir = NULL,
                            maxKcol = NULL, maxKrow = NULL,
                            reps = 1000, pFeature = 1, pItem = .8, seed = 9999,
                            clusterAlg = "km", distance = "euclidean",weight=TRUE,replace=TRUE) {
  
  if (is.null(facs) || !length(facs)) {
    n_factors <- ncol(tops)
    if (is.null(n_factors) || n_factors < 1) {
      stop("Cannot determine factor indices because `tops` has zero columns.")
    }
    facs <- seq_len(n_factors)
  }
  facs <- clean_factor_ids(facs)
  if (!length(facs)) {
    stop("No valid factor indices supplied to run_clustering.")
  }
  
  # Load dependency
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("ConsensusClusterPlus package is required.")
  }
  
  # Handle 'each' mode by recursive loop over all factors
  if (type == "each") {
    for (i in seq_len(ncol(tops))) {
      data <- run_clustering_internal(tops, data, gene_lists, color.lists, type = i,
                                      facs = facs, dir = dir,
                                      maxKcol = maxKcol, maxKrow = maxKrow,
                                      reps = reps, pFeature = pFeature, pItem = pItem,
                                      seed = seed, clusterAlg = clusterAlg, distance = distance,
                                      save = save, plot = plot, weight=weight, replace=replace)
    }
  }else if(type == "by tissue type"){
    # id tumor factors
    data <- run_clustering_internal(tops, data, gene_lists, color.lists, type = "tumor",
                                    facs = facs, dir = dir,
                                    maxKcol = maxKcol, maxKrow = maxKrow,
                                    reps = reps, pFeature = pFeature, pItem = pItem,
                                    seed = seed, clusterAlg = clusterAlg, distance = distance,
                                    save = save, plot = plot, weight=weight, replace=replace)
    # id stromal factors
    data <- run_clustering_internal(tops, data, gene_lists, color.lists, type = "stroma",
                                    facs = facs, dir = dir,
                                    maxKcol = maxKcol, maxKrow = maxKrow,
                                    reps = reps, pFeature = pFeature, pItem = pItem,
                                    seed = seed, clusterAlg = clusterAlg, distance = distance,
                                    save = save, plot = plot, weight=weight, replace=replace)
  }else{
    data <- run_clustering_internal(tops, data, gene_lists, color.lists, type = type,
                                    facs = facs, dir = dir,
                                    maxKcol = maxKcol, maxKrow = maxKrow,
                                    reps = reps, pFeature = pFeature, pItem = pItem,
                                    seed = seed, clusterAlg = clusterAlg, distance = distance,
                                    save = save, plot = plot, weight=weight, replace=replace)
  }
  
  return(data)
}
