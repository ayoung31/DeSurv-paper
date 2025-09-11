#' Construct a prepared clustering input
#'
#' Wraps the filtered expression matrix and metadata into a validated
#' \code{clustering_input} object for downstream clustering.
#'
#' @param Xtemp Numeric matrix (genes x kept samples).
#' @param samp_keeps Integer vector of kept-sample indices.
#' @param weightsItem Numeric vector of per-sample weights or \code{NULL}.
#' @param facs Integer vector of factor indices used to construct \code{Xtemp}.
#' @param name,title Character scalars identifying this analysis context.
#'
#' @return An object of class \code{clustering_input}.
#' @examples
#' \dontrun{
#' inp <- as_clustering_input(mat$Xtemp, mat$samp_keeps, mat$weightsItem,
#'                            facs = sel$facs, name = sel$name, title = sel$title)
#' }
#' @export
as_clustering_input <- function(Xtemp, samp_keeps, weightsItem, facs, name, title) {
  x <- list(
    Xtemp       = if (is.matrix(Xtemp)) Xtemp else as.matrix(Xtemp),
    samp_keeps  = as.integer(samp_keeps),
    weightsItem = weightsItem,
    facs        = as.integer(facs),
    name        = as.character(name),
    title       = as.character(title)
  )
  class(x) <- "clustering_input"
  validate_clustering_input(x)
}

#' Validate a \code{clustering_input}
#' @param x Object to validate.
#' @return The same object (invisibly) if valid, otherwise an error.
#' @keywords internal
#' @export
validate_clustering_input <- function(x) {
  stopifnot(inherits(x, "clustering_input"))
  if (!is.matrix(x$Xtemp) || !is.numeric(x$Xtemp)) stop("Xtemp must be a numeric matrix.")
  if (!is.integer(x$samp_keeps)) stop("samp_keeps must be integer.")
  if (!is.null(x$weightsItem) && !is.numeric(x$weightsItem)) stop("weightsItem must be numeric or NULL.")
  if (!is.integer(x$facs)) stop("facs must be integer.")
  if (length(x$name) != 1L || length(x$title) != 1L) stop("name/title must be length-1 character.")
  x
}

#' Pretty-print a \code{clustering_input}
#' @param x A \code{clustering_input} object.
#' @param ... Unused.
#' @return \code{invisible(x)}.
#' @export
print.clustering_input <- function(x, ...) {
  cat("<clustering_input>\n")
  cat(" name  :", x$name, "\n")
  cat(" title :", x$title, "\n")
  cat(" facs  :", paste(x$facs, collapse = ", "), "\n")
  cat(" dims  :", paste(dim(x$Xtemp), collapse = " x "), "\n")
  invisible(x)
}

#' Construct a clustering result object
#'
#' Creates a validated S3 object of class \code{clustering_result} holding CCP
#' outputs, assignments, and optional plot holders.
#'
#' @param name,title,facs,samp_keeps,samples,weightsItem,Xtemp_dim,clusCol,clusRow,assignments,plots
#'   Fields to populate the result; see return value for structure.
#' @return A \code{clustering_result}.
#' @export
new_clustering_result <- function(
    name, title, facs,
    samp_keeps, samples, weightsItem,
    Xtemp_dim, clusCol, clusRow,
    assignments, plots
) {
  x <- list(
    name        = name,
    title       = title,
    facs        = as.integer(facs),
    samp_keeps  = as.integer(samp_keeps),
    samples     = samples,
    weightsItem = weightsItem,
    Xtemp_dim   = as.integer(Xtemp_dim),
    clusCol     = clusCol,
    clusRow     = clusRow,
    assignments = assignments,
    plots       = plots
  )
  class(x) <- "clustering_result"
  validate_clustering_result(x)
}

#' Validate a \code{clustering_result}
#' @param x Object to validate.
#' @return The same object (invisibly) if valid, otherwise an error.
#' @keywords internal
#' @export
validate_clustering_result <- function(x) {
  stopifnot(inherits(x, "clustering_result"))
  if (!is.character(x$name)  || length(x$name)  != 1) stop("name must be length-1 character.")
  if (!is.character(x$title) || length(x$title) != 1) stop("title must be length-1 character.")
  if (!is.integer(x$facs)) stop("facs must be integer.")
  if (!is.integer(x$samp_keeps)) stop("samp_keeps must be integer.")
  if (!is.character(x$samples)) stop("samples must be character.")
  if (!is.null(x$weightsItem) && !is.numeric(x$weightsItem)) stop("weightsItem must be numeric or NULL.")
  if (!is.integer(x$Xtemp_dim) || length(x$Xtemp_dim) != 2L) stop("Xtemp_dim must be length-2 integer.")
  if (!is.list(x$clusCol)) stop("clusCol must be a list.")
  if (!is.null(x$clusRow) && !is.list(x$clusRow)) stop("clusRow must be a list or NULL.")
  if (!is.data.frame(x$assignments)) stop("assignments must be a data.frame.")
  need <- c("sample","K","cluster")
  if (!all(need %in% names(x$assignments))) stop("assignments must have columns: sample, K, cluster.")
  if (!is.list(x$plots)) stop("plots must be a list.")
  x
}

#' Pretty-print a \code{clustering_result}
#' @param x A \code{clustering_result} object.
#' @param ... Unused.
#' @return \code{invisible(x)}.
#' @export
print.clustering_result <- function(x, ...) {
  cat("<clustering_result>\n")
  cat(" name   :", x$name, "\n")
  cat(" title  :", x$title, "\n")
  cat(" facs   :", paste(x$facs, collapse = ", "), "\n")
  cat(" dims   :", paste(x$Xtemp_dim, collapse = " x "), "\n")
  cat(" samples:", length(x$samples), "kept\n")
  if (length(x$clusCol) >= 1) cat(" Ks     :", paste0(2L, "..", length(x$clusCol)), "\n")
  invisible(x)
}

#' Summarize a \code{clustering_result}
#'
#' @param object A \code{clustering_result}.
#' @param ... Unused.
#' @return An object of class \code{summary.clustering_result} containing:
#'   \itemize{
#'     \item \code{name}, \code{title}, \code{facs}, \code{dims}
#'     \item \code{ks}: data frame with columns \code{K} and \code{n_clusters_observed}
#'   }
#' @export
summary.clustering_result <- function(object, ...) {
  agg <- aggregate(cluster ~ K, object$assignments, function(z) length(unique(z)))
  names(agg)[2] <- "n_clusters_observed"
  structure(list(
    name  = object$name,
    title = object$title,
    facs  = object$facs,
    dims  = object$Xtemp_dim,
    ks    = agg
  ), class = "summary.clustering_result")
}

#' Coerce a \code{clustering_result} to assignments data.frame
#'
#' @param x A \code{clustering_result}.
#' @param ... Ignored.
#' @param row.names,optional Ignored.
#' @return The \code{assignments} data frame with \code{sample}, \code{K}, \code{cluster}.
#' @export
as.data.frame.clustering_result <- function(x, ..., row.names = NULL, optional = FALSE) {
  x$assignments
}

#' Create a container for multiple clustering results
#'
#' Wraps a list of \code{clustering_result} objects in class \code{clustering_results}
#' for convenient printing and downstream handling (e.g., when clustering each factor).
#'
#' @param x A list of \code{clustering_result} objects.
#' @return A \code{clustering_results} object.
#' @export
new_clustering_results <- function(x) {
  if (!is.list(x) || any(!vapply(x, inherits, logical(1), what = "clustering_result")))
    stop("new_clustering_results: x must be a list of clustering_result objects.")
  class(x) <- c("clustering_results", "list")
  x
}

#' Pretty-print a \code{clustering_results} container
#'
#' @param x A \code{clustering_results} object.
#' @param ... Unused.
#' @return \code{invisible(x)}.
#' @export
print.clustering_results <- function(x, ...) {
  cat("<clustering_results> (", length(x), " runs)\n", sep = "")
  nms <- names(x)
  for (i in seq_along(x)) {
    nm <- if (length(nms) && nzchar(nms[i])) nms[i] else paste0("[[", i, "]]")
    cat(" - ", nm, ": ", x[[i]]$name, " (facs: ", paste(x[[i]]$facs, collapse = ", "), ")\n", sep = "")
  }
  invisible(x)
}


#' Select factors by predefined type or custom indices
#'
#' Resolves which latent factors to use for downstream clustering. You can
#' specify a predefined \code{type} (e.g., "surv", "tumor", "stroma", or a
#' single integer index), or pass \code{factors} directly (with an optional
#' \code{label}). Returns factor indices plus a stable name/title for the run.
#'
#' @param results List with model outputs. Must include at least:
#'   \itemize{
#'     \item \code{model.params$k} (integer total factors)
#'     \item \code{ntop}, \code{top.type}, \code{alpha} (used in naming)
#'     \item \code{labels} (data frame with columns \code{factor}, \code{label},
#'       \code{p.adjust}) for tumor/stroma selection
#'     \item \code{fit_cox$beta} for \code{type="surv"}
#'   }
#' @param type One of \code{"surv"}, \code{"tumor"}, \code{"stroma"},
#'   a single integer factor index, or \code{NULL}. Ignored if \code{factors} is supplied.
#' @param factors Optional integer vector of factor indices. If provided, takes
#'   precedence over \code{type}.
#' @param label Optional character label used when \code{factors} is provided.
#'   If a single factor is given and \code{label} is missing, defaults to \code{"factor"}.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{facs}: integer vector of factor indices (sorted, unique)
#'     \item \code{name}: stable ID string for the run (used in filenames/keys)
#'     \item \code{title}: human-readable description for plotting
#'   }
#'
#' @examples
#' \dontrun{
#' sel <- select_factors(results, type = "tumor")
#' sel <- select_factors(results, factors = c(2, 5), label = "customset")
#' sel <- select_factors(results, type = 3)  # single factor by index
#' }
#' @export
select_factors <- function(results, type = NULL, factors = NULL, label = NULL) {
  .nm <- function(rslt, suffix_label, facs_vec) {
    paste0(
      "k=", rslt$model.params$k,
      "_top", rslt$ntop, rslt$top.type, "_",
      suffix_label,
      if (length(facs_vec)) paste0("_", paste0(facs_vec, collapse = "_")) else "",
      "_alpha=", rslt$alpha
    )
  }
  .ttl <- function(rslt, suffix_label, facs_vec) {
    paste0(
      "k=", rslt$model.params$k,
      " top ", rslt$ntop, rslt$top.type,
      " alpha=", rslt$alpha, " ",
      suffix_label,
      if (length(facs_vec)) paste0(" ", paste0(facs_vec, collapse = ", ")) else ""
    )
  }
  .sanitize <- function(fx, k) {
    if (is.null(fx) || length(fx) == 0) return(integer(0))
    fx <- unique(as.integer(fx))
    fx <- fx[is.finite(fx) & fx >= 1 & fx <= k]
    sort(fx)
  }
  
  k <- results$model.params$k
  labels_df <- results$labels
  
  # Custom factors take precedence
  if (!is.null(factors)) {
    facs <- .sanitize(factors, k)
    if (length(facs) == 0) warning("No valid factor indices in `factors` after sanitization.")
    if (length(facs) == 1) {
      label <- "factor"
    } else if (is.null(label) || !nzchar(label)) {
      label <- "custom"
    }
    return(list(facs = facs, name = .nm(results, label, facs), title = .ttl(results, label, facs)))
  }
  
  if (is.null(type)) stop("Provide `type` (surv/tumor/stroma/integer) or `factors` + optional `label`.")
  if (is.numeric(type)) {
    facs <- .sanitize(type, k)
    return(list(facs = facs, name = .nm(results, "factor", facs), title = .ttl(results, "factor", facs)))
  }
  
  type <- match.arg(tolower(type), c("surv","tumor","stroma"))
  
  if (type == "surv") {
    if (is.null(results$fit_cox$beta)) stop("results$fit_cox$beta required for type='surv'.")
    facs <- .sanitize(which(results$fit_cox$beta != 0), k)
    return(list(facs = facs, name = .nm(results, "survival", facs), title = .ttl(results, "survival", facs)))
  }
  
  if (type == "tumor") {
    if (!all(c("factor","label","p.adjust") %in% names(labels_df)))
      stop("results$labels must have factor, label, p.adjust for type='tumor'.")
    cf <- labels_df[grepl("classical", labels_df$label, TRUE), , drop = FALSE]
    bf <- labels_df[grepl("basal",     labels_df$label, TRUE), , drop = FALSE]
    facs <- integer(0)
    if (nrow(cf)) facs <- c(facs, cf$factor[which.min(cf$p.adjust)])
    if (nrow(bf)) facs <- c(facs, bf$factor[which.min(bf$p.adjust)])
    facs <- .sanitize(facs, k)
    return(list(facs = facs, name = .nm(results, "tumor", facs), title = .ttl(results, "tumor factors", facs)))
  }
  
  # stroma
  key_terms <- c("caf","stroma","PurISS","MS")
  idx <- Reduce(`|`, lapply(key_terms, function(kw) grepl(kw, labels_df$label, TRUE)))
  facs <- .sanitize(which(idx), k)
  list(facs = facs, name = .nm(results, "stroma", facs), title = .ttl(results, "stroma factors", facs))
}

#' Build filtered/weighted expression matrix for selected factors
#'
#' Filters \code{data$ex} (genes x samples) to genes associated with the given
#' factors, removes zero-variance samples, and (optionally) computes per-sample
#' weights to equalize dataset contributions.
#'
#' @param results Model results; must contain \code{tops}, where each column
#'   corresponds to a factor and holds a vector of top genes for that factor.
#' @param data List with:
#'   \itemize{
#'     \item \code{ex}: numeric matrix (genes x samples)
#'     \item \code{sampInfo}: data.frame with at least \code{dataset} column
#'     \item \code{samp_keeps}: integer indices of candidate samples
#'   }
#' @param facs Integer vector of selected factor indices.
#' @param weight Logical; if \code{TRUE}, compute \code{weightsItem} so that each
#'   dataset contributes equally (1 / (n_in_dataset * n_datasets)).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{Xtemp}: numeric matrix (genes x kept samples)
#'     \item \code{samp_keeps}: integer indices of samples retained (non-zero SD)
#'     \item \code{weightsItem}: numeric vector of per-sample weights or \code{NULL}
#'   }
#'
#' @examples
#' \dontrun{
#' mat <- build_expression_matrix(results, data, facs = c(2,5), weight = TRUE)
#' }
#' @export
build_expression_matrix <- function(results, data, facs, weight = FALSE) {
  ex         <- data$ex
  tops       <- results$tops
  samp_keeps <- data$samp_keeps
  
  genes      <- unlist(tops[, facs])
  keep_genes <- which(rownames(ex) %in% genes)
  Xtemp      <- ex[keep_genes, , drop = FALSE]
  if (nrow(Xtemp) == 0) stop("No matching genes found after filtering.")
  
  # filter samples with nonzero SD
  sds            <- apply(Xtemp, 2, sd)
  nonzero_samps  <- which(sds > 0)
  if (any(sds == 0)) warning("Some samples had zero variance for selected genes; removed.")
  samp_keeps     <- intersect(samp_keeps, nonzero_samps)
  Xtemp          <- Xtemp[, samp_keeps, drop = FALSE]
  
  # optional dataset weighting
  if (isTRUE(weight)) {
    temp    <- data$sampInfo[samp_keeps, ]
    ndset   <- length(unique(temp$dataset))
    weights <- dplyr::summarise(dplyr::group_by(temp, dataset),
                                n = dplyr::n(),
                                weight = 1/(n * ndset)) |>
      dplyr::ungroup()
    temp    <- dplyr::left_join(temp, weights, by = "dataset")
    weightsItem <- temp$weight
  } else {
    weightsItem <- NULL
  }
  
  list(Xtemp = Xtemp, samp_keeps = samp_keeps, weightsItem = weightsItem)
}


#' Internal: call ConsensusClusterPlus for items or features
#'
#' Thin wrapper around \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}} that
#' flips the matrix for gene clustering and disables CCP's internal plotting.
#'
#' @param mat Numeric matrix; items in columns for sample clustering.
#' @param maxK,reps,pItem,pFeature,seed,clusterAlg,distance See CCP docs.
#' @param cluster_by \code{"col"} to cluster samples (default) or \code{"row"} to cluster genes.
#' @param weightsItem Optional numeric vector of weights (as items or as features, per CCP).
#'
#' @return A list of CCP results (one per K).
#' @keywords internal
consensusclusterplus_run <- function(
    mat, maxK, reps, pItem, pFeature, seed,
    clusterAlg, distance, cluster_by = c("col","row"),
    weightsItem = NULL
) {
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE))
    stop("ConsensusClusterPlus package is required.")
  cluster_by <- match.arg(cluster_by)
  mat_ccp <- if (cluster_by == "row") t(as.matrix(mat)) else as.matrix(mat)
  
  pdf(file = NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  if (cluster_by == "row") {
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = mat_ccp, maxK = maxK, reps = reps, pItem = pItem, pFeature = pFeature,
      clusterAlg = clusterAlg, distance = distance, seed = seed,
      weightsFeature = weightsItem
    )
  } else {
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = mat_ccp, maxK = maxK, reps = reps, pItem = pItem, pFeature = pFeature,
      clusterAlg = clusterAlg, distance = distance, seed = seed,
      weightsItem = weightsItem
    )
  }
}


#' Run consensus clustering from a prepared input
#'
#' Consumes a \code{clustering_input} and returns a validated
#' \code{clustering_result} with clustering objects, tidy assignments, and
#' placeholders for plots. Does not modify any inputs.
#'
#' @param input A \code{clustering_input} object.
#' @param maxKcol,maxKrow Optional maxima for K (derived from data if \code{NULL}).
#' @param reps,pFeature,pItem,seed,clusterAlg,distance ConsensusClusterPlus controls.
#'
#' @return A \code{clustering_result} object.
#'
#' @examples
#' \dontrun{
#' res <- consensus_cluster(input, maxKcol = 6, maxKrow = 3,
#'                          reps = 1000, pFeature = 0.8, pItem = 1,
#'                          seed = 42, clusterAlg = "km", distance = "euclidean")
#' }
#' @export
consensus_cluster <- function(
    input,                      # clustering_input
    maxKcol, maxKrow,
    reps, pFeature, pItem, seed,
    clusterAlg, distance
) {
  validate_clustering_input(input)
  Xtemp <- input$Xtemp
  facs  <- input$facs
  
  if (is.null(maxKrow)) maxKrow <- max(2L, length(facs))
  if (is.null(maxKcol)) maxKcol <- max(2L, min(length(facs) + 2L, ncol(Xtemp)))
  
  clusCol <- consensusclusterplus_run(
    mat = Xtemp, maxK = maxKcol, reps = reps, pItem = pItem, pFeature = pFeature,
    seed = seed, clusterAlg = clusterAlg, distance = distance, cluster_by = "col",
    weightsItem = input$weightsItem
  )
  clusRow <- NULL
  if (maxKrow >= 2) {
    clusRow <- consensusclusterplus_run(
      mat = Xtemp, maxK = maxKrow, reps = reps, pItem = pItem, pFeature = pFeature,
      seed = seed, clusterAlg = clusterAlg, distance = distance, cluster_by = "row",
      weightsItem = input$weightsItem
    )
  }
  
  kept_samples <- colnames(Xtemp)
  assignments <- do.call(
    rbind,
    lapply(2:length(clusCol), function(k)
      data.frame(sample = kept_samples, K = k,
                 cluster = as.integer(clusCol[[k]]$consensusClass),
                 stringsAsFactors = FALSE))
  )
  
  new_clustering_result(
    name = input$name, title = input$title, facs = input$facs,
    samp_keeps = input$samp_keeps, samples = kept_samples, weightsItem = input$weightsItem,
    Xtemp_dim = dim(Xtemp), clusCol = clusCol, clusRow = clusRow,
    assignments = assignments, plots = list(survival = list(), heatmap = list())
  )
}


#' Cluster a dataset by predefined type or custom factors
#'
#' High-level user function. Selects factors (via \code{\link{select_factors}}),
#' builds a filtered/weighted expression matrix (via \code{\link{build_expression_matrix}}),
#' constructs a \code{\link{clustering_input}} and runs CCP (via \code{\link{consensus_cluster}}).
#' Returns a single \code{clustering_result} or a \code{clustering_results} container
#' (when \code{type="each"} or \code{type="by tissue type"}). Does not mutate \code{data}.
#'
#' @param results,data Inputs required by selection and matrix building.
#' @param type One of:
#'   \itemize{
#'     \item \code{"each"}: run clustering per factor (returns \code{clustering_results})
#'     \item \code{"by tissue type"}: runs tumor and stroma selections
#'     \item \code{"surv"}, \code{"tumor"}, \code{"stroma"}: predefined selections
#'     \item single integer factor index
#'   }
#'   If \code{factors} is provided, it takes precedence over \code{type} (except batch modes).
#' @param factors Optional integer vector of factor indices.
#' @param label Optional label when \code{factors} is supplied. Defaults to \code{"factor"}
#'   when a single factor is provided.
#' @param maxKcol,maxKrow,reps,pFeature,pItem,seed,clusterAlg,distance,weight
#'   Clustering controls. Defaults are defined here (the user-facing entry point).
#'
#' @return A \code{clustering_result} or \code{clustering_results}.
#'
#' @examples
#' \dontrun{
#' res  <- cluster_dataset(results, data, type = "tumor")
#' res2 <- cluster_dataset(results, data, factors = c(2,5), label = "customset")
#' all  <- cluster_dataset(results, data, type = "each")
#' }
#' @export
cluster_dataset <- function(
    results, data,
    type = "each", factors = NULL, label = NULL,
    maxKcol = NULL, maxKrow = NULL,
    reps = 1000, pFeature = 0.8, pItem = 1, seed = 9999,
    clusterAlg = "km", distance = "euclidean",
    weight = FALSE
) {
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("ConsensusClusterPlus package is required.")
  }
  if (is.null(results$tops)) stop("`results$tops` is required to enumerate factors.")
  
  build_input <- function(sel) {
    mat <- build_expression_matrix(results, data, facs = sel$facs, weight = weight)
    as_clustering_input(
      Xtemp = mat$Xtemp,
      samp_keeps = mat$samp_keeps,
      weightsItem = mat$weightsItem,
      facs = sel$facs,
      name = sel$name,
      title = sel$title
    )
  }
  
  run_sel <- function(sel) {
    inp <- build_input(sel)
    consensus_cluster(
      input = inp,
      maxKcol = maxKcol, maxKrow = maxKrow,
      reps = reps, pFeature = pFeature, pItem = pItem, seed = seed,
      clusterAlg = clusterAlg, distance = distance
    )
  }
  
  if (identical(type, "each")) {
    outs <- vector("list", ncol(results$tops))
    names(outs) <- paste0("factor_", seq_len(ncol(results$tops)))
    for (i in seq_len(ncol(results$tops))) {
      sel <- select_factors(results, type = i)
      outs[[i]] <- run_sel(sel)
    }
    outs <- outs[!vapply(outs, is.null, logical(1))]
    return(new_clustering_results(outs))
  }
  
  if (identical(type, "by tissue type")) {
    sel_t <- select_factors(results, type = "tumor")
    sel_s <- select_factors(results, type = "stroma")
    outs <- list(tumor = run_sel(sel_t), stroma = run_sel(sel_s))
    outs <- outs[!vapply(outs, is.null, logical(1))]
    return(new_clustering_results(outs))
  }
  
  sel <- select_factors(results, type = type, factors = factors, label = label)
  run_sel(sel)
}




#' Compare discovered clusters to known biological labels
#'
#' This function compares clustering results (from a factorization model)
#' against known biological or clinical labels (e.g. PurIST, DeCAF) 
#' using contingency tables and chi-square tests.
#'
#' @param data A list-like object containing sample metadata, where
#'   `data$sampInfo` is a data.frame with both known labels and cluster assignments.
#' @param results A list-like object containing model parameters:
#'   - `results$model.params$k`: number of clusters
#'   - `results$ntop`: number of top features used
#'   - `results$top.type`: type of feature selection
#'   - `results$alpha`: model regularization parameter
#'   - `results$model_save_dir`: base directory to save output (if saving)
#' @param type Character string specifying which comparison to run.
#'   One of `"all"`, `"tumor"`, `"stroma"`, or `"surv"`.
#' @param save_csv Logical, whether to save the summary results as a CSV file.
#'   Default is `TRUE`.
#'
#' @return A list with two components:
#'   \item{results}{A named list of results for each comparison, including:
#'     - `known`: known label name
#'     - `cluster`: cluster column name in `data$sampInfo`
#'     - `nclus`: number of clusters
#'     - `alpha`: alpha parameter used
#'     - `table`: raw contingency table
#'     - `table_prop`: row-normalized contingency table
#'     - `chisq`: result of `chisq.test()`}
#'   \item{summary}{A data.frame summarizing the chi-square results for all comparisons:
#'     - known, cluster, nclus, alpha, chisq_stat, df, pvalue, fdr}
#'
#' @export
compare_clusters_to_known_labels <- function(data, results, type = "all", save_csv = TRUE) {
  
  #------------------------------------------------------------
  # 1. Choose known labels and cluster labels based on `type`
  #------------------------------------------------------------
  if(type == "tumor"){
    known_labels <- "PurIST"
    pattern <- paste0("k=", results$model.params$k,
                      "_top", results$ntop, results$top.type, "_tumor",
                      "_alpha=", results$alpha, ".*clusters$")
  } else if(type == "stroma"){
    known_labels <- "DeCAF"
    pattern <- paste0("k=", results$model.params$k,
                      "_top", results$ntop, results$top.type, "_stroma",
                      "_alpha=", results$alpha, ".*clusters$")
  } else if(type == "surv"){
    known_labels <- c("PurIST", "DeCAF")
    pattern <- paste0("k=", results$model.params$k,
                      "_top", results$ntop, results$top.type, "_survival",
                      "_alpha=", results$alpha, ".*clusters$")
  } else {
    known_labels <- c("PurIST", "DeCAF")
    pattern <- paste0("k=", results$model.params$k,
                      "_top", results$ntop, results$top.type, "_factor.*",
                      "_alpha=", results$alpha, ".*clusters$")
  }
  
  cluster_labels <- grep(pattern, colnames(data$sampInfo), value = TRUE)
  
  #------------------------------------------------------------
  # 2. Initialize output containers
  #------------------------------------------------------------
  results_list <- list()
  summary_rows <- list()
  
  #------------------------------------------------------------
  # 3. Loop over all label × cluster combinations
  #------------------------------------------------------------
  if(length(cluster_labels) > 0){
    for(known in known_labels){
      for(cluster in cluster_labels){
        
        # extract metadata
        save_name <- sub(".*_(\\w+)_.*_(\\d+clusters)", "\\1_\\2", cluster)
        nclus <- as.numeric(sub(".*_(\\d+)clusters.*", "\\1", cluster))
        alpha <- sub(".*(alpha=.*)", "\\1", cluster)
        
        # contingency table
        tab <- table(data$sampInfo[[known]], data$sampInfo[[cluster]])
        tab2 <- prop.table(tab, 1)
        
        # chi-square test (suppress warnings about small expected counts)
        tst <- suppressWarnings(chisq.test(tab))
        
        # store detailed result
        key <- paste0(known, "_", save_name)
        results_list[[key]] <- list(
          known = known,
          cluster = cluster,
          nclus = nclus,
          alpha = alpha,
          table = tab,
          table_prop = tab2,
          chisq = tst
        )
        
        # store summary row
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          known = known,
          cluster = cluster,
          nclus = nclus,
          alpha = alpha,
          chisq_stat = as.numeric(tst$statistic),
          df = as.numeric(tst$parameter),
          pvalue = as.numeric(tst$p.value),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  #------------------------------------------------------------
  # 4. Assemble summary table and add FDR correction
  #------------------------------------------------------------
  summary_df <- if(length(summary_rows) > 0){
    do.call(rbind, summary_rows)
  } else {
    data.frame(known = character(),
               cluster = character(),
               nclus = integer(),
               alpha = character(),
               chisq_stat = numeric(),
               df = numeric(),
               pvalue = numeric(),
               stringsAsFactors = FALSE)
  }
  
  if(nrow(summary_df) > 0){
    summary_df$fdr <- p.adjust(summary_df$pvalue, method = "fdr")
  } else {
    summary_df$fdr <- numeric()
  }
  
  #------------------------------------------------------------
  # 5. Optionally save CSV
  #------------------------------------------------------------
  if(save_csv){
    save_dir <- file.path(results$model_save_dir, data$dataname, "clusters_vs_known_labels")
    if(!dir.exists(save_dir)){
      dir.create(save_dir, recursive = TRUE)
    }
    csv_file <- file.path(save_dir, paste0("summary_", type, ".csv"))
    write.csv(summary_df, csv_file, row.names = FALSE)
    message("Summary saved to: ", csv_file)
  }
  
  #------------------------------------------------------------
  # 6. Return both detailed results and summary table
  #------------------------------------------------------------
  return(list(results = results_list, summary = summary_df))
}



comp_class = function(data,pairs=NULL,results_a0,results_best,type='all'){
  library(pheatmap)
  
  if(results_a0$ntop != results_best$ntop){
    stop("Number of top genes for the two models must match")
  }
  ntop=basename(results_a0$model_save_dir)
  save_dir = dirname(dirname(results_a0$model_save_dir))
  mod1 = basename(dirname(results_a0$model_save_dir))
  mod2 = basename(dirname(results_best$model_save_dir))
  save_dir = file.path(save_dir,"comps",paste0(mod1,"_VS_",mod2),ntop,data$dataname)
  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive = TRUE)
  }
  
  if(type=="all") type="factor"
  if(type=="surv") type="survival"
  
  if(type=="tumor" | type=="stroma" | type=="survival"){
    pairs=data.frame(alpha.0_factor.name="",alpha.0_factor.number="",
                     alpha.best_factor.name="",alpha.0_facor.number="")
  }
  
  cols_a.0 <- grep(paste0("k=",results_a0$model.params$k,
                          "_top",results_a0$ntop,results_a0$top.type,"_",type,
                          ".*_alpha=0.*clusters$"),
                   colnames(data$sampInfo),value = TRUE)
  
  cols_a.best <- grep(paste0("k=",results_best$model.params$k,
                             "_top",results_best$ntop,results_best$top.type,"_",type,
                             ".*_alpha=best.*clusters$"),
                      colnames(data$sampInfo),value = TRUE)
  
  
  
  
  if(length(cols_a.0)>0 & length(cols_a.best)>0){
    nclus_a.0 <- unique(as.numeric(sub(".*_([0-9]+)\\clusters.*", "\\1", cols_a.0)))
    nclus_a.best <- unique(as.numeric(sub(".*_([0-9]+)\\clusters.*", "\\1", cols_a.best)))
    
    nclus_common = intersect(nclus_a.0,nclus_a.best)
    
    for(i in 1:nrow(pairs)){
      a.0_num = pairs$alpha.0_factor.number[i]
      a.best_num = pairs$alpha.best_factor.number[i]
      
      a.0_name = pairs$alpha.0_factor.name[i]
      a.best_name = pairs$alpha.best_factor.name[i]
      
      
      for(j in nclus_common){
        
        abest = paste0("k=",results_best$model.params$k,"_top",results_best$ntop,
                       results_best$top.type,"_",type,a.best_num,"_alpha=best","_",j,"clusters")
        a0 = paste0("k=",results_a0$model.params$k,"_top",results_a0$ntop,
                    results_a0$top.type,"_",type,a.0_num,"_alpha=0","_",j,"clusters")
        tab=table(data$sampInfo[[a0]],
                  data$sampInfo[[abest]])
        
        tab2 = prop.table(tab,2)
        
        tst = chisq.test(tab)
        
        
        save_name = file.path(save_dir,paste0("compare.clustering_",j,
                                              "clusters_alpha.0.",type,a.0_num,"_alpha.best.",type,a.best_num,".png"))
        
        title=paste0(j," Clusters for\n","alpha=0, ",type,a.0_num,
                     " ",a.0_name," (rows)\n","alpha=best, ",type,
                     a.best_num," ",a.best_name," (cols)\n","Chi-square test: statistic=",
                     round(tst$statistic,2),", df=",tst$parameter,", pvalue=",format(tst$p.value,scientific = TRUE,digits=3))
        
        
        breaks <- seq(0, 1, length.out = 100)
        colors <- colorRampPalette(c("white", "blue"))(99)
        png(save_name)
        print(pheatmap(tab2, display_numbers = 
                         matrix(paste0(tab,"\ncol pct=",round(tab2*100,1),"%"),nrow=j), 
                       cluster_rows = FALSE, cluster_cols = FALSE, 
                       breaks=breaks,colors=colors,
                       number_color = "black", fontsize_number = 12,
                       main=title))
        dev.off()
        
      }
      
    }
  }else{
    warning("no clusters of this type detected")
  }
  
  
}
