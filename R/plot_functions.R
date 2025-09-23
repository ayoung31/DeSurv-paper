#' Generate an ORA dot plot for a single factor
#'
#' This function generates a dot plot from over-representation analysis (ORA) 
#' results (e.g., from `clusterProfiler`) for a specific factor/component 
#' in the results object.
#'
#' @param results A results object containing:
#'   - `results$model_save_dir`: base directory for saving plots
#'   - `results$model.params$k`: number of factors/components
#'   - `results$alpha`: regularization parameter
#'   - `results$labels$label`: optional sample labels
#'   - an ORA results list (`ora_res_list` or `ora_res_list_dbr`)
#' @param factor Integer. Index of the factor/component to plot.
#' @param msigdbr Logical. If `TRUE`, use `ora_res_list_dbr`; 
#'   otherwise use `ora_res_list`.
#' @param p.adj Numeric. Maximum adjusted p-value cutoff for the color scale.
#'
#' @details
#' The function retrieves ORA results for the specified factor and generates 
#' a `clusterProfiler::dotplot`. If no results are found or the results are 
#' empty, no plot is returned.
#'
#' @return A `ggplot` object representing the dot plot, or `NULL` if the 
#'         ORA results are empty for the requested factor.
#'
#' @examples
#' # Generate dot plot for factor 3
#' # p <- ORA_plots(results, factor = 3, msigdbr = TRUE, p.adj = 0.05)
#' # print(p)
#'
#' @export
ORA_plots <- function(results, factor, msigdbr = FALSE, p.adj = 0.05) {
  
  if(msigdbr){
    list_name="ora_res_list_dbr"
  }else{
    list_name="ora_res_list"
  }
  
  #------------------------------------------------------------
  # 1. Prepare output directory and file names
  #------------------------------------------------------------
  save_dir <- file.path(results$model_save_dir, "dotplots")
  if(!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
  }
  
  suffix <- if(msigdbr) "_msigdbr" else ""
  plot_files <- file.path(save_dir, 
                          paste0("dotplot_factor", 1:results$model.params$k, suffix, ".png"))
  
  #------------------------------------------------------------
  # 2. Check for existing plots (skip if already generated)
  #------------------------------------------------------------
  # Only check for factors that actually have labels
  label_mask <- !is.na(results$labels$label)
  if(length(label_mask) > 0 && all(file.exists(plot_files[label_mask]))) {
    message("Dotplots already exist in ", save_dir)
    return(invisible(plot_files))
  }
  
  #------------------------------------------------------------
  # 3. Extract ORA results list
  #------------------------------------------------------------
  if(!list_name %in% names(results)) {
    stop("List name '", list_name, "' not found in results object.")
  }
  ora_res_list <- results[[list_name]]
  

  ora_res <- ora_res_list[[factor]]
  ora_df  <- as.data.frame(ora_res)
  
  if (!is.null(ora_res) && nrow(ora_df) > 0) {
    
    # dotplot
    p1 <- clusterProfiler::dotplot(
      ora_res,
      title = paste0("alpha=", results$alpha, " factor ", factor)
    ) +
      ggplot2::scale_fill_gradientn(
        name = "adjusted p-values",
        colours = c("red", "white", "blue"),
        limits = c(0, p.adj)
      )
  }

  return(p1)
}

#' Plot a heatmap for clustering inputs or results
#'
#' Generic for producing ggplot-based heatmaps from either a
#' \code{clustering_input} (pre-clustering) or a \code{clustering_result}
#' (post-clustering, with cluster strips).
#'
#' @param x An object to plot (see methods).
#' @param ... Passed to methods.
#' @return A \code{ggplot} object.
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' Heatmap for a \code{clustering_result}
#'
#' Draws an expression heatmap using the clustering structure from a
#' \code{clustering_result}. You can supply \code{Xtemp} directly (recommended),
#' or rebuild it from \code{results}+\code{data}.
#'
#' @param x A \code{clustering_result}.
#' @param k Integer K (>=2) specifying which sample clustering to display.
#' @param results,data Optional; used only if \code{Xtemp} is not provided to rebuild the matrix.
#' @param Xtemp Optional numeric matrix (genes x kept samples) aligned to the run.
#' @param palette Color palette for the expression scale.
#' @param scale_rows Logical; if \code{TRUE}, z-score by gene before plotting.
#'
#' @return A \code{ggplot} heatmap.
#' @export
plot_heatmap.clustering_result <- function(
    x, k,
    results = NULL, data = NULL,   # only needed if you want to rebuild Xtemp
    Xtemp = NULL,
    palette = c("darkblue","blue","white","red","darkred"),
    scale_rows = TRUE
) {
  if (missing(k) || length(k) != 1L || !is.numeric(k) || k < 2) stop("`k` must be a single integer >= 2.")
  if (length(x$clusCol) < k || is.null(x$clusCol[[k]])) {
    avail <- paste(seq_along(x$clusCol), collapse = ", ")
    stop(sprintf("clusCol[[%d]] not available. Available K: %s", k, avail))
  }
  
  if (is.null(Xtemp)) {
    if (is.null(results) || is.null(data)) stop("Provide `Xtemp`, or `results`+`data` to rebuild it.")
    rebuilt <- build_expression_matrix(results, data, facs = x$facs, weight = FALSE)
    Xtemp <- rebuilt$Xtemp
  }
  if (!is.matrix(Xtemp)) Xtemp <- as.matrix(Xtemp)
  common <- intersect(colnames(Xtemp), x$samples)
  if (!length(common)) stop("No overlap between Xtemp and run$samples.")
  Xtemp <- Xtemp[, intersect(x$samples, colnames(Xtemp)), drop = FALSE]
  
  if (isTRUE(scale_rows)) {
    cn <- colnames(Xtemp); Xtemp <- t(apply(Xtemp, 1, scale)); colnames(Xtemp) <- cn
  }
  
  dend_col <- as.dendrogram(x$clusCol[[k]]$consensusTree)
  samp_order <- stats::labels(dend_col)
  samp_order <- intersect(samp_order, colnames(Xtemp))
  if (!length(samp_order)) stop("No overlap between dendrogram labels and Xtemp columns.")
  Xtemp <- Xtemp[, samp_order, drop = FALSE]
  
  gene_anno <- NULL
  if (!is.null(x$clusRow)) {
    k_row <- min(length(x$facs), length(x$clusRow))
    if (!is.null(x$clusRow[[k_row]])) {
      dend_row <- as.dendrogram(x$clusRow[[k_row]]$consensusTree)
      gene_order <- stats::labels(dend_row)
      gcommon <- intersect(gene_order, rownames(Xtemp))
      if (length(gcommon)) {
        Xtemp <- Xtemp[gcommon, , drop = FALSE]
        gcl <- x$clusRow[[k_row]]$consensusClass
        gene_anno <- data.frame(Gene = names(gcl), GeneCluster = factor(as.integer(gcl)))
        gene_anno <- gene_anno[gene_anno$Gene %in% rownames(Xtemp), , drop = FALSE]
      }
    }
  }
  
  scl <- x$clusCol[[k]]$consensusClass
  sample_anno <- data.frame(Sample = names(scl), Cluster = factor(as.integer(scl)))
  sample_anno <- sample_anno[match(colnames(Xtemp), sample_anno$Sample), , drop = FALSE]
  
  df <- tibble::as_tibble(Xtemp, rownames = "Gene") |>
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  df$Sample <- factor(df$Sample, levels = colnames(Xtemp))
  df$Gene   <- factor(df$Gene,   levels = rownames(Xtemp))
  sample_anno$Sample <- factor(sample_anno$Sample, levels = colnames(Xtemp))
  if (!is.null(gene_anno)) gene_anno$Gene <- factor(gene_anno$Gene, levels = rownames(Xtemp))
  
  rng <- max(abs(df$Expression), na.rm = TRUE)
  min_val <- -rng; max_val <- rng
  
  p <- ggplot2::ggplot(df, ggplot2::aes(Sample, Gene, fill = Expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      values = scales::rescale(c(min_val, -3, 0, 3, max_val)),
      limits = c(min_val, max_val)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = sample_anno,
      ggplot2::aes(x = Sample, y = 0, fill = Cluster),
      inherit.aes = FALSE, height = 0.5
    ) +
    ggnewscale::new_scale_fill()
  
  if (!is.null(gene_anno) && nrow(gene_anno) > 0L) {
    p <- p + ggplot2::geom_tile(
      data = gene_anno,
      ggplot2::aes(x = 0, y = Gene, fill = GeneCluster),
      inherit.aes = FALSE, width = 8
    )
  }
  
  p +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(title = x$title, x = NULL, y = NULL)
}


#' Heatmap for a \code{clustering_input} (pre-clustering)
#'
#' Quick visualization of the prepared expression matrix; does not display
#' cluster strips since clustering has not yet been run.
#'
#' @param x A \code{clustering_input}.
#' @param palette Color palette for the expression scale.
#' @param scale_rows Logical; if \code{TRUE}, z-score by gene before plotting.
#'
#' @return A \code{ggplot} heatmap.
#' @export
plot_heatmap.clustering_input <- function(
    x,
    palette = c("darkblue","blue","white","red","darkred"),
    scale_rows = TRUE
) {
  X <- x$Xtemp
  if (isTRUE(scale_rows)) {
    cn <- colnames(X); X <- t(apply(X, 1, scale)); colnames(X) <- cn
  }
  df <- tibble::as_tibble(X, rownames = "Gene") |>
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  rng <- max(abs(df$Expression), na.rm = TRUE)
  ggplot2::ggplot(df, ggplot2::aes(Sample, Gene, fill = Expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      values = scales::rescale(c(-rng, -3, 0, 3, rng)),
      limits = c(-rng, rng)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(title = x$title, x = NULL, y = NULL)
}




#' Kaplan–Meier survival curves with Cox model annotation
#'
#' This function fits and plots Kaplan–Meier survival curves for clusters of
#' samples, overlays Cox proportional hazards results, and includes a risk table.
#' Colors are assigned automatically and kept consistent between the survival
#' curves and the risk table.
#'
#' @param data A list that must include a `sampInfo` data frame with columns:
#'   - `time`: follow-up time
#'   - `event`: censoring indicator (1 = event, 0 = censored)
#'   - a cluster column specified by `cluster_name`
#'   - logical columns `whitelist` and `pdac` for filtering.
#' @param cluster_name A string giving the name of the cluster column inside
#'   `data$sampInfo`.
#' @param title A string for the plot title.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters samples to those with `whitelist = TRUE` and `pdac = TRUE`.
#'   \item Drops rows with missing values in `time`, `event`, or `cluster_name`.
#'   \item Fits Kaplan–Meier curves and a Cox model using \pkg{survival}.
#'   \item Chooses a palette using \pkg{ggplot2}'s default discrete hues,
#'         ensuring that colors match between curves and the risk table.
#'   \item Annotates the survival plot with the BIC of the Cox model and hazard
#'         ratios for each cluster relative to the reference (first factor level).
#'   \item Returns a combined figure (survival plot + risk table) using
#'         \pkg{ggpubr}.
#' }
#'
#' @return A `ggpubr::ggarrange` object containing:
#'   \itemize{
#'     \item the survival plot (Kaplan–Meier curves with annotation), and
#'     \item the corresponding risk table with matching colors.
#'   }
#'
#' @examples
#' \dontrun{
#' # Suppose `mydata` is a list with `sampInfo` including time, event, clusters:
#' p <- plot_survival(mydata, cluster_name = "tumor_clusters", title = "Survival by Cluster")
#' print(p)
#' # To save:
#' ggplot2::ggsave("survival_plot.png", p, width = 8, height = 6, dpi = 300)
#' }
#'
#' @import dplyr tidyr survival survminer ggpubr scales
#' @export
plot_survival <- function(data, cluster_name, title) {
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(ggpubr)
  library(scales)  # hue_pal()
  
  if (!"sampInfo" %in% names(data)) stop("data must contain 'sampInfo'")
  
  # Filter & basic checks
  sampInfo <- data$sampInfo %>%
    dplyr::filter(.data$whitelist, .data$pdac)
  
  if (!all(c("time", "event", cluster_name) %in% names(sampInfo))) {
    stop(sprintf(
      "Missing required columns. Needed: time, event, %s", cluster_name
    ))
  }
  
  # Drop rows with missing essentials
  before_n <- nrow(sampInfo)
  sampInfo <- sampInfo %>%
    dplyr::select(time, event, !!rlang::sym(cluster_name), dplyr::everything()) %>%
    dplyr::rename(cluster = !!rlang::sym(cluster_name)) %>%
    dplyr::filter(!is.na(time), !is.na(event), !is.na(cluster))
  removed_n <- before_n - nrow(sampInfo)
  if (removed_n > 0) {
    warning(sprintf("Removed %d rows with NA in time/event/cluster.", removed_n))
  }
  
  if (nrow(sampInfo) == 0) stop("No rows remain after filtering/NA removal.")
  if (dplyr::n_distinct(sampInfo$cluster) < 2) {
    stop("Need at least 2 clusters to draw KM curves.")
  }
  
  # Stable factor (explicit levels)
  sampInfo <- sampInfo %>% mutate(cluster = factor(cluster))
  levels_cluster <- levels(sampInfo$cluster)
  nclass <- length(levels_cluster)
  
  # Fit models (ties='efron' by default)
  fit  <- survfit(Surv(time, event) ~ cluster, data = sampInfo)
  fit2 <- coxph(Surv(time, event) ~ cluster, data = sampInfo)
  bic  <- BIC(fit2)
  
  # Build KM plot + risk table with a deterministic palette
  # First call without palette to learn exact 'strata' labels used by survminer
  splot0 <- ggsurvplot(
    fit,
    data = sampInfo,
    risk.table = TRUE,
    pval = TRUE
  )
  
  # strata labels look like "cluster=Level"
  strata_levels <- unique(splot0$data.survtable$strata)
  # Deterministic default palette from ggplot2 hues
  pal_vec <- hue_pal(h.start = 15, direction = 1)(length(strata_levels))
  pal_named <- setNames(pal_vec, strata_levels)
  
  # Rebuild ggsurvplot with our palette in the correct order
  splot <- ggsurvplot(
    fit,
    data       = sampInfo,
    risk.table = TRUE,
    pval       = TRUE,
    palette    = unname(pal_named[strata_levels])
  )
  
  # HR labels with correct reference
  csum <- summary(fit2)
  # Coefficient rownames are like "clusterLevel"
  coef_df <- data.frame(
    term   = rownames(csum$coefficients),
    HR     = csum$coefficients[, "exp(coef)"],
    check.names = FALSE
  )
  # Extract level names from terms
  coef_df$level <- sub("^cluster", "", coef_df$term)
  ref_level <- levels_cluster[1]
  hr_lines <- paste0(coef_df$level, " vs ", ref_level, ": HR = ", round(coef_df$HR, 2))
  ann_text <- paste0("BIC = ", round(bic, 2),
                     if (nrow(coef_df) > 0) paste0("\n", paste(hr_lines, collapse = "\n")) else "")
  
  # Risk table (transpose; color rows by strata)
  rt_wide <- splot$data.survtable %>%
    select(strata, time, n.risk) %>%
    tidyr::pivot_wider(names_from = strata, values_from = n.risk, id_cols = time)
  
  rt_t <- t(rt_wide)
  row_names <- rownames(rt_t)
  
  # First row is 'time'; keep it black. Others match by exact strata name.
  row_colors <- rep("black", length(row_names))
  if (length(row_names) > 1) {
    match_idx <- match(row_names[-1], names(pal_named))
    row_colors[-1] <- pal_named[match_idx]
    # Any unmatched (shouldn't happen) fall back to black
    row_colors[is.na(row_colors)] <- "black"
  }
  
  risk.table <- ggtexttable(
    rt_t,
    theme = ttheme(
      tbody.style = tbody_style(
        color = row_colors,
        fill  = "white",
        face  = "bold"
      )
    )
  )
  
  # Compose final plot with annotation placed relative to observed time range
  tmax <- max(splot$data.survtable$time, na.rm = TRUE)
  surv.plot <- splot$plot +
    labs(title = title) +
    annotate(
      "text",
      x = tmax - 0.05 * tmax,
      y = 0.85,
      label = ann_text,
      size = 4.5,
      hjust = 1
    )
  
  final_plot <- ggarrange(surv.plot, risk.table, nrow = 2, heights = c(3, 1))
  return(final_plot)
}
