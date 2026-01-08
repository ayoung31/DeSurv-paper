prepare_bo_history <- function(history_path) {
  bo_history <- read.csv(history_path, stringsAsFactors = FALSE)
  bo_history <- bo_history[order(bo_history$eval_id), , drop = FALSE]
  bo_history <- bo_history[bo_history$status == "ok" & !is.na(bo_history$mean_cindex), , drop = FALSE]
  bo_history$evaluation <- seq_len(nrow(bo_history))
  bo_history$stage <- factor(
    bo_history$stage,
    levels = c("init", "bo"),
    labels = c("Initial design", "BO iteration")
  )
  bo_history
}

make_panel <- function(history_df, labels, include_alpha = TRUE, cindex_label = "Best CV C-index across k") {
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  if (is.na(alpha_col)) {
    history_df$alpha_fixed <- 0
    alpha_col <- "alpha_fixed"
  }

  best_per_k <- history_df %>%
    dplyr::group_by(.data[[k_col]]) %>%
    dplyr::arrange(desc(mean_cindex), .data[[alpha_col]]) %>%
    dplyr::mutate(
      c_se = stats::sd(mean_cindex, na.rm = TRUE) /
        sqrt(sum(!is.na(mean_cindex)))
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      k = .data[[k_col]],
      c_best = mean_cindex,
      alpha_best = .data[[alpha_col]],
      c_se = c_se
    ) %>%
    dplyr::mutate(
      k = factor(k, levels = sort(unique(k)))
    )

  x_lab <- if (include_alpha) NULL else "Rank k"
  p_k_cindex <- ggplot2::ggplot(best_per_k, ggplot2::aes(x = k, y = c_best, group = 1)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = c_best - c_se, ymax = c_best + c_se),
      width = 0.15,
      color = "#1f78b4"
    ) +
    ggplot2::geom_line(color = "#1f78b4", linewidth = 0.6) +
    ggplot2::geom_point(color = "#1f78b4", size = 2.2) +
    ggplot2::labs(x = x_lab, y = cindex_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  if (include_alpha) {
    p_k_alpha <- ggplot2::ggplot(best_per_k, ggplot2::aes(x = k, y = alpha_best, group = 1)) +
      ggplot2::geom_line(color = "#33a02c", linewidth = 0.6) +
      ggplot2::geom_point(color = "#33a02c", size = 2.2) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Rank k", y = "Selected alpha") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "none")

    cowplot::plot_grid(
      p_k_cindex,
      p_k_alpha,
      ncol = 1,
      align = "v",
      rel_heights = c(1, 1),
      labels = labels
    )
  } else {
    p_k_cindex + ggplot2::labs(tag = labels[1])
  }
}

make_cindex_violin <- function(history_df, label) {
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  k_levels <- sort(unique(history_df[[k_col]]))
  history_df <- history_df %>%
    dplyr::mutate(k = factor(.data[[k_col]], levels = k_levels))

  ggplot2::ggplot(history_df, ggplot2::aes(x = k, y = mean_cindex)) +
    ggplot2::geom_violin(fill = "#a6cee3", color = NA, scale = "width", alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.4, alpha = 0.7, fill = "#1f78b4", color = "#1f78b4") +
    ggplot2::labs(
      x = "Rank k",
      y = "CV C-index distribution",
      tag = label
    ) +
    ggplot2::theme_minimal(base_size = 9)
}

normalize_gp_params <- function(param_df, bounds_df) {
  for (param in names(param_df)) {
    bound_row <- bounds_df[bounds_df$parameter == param, , drop = FALSE]
    if (!nrow(bound_row)) {
      next
    }
    lower <- as.numeric(bound_row$lower[[1]])
    upper <- as.numeric(bound_row$upper[[1]])
    scale_type <- bound_row$scale[[1]]
    if (!is.na(scale_type) && identical(scale_type, "log10")) {
      param_df[[param]] <- log10(param_df[[param]])
      lower <- log10(lower)
      upper <- log10(upper)
    }
    param_df[[param]] <- (param_df[[param]] - lower) / (upper - lower)
  }
  param_df
}

extract_gp_curve <- function(bo_results, history_df) {
  runs <- bo_results[["runs"]]
  last_run <- runs[[length(runs)]]
  km_fit <- last_run[["km_fit"]]
  bounds <- last_run[["bounds"]]
  param_names <- colnames(km_fit@X)
  if (is.null(param_names)) {
    stop("GP design matrix has no column names.")
  }

  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  best_per_k <- history_df %>%
    dplyr::group_by(.data[[k_col]]) %>%
    dplyr::arrange(desc(mean_cindex)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(k_numeric = as.numeric(.data[[k_col]])) %>%
    dplyr::arrange(k_numeric)

  newdata_actual <- best_per_k[, param_names, drop = FALSE]
  newdata_scaled <- normalize_gp_params(newdata_actual, bounds)

  preds <- DiceKriging::predict(
    km_fit,
    newdata = newdata_scaled,
    type = "UK",
    se.compute = TRUE,
    cov.compute = FALSE
  )

  tibble::tibble(
    k = best_per_k$k_numeric,
    mean = preds$mean,
    lower = preds$mean - 1.96 * preds$sd,
    upper = preds$mean + 1.96 * preds$sd
  )
}

make_gp_curve_plot <- function(curve_df, label) {
  ggplot2::ggplot(curve_df, ggplot2::aes(x = k, y = mean, group = 1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#b2df8a",
      alpha = 0.4
    ) +
    ggplot2::geom_line(color = "#33a02c", linewidth = 0.7) +
    ggplot2::geom_point(color = "#33a02c", size = 2) +
    ggplot2::scale_x_continuous(breaks = curve_df$k) +
    ggplot2::labs(
      x = "Rank k",
      y = "GP-predicted CV C-index",
      tag = label
    ) +
    ggplot2::theme_minimal(base_size = 9)
}

save_fig_bo <- function(bo_history_path, bo_history_alpha0_path, bo_results_supervised,
                        bo_results_alpha0, fit_std, path, width = 6, height = 5.5) {
  bo_history_supervised <- prepare_bo_history(bo_history_path)
  bo_history_alpha0 <- prepare_bo_history(bo_history_alpha0_path)

  p_panelA <- make_panel(bo_history_supervised, labels = c("A", "B"))
  p_cindex_violin <- make_cindex_violin(bo_history_supervised, label = "C")
  p_cindex_violin_alpha0 <- make_cindex_violin(bo_history_alpha0, label = "D")

  gp_curve_supervised <- extract_gp_curve(bo_results_supervised, bo_history_supervised)
  gp_curve_alpha0 <- extract_gp_curve(bo_results_alpha0, bo_history_alpha0)

  p_gp_supervised <- make_gp_curve_plot(gp_curve_supervised, label = "E")
  p_gp_alpha0 <- make_gp_curve_plot(gp_curve_alpha0, label = "F")

  plots_std <- plot(
    fit_std,
    what = c("cophenetic", "dispersion", "evar", "residuals", "silhouette", "sparseness"),
    main = NULL,
    xlab = "Rank (k)"
  ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(2, 12, by = 2)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  panel_right <- cowplot::plot_grid(
    p_cindex_violin,
    p_cindex_violin_alpha0,
    ncol = 1
  )

  panel_gp <- cowplot::plot_grid(
    p_gp_supervised,
    p_gp_alpha0,
    ncol = 2
  )

  panel_top <- cowplot::plot_grid(
    p_panelA,
    panel_right,
    ncol = 2,
    rel_widths = c(1.2, 0.8)
  )

  fig <- cowplot::plot_grid(
    panel_top,
    panel_gp,
    plots_std,
    ncol = 1,
    rel_heights = c(1.8, 0.8, 1),
    labels = c("", "", "G")
  )

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, fig, width = width, height = height, units = "in")
  path
}

save_fig_bio <- function(ora_analysis, fit_desurv, tops_desurv, top_genes_ref,
                         path, width = 7, height = 4.5) {
  ora_results <- ora_analysis$enrich_GO
  if (is.null(ora_results) || !length(ora_results)) {
    stop("ORA results are missing.")
  }

  p1 <- vector("list", length(ora_results))
  for (i in seq_along(ora_results)) {
    ORA_GO <- ora_results[[i]]
    if (nrow(ORA_GO@result) > 0) {
      p1[[i]] <- enrichplot::dotplot(ORA_GO, showCategory = 10, label_format = 100) +
        ggplot2::theme_minimal(base_size = 7) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.margin = ggplot2::margin(0, 5, 0, 10)
        ) +
        ggplot2::scale_y_discrete(
          labels = function(x) stringr::str_trunc(x, 30)
        ) +
        ggplot2::scale_size(range = c(0.5, 2))
    } else {
      p1[[i]] <- NULL
    }
  }

  p1 <- p1[!vapply(p1, is.null, logical(1))]
  if (length(p1) < 3) {
    stop("Need at least 3 ORA dotplots to build the ORA panel.")
  }

  all_dat <- do.call(rbind, lapply(p1[1:3], function(p) p$data))

  x_max <- max(all_dat$GeneRatio, na.rm = TRUE)
  col_min <- min(all_dat$p.adjust, na.rm = TRUE)
  col_max <- max(all_dat$p.adjust, na.rm = TRUE)
  size_min <- min(all_dat$Count, na.rm = TRUE)
  size_max <- max(all_dat$Count, na.rm = TRUE)

  make_equal_scales <- function(p) {
    p +
      ggplot2::scale_x_continuous(limits = c(0, x_max)) +
      ggplot2::scale_color_viridis_c(
        limits = c(col_min, col_max),
        direction = -1,
        name = "Adjusted p-value"
      ) +
      ggplot2::scale_size_continuous(
        limits = c(size_min, size_max),
        name = "Gene count"
      ) +
      ggplot2::theme(legend.position = "none")
  }

  p1_eq <- make_equal_scales(p1[[1]])
  p2_eq <- make_equal_scales(p1[[2]])
  p3_eq <- make_equal_scales(p1[[3]])

  p_for_legend <- p1[[1]] +
    ggplot2::scale_x_continuous(limits = c(0, x_max)) +
    ggplot2::scale_color_viridis_c(
      limits = c(col_min, col_max),
      direction = -1,
      name = "Adjusted p-value"
    ) +
    ggplot2::scale_size_continuous(
      limits = c(size_min, size_max),
      name = "Gene count"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 6)
    )

  g <- ggplot2::ggplotGrob(p_for_legend)
  legend_shared <- g$grobs[which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")][[1]]

  dplots2 <- cowplot::plot_grid(
    p1_eq,
    p2_eq,
    p3_eq,
    legend_shared,
    nrow = 2,
    align = "v",
    axis = "tblr",
    labels = c("A.", "B.", "C.")
  )

  if (is.null(top_genes_ref) || !length(top_genes_ref)) {
    stop("Reference gene signatures are missing.")
  }

  top_genes_local <- top_genes_ref
  top_genes_local$deCAF <- list(
    proCAF = c("IGFL2", "NOX4", "VSNL1", "BICD1", "NPR3", "ETV1", "ITGA11", "CNIH3", "COL11A1"),
    restCAF = c("CHRDL1", "OGN", "PI16", "ANK2", "ABCA8", "TGFBR3", "FBLN5", "SCARA5", "KIAA1217")
  )
  if (length(top_genes_local) >= 3) names(top_genes_local)[3] <- "Moffitt"
  if (length(top_genes_local) >= 4) names(top_genes_local)[4] <- "Moffitt"
  if (length(top_genes_local) >= 13) names(top_genes_local)[13] <- "SCISSORS"
  if (length(top_genes_local) >= 16) names(top_genes_local)[16] <- "SCISSORS"
  if (length(top_genes_local) >= 12) names(top_genes_local)[12] <- "Elyada"

  temp <- purrr::list_flatten(top_genes_local)
  ref_sigs <- temp[
    !startsWith(names(temp), "Bailey") &
      !grepl("peri", names(temp)) &
      !startsWith(names(temp), "DECODER") &
      !startsWith(names(temp), "MSI") &
      !startsWith(names(temp), "PurISS")
  ]

  W <- fit_desurv$W
  tops <- get_top_genes(W, 100)$top_genes
  W <- W[unlist(tops), , drop = FALSE]
  common_genes <- Reduce(intersect, list(rownames(W), unique(unlist(ref_sigs))))
  W <- W[common_genes, , drop = FALSE]

  cor_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat_adj <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )

  for (j in seq_len(ncol(W))) {
    wj <- W[, j]
    for (k in seq_along(ref_sigs)) {
      vk <- as.numeric(common_genes %in% ref_sigs[[k]])
      cor_mat[k, j] <- stats::cor(wj, vk, method = "spearman")
      p_mat[k, j] <- stats::cor.test(wj, vk, method = "spearman")$p.value
    }
    p_mat_adj[, j] <- stats::p.adjust(p_mat[, j], method = "BH")
  }

  keep <- vapply(seq_len(nrow(cor_mat)), function(j) {
    !any(is.na(cor_mat[j, ])) & sum(p_mat_adj[j, ] < 0.05) > 0
  }, logical(1))
  mat <- cor_mat[which(keep), , drop = FALSE]

  my_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  ph <- pheatmap::pheatmap(
    mat,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = seq(-0.4, 0.4, length.out = 101),
    fontsize = 6,
    silent = TRUE,
    fontsize_number = 18,
    treeheight_row = 0
  )
  ph_grob <- ph$gtable
  pheat <- cowplot::plot_grid(NULL, cowplot::ggdraw(ph_grob), nrow = 2, rel_heights = c(0.25, 4))

  fig <- cowplot::plot_grid(dplots2, pheat, ncol = 2, rel_widths = c(2, 1), labels = c("", "D."))

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, fig, width = width, height = height, units = "in")
  path
}

get_vam_scores <- function(sc, desurv_genesets) {
  DefaultAssay(sc) <- "RNA"

  gene_ids <- rownames(sc)
  gs_collection <- createGeneSetCollection(
    gene.ids = gene_ids,
    gene.set.collection = desurv_genesets,
    min.size = 5
  )

  sc <- vamForSeurat(
    seurat.data = sc,
    gene.set.collection = gs_collection,
    center = FALSE,
    gamma = TRUE,
    sample.cov = FALSE,
    return.dist = FALSE
  )
  DefaultAssay(sc) <- "VAMcdf"
  sc
}

save_fig_sc <- function(tops_desurv, sc_all_path, sc_caf_path, sc_tum_path,
                        path, width = 7.5, height = 8) {
  if (!file.exists(sc_all_path)) {
    stop("Missing scRNA-seq file: ", sc_all_path)
  }
  if (!file.exists(sc_caf_path)) {
    stop("Missing scRNA-seq file: ", sc_caf_path)
  }
  if (!file.exists(sc_tum_path)) {
    stop("Missing scRNA-seq file: ", sc_tum_path)
  }

  desurv_genesets <- as.list(tops_desurv$top_genes)
  if (length(desurv_genesets) < 3) {
    stop("Need at least 3 DeSurv factors to build the scRNA-seq figure.")
  }

  factor_names <- paste0("DeSurv Factor ", seq_along(desurv_genesets))
  names(desurv_genesets) <- factor_names
  features_to_plot <- factor_names[seq_len(3)]

  sc_all <- readRDS(sc_all_path)
  sc_caf <- readRDS(sc_caf_path)
  sc_tum <- readRDS(sc_tum_path)

  sc_all <- get_vam_scores(sc_all, desurv_genesets)
  sc_caf <- get_vam_scores(sc_caf, desurv_genesets)
  sc_tum <- get_vam_scores(sc_tum, desurv_genesets)

  p_list_all <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_all,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data",
      max.cutoff = "q95"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-15, 0, -5, 0)
      ) +
      ggplot2::labs(color = "VAM\nscore")
    fp[[1]]
  })

  legend <- cowplot::get_legend(
    p_list_all[[1]] +
      ggplot2::theme(
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 6),
        legend.title = ggplot2::element_text(size = 6)
      )
  )
  p_list_all[[1]] <- p_list_all[[1]] + ggplot2::theme(legend.position = "none")
  p_list_all[[2]] <- p_list_all[[2]] + ggplot2::theme(legend.position = "none")
  p_list_all[[3]] <- p_list_all[[3]] + ggplot2::theme(legend.position = "none")

  p_list_caf <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_caf,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data",
      max.cutoff = "q95"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-20, 0, -5, 0),
        legend.position = "none"
      )
    fp[[1]]
  })

  p_list_tum <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_tum,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-20, 0, -5, 0),
        legend.position = "none"
      )
    fp[[1]]
  })

  scplot_all <- DimPlot(
    sc_all,
    group.by = "label_broad",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("All Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  scplot_caf <- DimPlot(
    sc_caf,
    group.by = "label",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("CAF Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  scplot_tum <- DimPlot(
    sc_tum,
    group.by = "label",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("PDAC Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  labelmid <- grid::textGrob("DeSurv\nfactor 1", gp = grid::gpar(fontsize = 8, fontface = "bold"))
  labelmid2 <- grid::textGrob("DeSurv\nfactor 2", gp = grid::gpar(fontsize = 8, fontface = "bold"))
  labelbottom <- grid::textGrob("DeSurv\nfactor 3", gp = grid::gpar(fontsize = 8, fontface = "bold"))

  top <- cowplot::plot_grid(
    scplot_all[[1]],
    scplot_caf[[1]],
    scplot_tum[[1]],
    ncol = 3,
    rel_widths = c(1, 1, 1)
  )
  mid <- cowplot::plot_grid(
    labelmid,
    p_list_all[[1]],
    p_list_caf[[1]],
    p_list_tum[[1]],
    ncol = 4,
    rel_widths = c(0.5, 1, 1, 1)
  )
  mid2 <- cowplot::plot_grid(
    labelmid2,
    p_list_all[[2]],
    p_list_caf[[2]],
    p_list_tum[[2]],
    ncol = 4,
    rel_widths = c(0.5, 1, 1, 1)
  )
  bottom <- cowplot::plot_grid(
    labelbottom,
    p_list_all[[3]],
    p_list_caf[[3]],
    p_list_tum[[3]],
    ncol = 4,
    rel_widths = c(0.5, 1, 1, 1)
  )

  add_outline <- function(p) {
    p + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6))
  }

  full <- cowplot::plot_grid(
    add_outline(mid),
    NULL,
    add_outline(mid2),
    NULL,
    add_outline(bottom),
    nrow = 5,
    labels = c("B.", "", "C.", "", "D."),
    rel_heights = c(1, 0.05, 1, 0.05, 1)
  )

  vam_mat <- t(as.matrix(GetAssayData(sc_all, assay = "VAMcdf", layer = "data")))
  sc_all@meta.data[, colnames(vam_mat)] <- vam_mat[rownames(sc_all@meta.data), ]

  ct_col <- "label_fine"
  if (!ct_col %in% colnames(sc_all@meta.data)) {
    stop("Expected a column named 'label_fine' in sc_all@meta.data.")
  }

  avg_scores <- sc_all@meta.data %>%
    dplyr::group_by(.data[[ct_col]]) %>%
    dplyr::summarise(
      dplyr::across(features_to_plot, ~ mean(.x, na.rm = TRUE))
    )

  mat <- as.matrix(avg_scores[, -1, drop = FALSE])
  rownames(mat) <- avg_scores[[ct_col]]

  mat_capped <- mat
  upper <- stats::quantile(mat, 0.99, na.rm = TRUE)
  lower <- stats::quantile(mat, 0.01, na.rm = TRUE)
  mat_capped[mat_capped > upper] <- upper
  mat_capped[mat_capped < lower] <- lower

  col_fun <- viridis::viridis(256, option = "D")

  ht <- pheatmap::pheatmap(
    mat_capped,
    color = col_fun,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    legend = FALSE,
    silent = TRUE,
    fontsize = 6
  )

  gght <- cowplot::plot_grid(NULL, ggplotify::as.ggplot(ht$gtable), nrow = 2, rel_heights = c(0.25, 4))
  leg <- cowplot::plot_grid(NULL, legend, nrow = 2, rel_heights = c(3, 2))
  main <- cowplot::plot_grid(
    full,
    NULL,
    gght,
    leg,
    ncol = 4,
    rel_widths = c(8, 0.25, 3, 0.8),
    labels = c("", "", "E.", ""),
    label_x = c(0, 0, -0.09, 0)
  )

  fig <- cowplot::plot_grid(top, main, nrow = 2, rel_heights = c(1.5, 3))

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, fig, width = width, height = height, units = "in")
  path
}

select_desurv_factors <- function(fit, n = 2) {
  beta <- fit$beta
  if (is.null(beta)) {
    stop("DeSurv fit has no beta coefficients.")
  }
  beta <- as.numeric(beta)
  if (!length(beta)) {
    stop("DeSurv fit has empty beta coefficients.")
  }
  ord <- order(abs(beta), decreasing = TRUE)
  ord[seq_len(min(n, length(ord)))]
}

subset_validation_data <- function(dataset) {
  if (!is.null(dataset$samp_keeps) && length(dataset$samp_keeps)) {
    keeps <- dataset$samp_keeps
    if (is.logical(keeps)) {
      keeps <- which(keeps)
    } else if (is.character(keeps)) {
      keeps <- match(keeps, colnames(dataset$ex))
      keeps <- keeps[!is.na(keeps)]
    }
    if (length(keeps)) {
      dataset$ex <- dataset$ex[, keeps, drop = FALSE]
      dataset$sampInfo <- dataset$sampInfo[keeps, , drop = FALSE]
    }
  }
  dataset
}

save_fig_clus <- function(data_val_filtered, fit_desurv, tops_desurv, path,
                          width = 7.5, height = 4.5) {
  if (is.null(data_val_filtered) || !length(data_val_filtered)) {
    stop("No validation data supplied for clustering.")
  }

  factors <- select_desurv_factors(fit_desurv, n = 2)
  if (length(factors) < 2) {
    stop("Need at least two factors for clustering.")
  }

  clus_entries <- mapply(
    function(dataset, name) {
      dataset <- subset_validation_data(dataset)
      if (is.null(dataset$sampInfo$dataset)) {
        dataset$sampInfo$dataset <- name
      }
      run_clustering(
        tops = tops_desurv$top_genes,
        data = dataset,
        gene_lists = list(),
        type = "bas/clas",
        facs = factors,
        maxKcol = 3,
        maxKrow = length(factors),
        plot = FALSE,
        save = FALSE
      )
    },
    dataset = data_val_filtered,
    name = names(data_val_filtered),
    SIMPLIFY = FALSE
  )

  cluster_col <- paste0("facs_", paste(factors, collapse = "_"), "_with_3clusters")
  for (i in seq_along(clus_entries)) {
    sampInfo <- clus_entries[[i]]$data$sampInfo
    sampInfo$samp_cluster <- sampInfo[[cluster_col]]
    clus_entries[[i]]$data$sampInfo <- sampInfo
  }

  sinfos <- lapply(clus_entries, function(entry) entry$data$sampInfo)
  sampInfo <- dplyr::bind_rows(sinfos)
  sampInfo <- sampInfo[!is.na(sampInfo$samp_cluster), , drop = FALSE]

  fit <- survival::survfit(survival::Surv(time, event) ~ samp_cluster, data = sampInfo)
  p <- survminer::ggsurvplot(
    fit,
    data = sampInfo,
    pval = TRUE,
    risk.table = TRUE,
    legend.labs = c("1", "2", "3"),
    legend.title = "Cluster",
    pval.coord = c(0, 0.05),
    pval.size = 2.5,
    censor.size = 2,
    risk.table.fontsize = 2
  )

  ex_list <- list()
  genes_list <- list()
  samps_list <- list()
  tops <- unlist(tops_desurv$top_genes[, factors, drop = FALSE])
  for (i in seq_along(clus_entries)) {
    data_entry <- clus_entries[[i]]$data
    ex_mat <- as.matrix(data_entry$ex)
    storage.mode(ex_mat) <- "numeric"
    ex_list[[i]] <- ex_mat
    genes_list[[i]] <- intersect(rownames(ex_list[[i]]), tops)
    samps_list[[i]] <- data_entry$sampInfo
  }

  genes <- Reduce(intersect, genes_list)
  if (!length(genes)) {
    stop("No shared genes available across validation datasets.")
  }

  ex_list <- lapply(ex_list, function(x) t(scale(t(x[genes, , drop = FALSE]))))
  X <- do.call("cbind", ex_list)
  X <- X[genes, , drop = FALSE]
  sampInfo <- do.call("rbind", samps_list)

  n_samples <- nrow(sampInfo)
  get_meta_col <- function(name) {
    if (!name %in% names(sampInfo)) {
      return(rep(NA_character_, n_samples))
    }
    vals <- sampInfo[[name]]
    if (length(vals) != n_samples) {
      return(rep(NA_character_, n_samples))
    }
    vals
  }
  col_anno <- data.frame(
    DeSurv_cluster = as.factor(sampInfo$samp_cluster),
    DeCAF = get_meta_col("DeCAF"),
    PurIST = get_meta_col("PurIST"),
    dataset = get_meta_col("dataset")
  )
  col_anno$DeCAF <- ifelse(col_anno$DeCAF == "permCAF", "proCAF", col_anno$DeCAF)
  rownames(col_anno) <- colnames(X)
  X <- X[, order(col_anno$DeSurv_cluster)]
  col_anno <- col_anno[order(col_anno$DeSurv_cluster), , drop = FALSE]

  row_anno <- data.frame(gene = genes)
  row_anno$`DeSurv factor` <- as.factor(
    ifelse(genes %in% tops_desurv$top_genes[, factors[1]], factors[1], factors[2])
  )
  rownames(row_anno) <- row_anno$gene
  row_anno$gene <- NULL
  X <- X[order(row_anno$`DeSurv factor`), , drop = FALSE]
  row_anno <- row_anno[order(row_anno$`DeSurv factor`), , drop = FALSE]

  annotation_colors <- list(
    `DeSurv factor` = stats::setNames(
      c("lightgrey", "black"),
      as.character(factors[1:2])
    ),
    DeSurv_cluster = c(
      `1` = "#F8766D",
      `2` = "#00BA38",
      `3` = "#619CFF"
    ),
    DeCAF = c(
      proCAF = "violetred2",
      restCAF = "cyan4"
    ),
    PurIST = c(
      `Basal-like` = "orange",
      Classical = "blue"
    )
  )
  if (!all(is.na(col_anno$dataset))) {
    annotation_colors$dataset <- c(
      Dijk = "slateblue1",
      Moffitt_GEO_array = "springgreen4",
      PACA_AU_array = "yellow3",
      PACA_AU_seq = "coral",
      Puleo_array = "dodgerblue3"
    )
  }

  min_val <- -2
  max_val <- 3.5
  ncolors <- 500
  my_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(ncolors)

  breaks_centered <- c(
    seq(min_val, 0, length.out = ceiling(ncolors / 2) + 1),
    seq(0, max_val, length.out = floor(ncolors / 2) + 1)[-1]
  )

  ph <- pheatmap::pheatmap(
    X,
    annotation_col = col_anno,
    annotation_row = row_anno,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = breaks_centered,
    show_colnames = FALSE,
    annotation_names_row = FALSE,
    show_rownames = FALSE,
    silent = TRUE,
    fontsize = 6
  )

  ph_grob <- ph$gtable
  pheat <- cowplot::ggdraw(ph_grob) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, -10, 0, 0))

  sampInfo$subtype <- ifelse(
    sampInfo$PurIST == "Basal-like",
    "Basal-like",
    ifelse(sampInfo$DeCAF == "restCAF", "Classical + restCAF", "Classical + proCAF")
  )
  tbl <- table(sampInfo$samp_cluster, sampInfo$subtype)

  gt_tbl <- gt::gt(as.data.frame.matrix(tbl), rownames_to_stub = TRUE) |>
    gt::cols_width(gt::everything() ~ gt::px(120)) |>
    gt::tab_options(table.font.size = gt::px(18)) |>
    gt::cols_align(align = "center", columns = gt::everything()) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#F8766D"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = 1, columns = `Basal-like`)
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#00BA38"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = 2, columns = `Classical + restCAF`)
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#619CFF"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = 3, columns = `Classical + proCAF`)
    )

  table_path <- tempfile(fileext = ".png")
  gt::gtsave(gt_tbl, table_path)
  img <- magick::image_read(table_path)
  table_plot <- cowplot::ggdraw() + cowplot::draw_image(img)

  sp <- p$plot +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = ggplot2::margin(20, 5, 0, 0)
    )
  st <- p$table +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 40, 0, 10))

  right <- cowplot::plot_grid(
    sp,
    st,
    NULL,
    nrow = 3,
    rel_heights = c(4, 2, 1),
    labels = c("C.", "D."),
    axis = "tblr"
  )
  tab <- cowplot::plot_grid(table_plot, NULL, ncol = 2, rel_widths = c(3, 2))
  left <- cowplot::plot_grid(
    NULL,
    pheat,
    NULL,
    NULL,
    NULL,
    tab,
    rel_heights = c(4, 0.05, 1.5),
    ncol = 2,
    labels = c("A.", "", "B.", ""),
    rel_widths = c(0.1, 1)
  )

  fig <- cowplot::plot_grid(left, right, rel_widths = c(3, 2))

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, fig, width = width, height = height, units = "in")
  path
}
