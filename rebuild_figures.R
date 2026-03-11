#!/usr/bin/env Rscript
# rebuild_figures.R — Fix ggplot2 version mismatch (3.x → 4.x)
#
# The targets store contains ggplot objects serialized with ggplot2 3.x on HPC.
# Local ggplot2 4.x cannot ggplot_build() these objects. This script:
#   1. Reads each serialized figure object
#   2. Extracts the underlying data or grobs
#   3. Rebuilds fresh ggplot2 4.x-compatible objects
#   4. Saves back to the store
#
# Run once after syncing store objects from HPC:
#   Rscript rebuild_figures.R

library(ggplot2)
library(cowplot)
library(ggrepel)
library(viridis)
library(scales)

store_dir <- file.path(
  "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main", "objects"
)

read_obj <- function(name) {
  readRDS(file.path(store_dir, name))
}

save_obj <- function(obj, name) {
  saveRDS(obj, file.path(store_dir, name))
  cat("  Saved:", name, "\n")
}

# Helper: extract gtable grob from a cowplot::ggdraw() wrapper
extract_grob <- function(ggdraw_obj) {
  ggdraw_obj$layers[[1]]$geom_params$grob
}

# Helper: recursively set absolute font sizes in pheatmap grobs (idempotent)
# Sets all text grob fontsizes to `size` (or `title_size` for title elements)
set_grob_text_size <- function(g, size = 7, title_size = 9) {
  if (inherits(g, "text") && !is.null(g$gp$fontsize)) {
    is_title <- !is.null(g$rot) && g$rot == 0 &&
                length(g$label) == 1 && nchar(g$label) <= 15 &&
                !grepl(":", g$label, fixed = TRUE)
    g$gp$fontsize <- if (is_title) title_size else size
  }
  if (!is.null(g$children)) {
    for (i in seq_along(g$children)) g$children[[i]] <- set_grob_text_size(g$children[[i]], size, title_size)
  }
  if (!is.null(g$grobs)) {
    for (i in seq_along(g$grobs)) g$grobs[[i]] <- set_grob_text_size(g$grobs[[i]], size, title_size)
  }
  g
}

# Helper: set fonts + rotate bottom column labels + keep row labels horizontal
# Used for main-text heatmaps where row labels need maximum room
set_heatmap_text <- function(g, label_size = 7, title_size = 9, col_angle = 315) {
  if (inherits(g, "text") && !is.null(g$gp$fontsize)) {
    nlbl <- length(g$label)
    has_colon <- any(grepl(":", g$label, fixed = TRUE))
    if (nlbl == 1 && !has_colon && nchar(g$label) <= 15) {
      g$gp$fontsize <- title_size
    } else if (nlbl > 1 && has_colon) {
      # Row labels (gene programs): horizontal
      g$gp$fontsize <- label_size
      g$rot <- 0
      g$just <- c("left", "center")
    } else if (nlbl > 1 && !has_colon) {
      # Column labels (factor names): angled
      g$gp$fontsize <- label_size
      g$rot <- col_angle
      g$just <- c("right", "center")
    }
  }
  if (!is.null(g$children)) {
    for (i in seq_along(g$children)) g$children[[i]] <- set_heatmap_text(g$children[[i]], label_size, title_size, col_angle)
  }
  if (!is.null(g$grobs)) {
    for (i in seq_along(g$grobs)) g$grobs[[i]] <- set_heatmap_text(g$grobs[[i]], label_size, title_size, col_angle)
  }
  g
}

# Helper: narrow the matrix cells in a pheatmap gtable to give row labels more room
narrow_matrix <- function(g, matrix_frac = 0.30) {
  if (inherits(g, "gtable") && !is.null(g$layout) && "matrix" %in% g$layout$name) {
    g$widths[3] <- unit(matrix_frac, "npc")
    g$widths[4] <- unit(1 - matrix_frac, "npc") - unit(5, "bigpts")
    return(g)
  }
  if (!is.null(g$children)) {
    for (i in seq_along(g$children)) g$children[[i]] <- narrow_matrix(g$children[[i]], matrix_frac)
  }
  if (!is.null(g$grobs)) {
    for (i in seq_along(g$grobs)) g$grobs[[i]] <- narrow_matrix(g$grobs[[i]], matrix_frac)
  }
  g
}

# Helper: remove a named row from a pheatmap grob (matrix rects + row labels)
remove_heatmap_row <- function(g, label_pattern) {
  if (inherits(g, "gtable") && !is.null(g$layout) && "matrix" %in% g$layout$name) {
    rn_idx <- which(g$layout$name == "row_names")
    rn_grob <- g$grobs[[rn_idx]]

    target <- grep(label_pattern, rn_grob$label, fixed = TRUE)
    if (!length(target)) return(g)

    target_y <- as.numeric(rn_grob$y[target])
    n_old <- length(rn_grob$label)
    n_new <- n_old - 1

    # Update row labels
    keep_lbl <- setdiff(seq_len(n_old), target)
    rn_grob$label <- rn_grob$label[keep_lbl]
    if (length(rn_grob$x) > 1) rn_grob$x <- rn_grob$x[keep_lbl]
    rn_grob$y <- unit(rev(seq(1, 2*n_new - 1, by = 2) / (2 * n_new)), "npc")
    g$grobs[[rn_idx]] <- rn_grob

    # Update matrix rects
    mat_idx <- which(g$layout$name == "matrix")
    rects <- g$grobs[[mat_idx]]$children[[1]]
    rect_ys <- as.numeric(rects$y)
    keep_cells <- which(abs(rect_ys - target_y) >= 0.001)
    old_kept_ys <- rect_ys[keep_cells]
    old_unique_ys <- sort(unique(old_kept_ys), decreasing = TRUE)
    new_cell_ys <- rev(seq(1, 2*n_new - 1, by = 2) / (2 * n_new))
    y_map <- setNames(new_cell_ys, as.character(old_unique_ys))
    new_ys <- y_map[as.character(round(old_kept_ys, 8))]
    if (any(is.na(new_ys))) {
      for (j in which(is.na(new_ys)))
        new_ys[j] <- new_cell_ys[which.min(abs(old_unique_ys - old_kept_ys[j]))]
    }
    rects$x <- rects$x[keep_cells]
    rects$y <- unit(as.numeric(new_ys), "npc")
    rects$height <- unit(1/n_new, "npc")
    rects$gp$fill <- rects$gp$fill[keep_cells]
    g$grobs[[mat_idx]]$children[[1]] <- rects
    return(g)
  }
  if (!is.null(g$children)) for (i in seq_along(g$children)) g$children[[i]] <- remove_heatmap_row(g$children[[i]], label_pattern)
  if (!is.null(g$grobs)) for (i in seq_along(g$grobs)) g$grobs[[i]] <- remove_heatmap_row(g$grobs[[i]], label_pattern)
  g
}

cat("=== Rebuilding figure objects for ggplot2", as.character(packageVersion("ggplot2")), "===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Gene overlap heatmaps — extract grob from ggdraw, rewrap
# ─────────────────────────────────────────────────────────────────────────────
cat("1. Gene overlap heatmaps (grob extraction + rewrap)\n")

# Main text heatmaps: angled bottom labels + narrowed matrix for label room
for (name in c("fig_gene_overlap_heatmap_desurv_tcgacptac",
               "fig_gene_overlap_heatmap_std_desurvk_tcgacptac")) {
  if (!file.exists(file.path(store_dir, name))) {
    cat("  SKIP (not in store):", name, "\n")
    next
  }
  obj <- read_obj(name)
  grob <- extract_grob(obj$plot)
  grob <- remove_heatmap_row(grob, "Bailey: Not Unique")
  grob <- set_heatmap_text(grob, label_size = 7, title_size = 9, col_angle = 315)
  grob <- narrow_matrix(grob, matrix_frac = 0.40)
  obj$plot <- cowplot::ggdraw(grob)
  save_obj(obj, name)
}

# Supplement heatmaps: simple font fix only (more columns, less crowded)
for (name in c("fig_gene_overlap_heatmap_desurv_alpha0_tcgacptac",
               "fig_gene_overlap_heatmap_std_elbowk_tcgacptac")) {
  if (!file.exists(file.path(store_dir, name))) {
    cat("  SKIP (not in store):", name, "\n")
    next
  }
  obj <- read_obj(name)
  grob <- extract_grob(obj$plot)
  grob <- set_grob_text_size(grob, size = 7, title_size = 9)
  obj$plot <- cowplot::ggdraw(grob)
  save_obj(obj, name)
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. Desurv-vs-NMF W-matrix correlation heatmap — extract grob, rewrap
# ─────────────────────────────────────────────────────────────────────────────
cat("\n2. W-matrix correlation heatmap (grob extraction + rewrap)\n")

obj <- read_obj("fig_desurv_std_correlation_tcgacptac")
grob <- extract_grob(obj$plot)
obj$plot <- cowplot::ggdraw(grob)
# Legend is already a gTree grob
save_obj(obj, "fig_desurv_std_correlation_tcgacptac")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Variance vs survival scatterplot — rebuild from data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n3. Variance vs survival scatterplot (rebuild from data)\n")

obj <- read_obj("fig_variation_explained_tcgacptac")

# Backup path for the extracted data (survives across rebuilds)
backup_path <- file.path(store_dir, ".df_plot_variation_explained_backup.rds")

# Try to extract df_plot from original HPC ggplot (only works on first run)
df_plot <- tryCatch(obj$plot_env$df_plot, error = function(e) NULL)

if (!is.null(df_plot) && is.data.frame(df_plot) && nrow(df_plot) > 0) {
  # First run: save backup for future idempotent rebuilds
  saveRDS(df_plot, backup_path)
  cat("  Saved df_plot backup (", nrow(df_plot), " rows)\n")
} else if (file.exists(backup_path)) {
  # Subsequent runs: use saved backup
  df_plot <- readRDS(backup_path)
  cat("  Using saved df_plot backup (", nrow(df_plot), " rows)\n")
} else if (inherits(obj, "ggplot") && is.data.frame(obj$data) && nrow(obj$data) > 0) {
  # Already rebuilt with good data: skip
  cat("  Already rebuilt with embedded data. Skipping.\n")
  df_plot <- NULL
} else {
  cat("  WARNING: Cannot extract df_plot and no backup exists.\n")
  cat("  Re-sync fig_variation_explained_tcgacptac from HPC and rebuild.\n")
  df_plot <- NULL
}

if (!is.null(df_plot)) {
  new_plot <- ggplot(df_plot,
         aes(x = variance_explained, y = delta_loglik,
             label = factor_label, color = method)) +
    geom_point(size = 4) +
    geom_text_repel(
      size = 4, max.overlaps = Inf, box.padding = 0.6,
      point.padding = 0.4, segment.size = 0.3, force = 2
    ) +
    scale_color_manual(values = c("NMF" = "red", "DeSurv" = "blue")) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Conditional variance explained\n(semi-partial R\u00b2)",
      y = expression(atop(Delta ~ "partial log-likelihood",
                          "(full vs. k-1 factor model)")),
      color = "Method"
    ) +
    theme_classic(base_size = 10)

  save_obj(new_plot, "fig_variation_explained_tcgacptac")
} else {
  cat("  SKIPPED: fig_variation_explained_tcgacptac not rebuilt\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. BO heatmap — rebuild from data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n4. BO heatmap (rebuild from data)\n")

obj <- read_obj("fig_bo_heat_tcgacptac")
curve_data <- obj$data

new_bo <- ggplot(curve_data, aes(x = k, y = alpha, fill = mean)) +
  geom_tile(color = NA) +
  scale_x_continuous(breaks = seq(2, 12, by = 1)) +
  scale_fill_viridis_c(
    name = "CV C-index",
    option = "D",
    guide = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.4, "cm"))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
  ) +
  labs(x = "Factorization rank (k)", y = "Supervision strength") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(color = "black")
  )

save_obj(new_bo, "fig_bo_heat_tcgacptac")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Forest plot — rebuild from data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n5. Forest plot (rebuild from data)\n")

obj <- read_obj("fig_hr_forest_tcgacptac")
df_forest <- obj$data

pd <- position_dodge(width = 0.6)

new_forest <- ggplot(df_forest,
       aes(x = HR, y = factor_name, color = dataset, group = dataset)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.5, color = "grey60") +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.25, linewidth = 0.8, position = pd
  ) +
  geom_point(size = 2.8, position = pd) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4),
    labels = c("0.5", "1", "2", "4")
  ) +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(
    xlim = c(min(df_forest$lower), max(df_forest$upper) * 1.4)
  ) +
  labs(x = "Hazard ratio (95% CI)", y = NULL) +
  facet_wrap(~method)

save_obj(new_forest, "fig_hr_forest_tcgacptac")

# ─────────────────────────────────────────────────────────────────────────────
# 6. KM survival curves — rebuild from curve data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n6. KM survival curves (rebuild from curve data)\n")

rebuild_km <- function(name) {
  obj <- read_obj(name)

  # Skip if already rebuilt (rebuilt plots have 3 layers: step, point, text)
  if (length(obj$plot$layers) <= 3 && !is.null(obj$table)) {
    cat("  Already rebuilt:", name, "\n")
    return(invisible(NULL))
  }

  curve_df <- obj$data.survplot
  table_df <- obj$data.survtable

  # Extract annotation text and parameters from old plot
  old_env <- obj$plot$plot_env
  pms <- get("pms", envir = old_env)
  break_by <- pms$break.time.by
  palette <- pms$palette
  legend_labs <- pms$legend.labs

  # Get HR annotation from layer 4
  annot_label <- obj$plot$layers[[4]]$aes_params$label
  annot_data <- obj$plot$layers[[4]]$data
  annot_x <- annot_data$x
  annot_y <- annot_data$y
  annot_hjust <- obj$plot$layers[[4]]$aes_params$hjust
  annot_size <- obj$plot$layers[[4]]$aes_params$size

  # Map factor to strata labels
  curve_df$group <- factor(curve_df$factor,
                           levels = c(0, 1),
                           labels = legend_labs)

  # Identify censoring rows
  censor_df <- curve_df[curve_df$n.censor > 0, ]

  x_max <- max(curve_df$time, na.rm = TRUE)

  # Build survival curve plot
  p <- ggplot(curve_df, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 0.8) +
    geom_point(data = censor_df, aes(x = time, y = surv, color = group),
               shape = 3, size = 2, show.legend = FALSE) +
    annotate("text", x = annot_x, y = annot_y,
             hjust = annot_hjust, size = annot_size,
             label = annot_label) +
    scale_color_manual(values = setNames(palette, legend_labs), name = NULL) +
    scale_x_continuous(breaks = seq(0, x_max, by = break_by)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Time (months)", y = "Survival probability") +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom")

  # Build risk table
  table_df$group <- factor(table_df$factor,
                           levels = c(0, 1),
                           labels = legend_labs)
  t <- ggplot(table_df, aes(x = time, y = group, label = n.risk)) +
    geom_text(size = 2.5) +
    scale_x_continuous(breaks = seq(0, x_max, by = break_by)) +
    labs(x = "Time (months)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 8, hjust = 0)
    ) +
    ggtitle("Number at risk")

  # Preserve the list structure expected by the Rmd
  obj$plot <- p
  obj$table <- t
  save_obj(obj, name)
}

rebuild_km("fig_median_survival_desurv_tcgacptac")
rebuild_km("fig_median_survival_std_desurvk_tcgacptac")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Simulation figures — rebuild from embedded data using sim_figs.R
# ─────────────────────────────────────────────────────────────────────────────
cat("\n7. Simulation figures (rebuild from embedded data)\n")

source("sim_figs.R")

sim_figs_list <- read_obj("sim_figs_by_scenario")

for (i in seq_along(sim_figs_list)) {
  entry <- sim_figs_list[[i]]
  sid <- entry$scenario_id
  aid <- entry$analysis_id
  cat(sprintf("  Element %d: %s / %s\n", i, sid, aid))

  # Extract data from old ggplot $data slots and rebuild
  results_data <- entry$cindex_box$data

  # Rebuild cindex_box
  sim_figs_list[[i]]$cindex_box <- plot_sim_metric_box(
    results_data, "cindex", "Test C-index", "C-index", base_size = 12
  )

  # Rebuild precision_box
  if (!is.null(entry$precision_box)) {
    sim_figs_list[[i]]$precision_box <- plot_sim_metric_box(
      results_data, "precision",
      "Best precision (mean across lethal factors)", "Precision",
      base_size = 12
    )
  }

  # Rebuild k_hist from its own data (has true_k columns)
  k_data <- entry$k_hist$data
  true_k_tbl <- if ("true_k.y" %in% names(k_data)) {
    unique(k_data[, c("scenario_panel", "true_k.y"), drop = FALSE]) |>
      setNames(c("scenario_panel", "true_k"))
  } else if ("true_k" %in% names(k_data)) {
    unique(k_data[, c("scenario_panel", "true_k"), drop = FALSE])
  } else {
    data.frame(scenario_panel = character(), true_k = integer())
  }
  true_k_tbl <- true_k_tbl[!is.na(true_k_tbl$true_k), , drop = FALSE]
  sim_figs_list[[i]]$k_hist <- plot_sim_k_hist(k_data, true_k_tbl, base_size = 12)

  # Rebuild precision_breakdown if present (mixed scenario)
  if (!is.null(entry$precision_breakdown)) {
    sim_figs_list[[i]]$precision_breakdown <- plot_mixed_precision_breakdown(
      results_data, base_size = 12
    )
  }

  # Rebuild matched_beta_box if present (mixed scenario)
  if (!is.null(entry$matched_beta_box)) {
    sim_figs_list[[i]]$matched_beta_box <- plot_matched_factor_beta(
      results_data, base_size = 12
    )
  }
}

save_obj(sim_figs_list, "sim_figs_by_scenario")

# ─────────────────────────────────────────────────────────────────────────────
# 8. NMF diagnostic plots (supplement) — rebuild from embedded data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n8. NMF diagnostic plots (rebuild from data)\n")

rebuild_nmf_diagnostic <- function(name) {
  obj <- read_obj(name)
  d <- obj$data
  metric <- unique(d$Measure)

  # Shorten silhouette.* facet labels before plot creation
  if (length(unique(d$variable)) > 1) {
    label_map <- c(
      "silhouette.coef"      = "Coeff.",
      "silhouette.basis"     = "Basis",
      "silhouette.consensus" = "Consens.",
      # Handle already-cleaned labels from prior rebuild
      "Coef"      = "Coeff.",
      "Basis"     = "Basis",
      "Consensus" = "Consens."
    )
    m <- match(d$variable, names(label_map))
    d$variable[!is.na(m)] <- label_map[m[!is.na(m)]]
  }

  p <- ggplot(d, aes(x = rank, y = value)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor.x = element_blank()) +
    scale_x_continuous(breaks = seq(2, 12, by = 2)) +
    labs(x = "Rank (k)")

  if (metric == "residuals") {
    p <- p + scale_y_continuous(
      labels = scales::label_number(scale = 1e-10, accuracy = 0.1),
      name = expression("Reconstruction error" ~ (x10^10))
    )
  } else {
    p <- p + labs(y = metric)
  }

  # Add facet if silhouette (multiple variables)
  if (length(unique(d$variable)) > 1) {
    p <- p + facet_wrap(~ variable)
  }

  save_obj(p, name)
}

for (nm in c("fig_residuals_tcgacptac", "fig_cophenetic_tcgacptac", "fig_silhouette_tcgacptac")) {
  if (file.exists(file.path(store_dir, nm))) {
    rebuild_nmf_diagnostic(nm)
  } else {
    cat("  SKIP (not in store):", nm, "\n")
  }
}

cat("\n=== Done. All main-text figures rebuilt. ===\n")
