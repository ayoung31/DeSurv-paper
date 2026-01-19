suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

sim_method_colors <- c(
  "DeSurv" = "#1f78b4",
  "alpha=0" = "#e31a1c"
)

sim_pub_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.9, "lines")
    )
}

normalize_scenario_id <- function(scenario_id) {
  if (is.null(scenario_id) || !length(scenario_id)) {
    return(NULL)
  }
  if (length(scenario_id) > 1) {
    stop(
      "Expected a single scenario_id, got: ",
      paste(scenario_id, collapse = ", "),
      call. = FALSE
    )
  }
  as.character(scenario_id[[1]])
}

mean_lethal_metric <- function(tbl, metric) {
  if (is.null(tbl) || !nrow(tbl)) {
    return(NA_real_)
  }
  vals <- tbl[[metric]]
  if (!length(vals)) {
    return(NA_real_)
  }
  mean(vals, na.rm = TRUE)
}

build_sim_fig_data <- function(sim_results_table,
                               analysis_id_base = NULL,
                               scenario_id = NULL) {
  results <- sim_results_table
  scen_id <- normalize_scenario_id(scenario_id)
  anal_id = analysis_id_base
  if (!is.null(analysis_id_base)) {
    results <- results %>%
      dplyr::filter(analysis_id_base == anal_id)
  }
  if (!is.null(scenario_id)) {
    results <- results %>%
      dplyr::filter(scenario_id == scen_id)
  }

  scenario_col <- if ("scenario_id" %in% names(results)) "scenario_id" else "scenario"
  results <- results %>%
    dplyr::mutate(
      method = dplyr::if_else(alpha == 0, "alpha=0", "DeSurv"),
      method = factor(method, levels = c("DeSurv", "alpha=0")),
      precision = purrr::map_dbl(
        lethal_factor_metrics,
        mean_lethal_metric,
        metric = "best_precision"
      ),
      scenario_panel = .data[[scenario_col]]
    ) %>%
    dplyr::mutate(
      scenario_panel = factor(scenario_panel, levels = unique(scenario_panel))
    )

  true_k_tbl <- if ("true_k" %in% names(results) && nrow(results)) {
    results %>%
      dplyr::filter(!is.na(true_k)) %>%
      dplyr::group_by(scenario_panel) %>%
      dplyr::summarise(true_k = dplyr::first(true_k), .groups = "drop")
  } else {
    tibble::tibble(scenario_panel = character(), true_k = integer())
  }

  k_plot_data <- results %>%
    dplyr::filter(!is.na(k)) %>%
    dplyr::left_join(true_k_tbl, by = "scenario_panel") %>%
    dplyr::mutate(k = as.integer(k))

  list(results = results, true_k_tbl = true_k_tbl, k_plot_data = k_plot_data)
}

plot_sim_k_hist <- function(k_plot_data, true_k_tbl, base_size = 12) {
  ggplot2::ggplot(k_plot_data, ggplot2::aes(x = k, fill = method)) +
    ggplot2::geom_histogram(
      binwidth = 1,
      boundary = 0.5,
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      data = true_k_tbl,
    ggplot2::aes(xintercept = true_k),
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
    ggplot2::facet_grid(scenario_panel ~ method) +
    ggplot2::scale_fill_manual(values = sim_method_colors, guide = "none") +
    ggplot2::scale_x_continuous(breaks = sort(unique(k_plot_data$k))) +
    ggplot2::labs(
      title = "Selected k distribution",
      x = "Selected k",
      y = "Count"
    ) +
    sim_pub_theme(base_size = base_size)
}

plot_sim_metric_box <- function(results, metric, title, ylab, base_size = 12) {
  plot_data <- results %>%
    dplyr::filter(!is.na(.data[[metric]]), !is.na(method))
  if(nrow(plot_data) > 0){
    ggplot2::ggplot(plot_data, ggplot2::aes(x = method, y = .data[[metric]], fill = method)) +
      ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = method),
        width = 0.15,
        size = 0.8,
        alpha = 0.35
      ) +
      ggplot2::facet_wrap(~ scenario_panel) +
      ggplot2::scale_fill_manual(values = sim_method_colors) +
      ggplot2::scale_color_manual(values = sim_method_colors, guide = "none") +
      ggplot2::labs(
        title = title,
        x = NULL,
        y = ylab,
        fill = "Method"
      ) +
      sim_pub_theme(base_size = base_size)
  }else{
    NULL
  }
  
}

build_sim_figs <- function(sim_results_table,
                           analysis_id_base = NULL,
                           scenario_id = NULL,
                           base_size = 12) {
  plot_data <- build_sim_fig_data(
    sim_results_table,
    analysis_id_base = analysis_id_base,
    scenario_id = scenario_id
  )

  list(
    k_hist = plot_sim_k_hist(
      plot_data$k_plot_data,
      plot_data$true_k_tbl,
      base_size = base_size
    ),
    cindex_box = plot_sim_metric_box(
      plot_data$results,
      "cindex",
      "Test C-index",
      "C-index",
      base_size = base_size
    ),
    precision_box = plot_sim_metric_box(
      plot_data$results,
      "precision",
      "Best precision (mean across lethal factors)",
      "Precision",
      base_size = base_size
    ),
    analysis_id = analysis_id_base,
    scenario_id = scenario_id
  )
}

save_sim_plot <- function(plot, path, width = 6.5, height = 4.5) {
  if(!is.null(plot)){
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    path
  }
  
}

sim_fig_basename <- function(filename) {
  tools::file_path_sans_ext(filename)
}

sim_fig_suffix <- function(scenario_id) {
  if (is.null(scenario_id) || is.na(scenario_id) || !nzchar(scenario_id)) {
    return("combined")
  }
  scenario_id
}

save_sim_figs <- function(plots,
                          sim_dir,
                          figure_configs) {
  scenario_id = plots$scenario_id
  analysis_id = plots$analysis_id
  out_dir <- sim_dir
  k_hist_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_selected_k_hist", scenario_id,analysis_id)
  )
  cindex_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_cindex_boxplot",scenario_id,analysis_id)
  )
  precision_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_precision_boxplot",scenario_id,analysis_id)
  )

  c(
    save_sim_plot(
      plots$k_hist,
      k_hist_path
    ),
    save_sim_plot(
      plots$cindex_box,
      cindex_path
    ),
    save_sim_plot(
      plots$precision_box,
      precision_path
    )
  )
}

build_sim_figs_by_scenario <- function(sim_results_table,
                                      base_size = 12) {
  sim_results_table$analysis_id_base = sub("_alpha0","",sim_results_table$analysis_id)
  scenario_ids <- sim_results_table$scenario_id
  scenario_ids <- as.character(scenario_ids)
  scenario_ids <- unique(scenario_ids)
  scenario_ids <- scenario_ids[!is.na(scenario_ids) & nzchar(scenario_ids)]
  scenario_ids <- sort(scenario_ids)
  
  analysis_ids = sim_results_table$analysis_id
  analysis_ids = as.character(analysis_ids)
  analysis_ids = unique(analysis_ids)
  analysis_ids = analysis_ids[!is.na(analysis_ids) & nzchar(analysis_ids)]
  analysis_ids = sort(analysis_ids)
  
  analysis_ids_base = sub("_alpha0","",analysis_ids)
  analysis_ids_base = unique(analysis_ids_base)
  
  plots = NULL
  for(sid in scenario_ids){
    for(aid in analysis_ids_base){
      p <- build_sim_figs(
        sim_results_table,
        analysis_id_base = aid,
        scenario_id = sid,
        base_size = base_size
      )
      plots = append(plots,list(p))
    }
  }
  
  plots
}

save_sim_figs_by_scenario = function(sim_figs,sim_dir,figure_configs){
  paths = vector()
  for(i in 1:length(sim_figs)){
    paths = c(paths,save_sim_figs(sim_figs[[i]],sim_dir,figure_configs))
  }
  paths
}

if (interactive()) {
  library(targets)
  tar_load_globals(script = "_targets_sims.R")
  tar_load(sim_results_table)
  sim_plots <- build_sim_figs(sim_results_table)
}
