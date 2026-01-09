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

build_sim_fig_data <- function(sim_results_table, sim_analysis_result) {
  results <- sim_results_table %>%
    dplyr::mutate(
      method = dplyr::if_else(alpha == 0, "alpha=0", "DeSurv"),
      method = factor(method, levels = c("DeSurv", "alpha=0")),
      precision = purrr::map_dbl(
        lethal_factor_metrics,
        mean_lethal_metric,
        metric = "best_precision"
      ),
      scenario = factor(scenario, levels = unique(scenario))
    )

  analysis_tbl <- dplyr::bind_rows(sim_analysis_result)
  true_k_tbl <- if (nrow(analysis_tbl)) {
    analysis_tbl %>%
      dplyr::mutate(true_k = purrr::map_int(truth, ~ get_simulation_k(list(simulation = .x)))) %>%
      dplyr::group_by(scenario_id, scenario) %>%
      dplyr::summarise(true_k = dplyr::first(true_k), .groups = "drop")
  } else {
    tibble::tibble(scenario_id = character(), scenario = character(), true_k = integer())
  }

  k_plot_data <- results %>%
    dplyr::filter(!is.na(k)) %>%
    dplyr::left_join(true_k_tbl, by = c("scenario_id", "scenario")) %>%
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
    ggplot2::facet_grid(scenario ~ method) +
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

  ggplot2::ggplot(plot_data, ggplot2::aes(x = method, y = .data[[metric]], fill = method)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = method),
      width = 0.15,
      size = 0.8,
      alpha = 0.35
    ) +
    ggplot2::facet_wrap(~ scenario) +
    ggplot2::scale_fill_manual(values = sim_method_colors) +
    ggplot2::scale_color_manual(values = sim_method_colors, guide = "none") +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = ylab,
      fill = "Method"
    ) +
    sim_pub_theme(base_size = base_size)
}

build_sim_figs <- function(sim_results_table, sim_analysis_result, base_size = 12) {
  plot_data <- build_sim_fig_data(sim_results_table, sim_analysis_result)

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
    )
  )
}

save_sim_plot <- function(plot, path, width = 6.5, height = 4.5) {
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

if (interactive()) {
  library(targets)
  tar_load_globals(script = "_targets_sims.R")
  tar_load(sim_results_table)
  tar_load(sim_analysis_result)
  sim_plots <- build_sim_figs(sim_results_table, sim_analysis_result)
}
