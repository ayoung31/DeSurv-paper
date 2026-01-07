library(targets)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(tibble)

tar_load_globals(script = "_targets_sims.R")
tar_load(sim_results_table)
tar_load(sim_analysis_result)

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

# Metric definitions:
# - cindex: test concordance index computed on held-out samples.
# - precision/f1/concentration: mean of per-lethal-factor best values (precision, F1, recall) from top-gene recovery.
# - marker_ari: adjusted Rand index comparing true vs learned marker assignments;
#   genes are assigned to factor labels based on true/learned marker sets, using
#   the corresponding W matrices to break ties when needed.
# - marker_recall_all: global recall of true marker genes; union of learned top genes
#   divided by union of true marker genes across all programs.
# - purity_mean_nonzero: mean purity among learned factors with |beta| > 0; purity is
#   max overlap with any true marker set divided by the learned factor's top genes.
# - purity_mean: mean purity across all learned factors (including beta = 0).
# - n_learned_with_lethal_markers: number of learned factors that overlap at least one
#   true lethal marker set (nonzero survival effect in the simulator).
results <- sim_results_table %>%
  mutate(
    alpha_group = if_else(alpha == 0, "alpha=0", "alpha>0"),
    alpha_group = factor(alpha_group, levels = c("alpha=0", "alpha>0")),
    precision = map_dbl(lethal_factor_metrics, mean_lethal_metric, metric = "best_precision"),
    precision_nonzero = map_dbl(
      lethal_factor_metrics,
      mean_lethal_metric,
      metric = "best_precision_nonzero"
    ),
    f1 = map_dbl(lethal_factor_metrics, mean_lethal_metric, metric = "best_f1"),
    concentration = map_dbl(lethal_factor_metrics, mean_lethal_metric, metric = "concentration")
  )

plot_metric <- function(df, metric, title, ylab = NULL) {
  ggplot(df, aes(x = alpha_group, y = .data[[metric]], fill = alpha_group)) +
    geom_violin(width = 0.8, alpha = 0.6, trim = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
    facet_wrap(~ scenario, scales = "free_y") +
    labs(title = title, x = NULL, y = ylab %||% metric) +
    theme_bw() +
    theme(legend.position = "none")
}

plot_cindex <- plot_metric(results, "cindex", "Test c-index by alpha", "C-index")
plot_precision <- plot_metric(
  results,
  "precision",
  "Best precision (mean across lethal factors)",
  "Precision"
)

extra_metric_labels <- c(
  marker_ari = "Marker ARI",
  marker_recall_all = "Global marker recall",
  purity_mean = "Purity mean (all factors)",
  purity_mean_nonzero = "Purity mean (nonzero)",
  n_learned_with_lethal_markers = "Learned factors with lethal markers",
  f1 = "Best F1 (mean across lethal factors)"
)
extra_metric_labels <- extra_metric_labels[names(extra_metric_labels) %in% names(results)]
extra_metric_plots <- purrr::imap(
  extra_metric_labels,
  function(title, metric) plot_metric(results, metric, title)
)

analysis_tbl <- dplyr::bind_rows(sim_analysis_result)
true_k_tbl <- analysis_tbl %>%
  mutate(true_k = purrr::map_int(truth, ~ get_simulation_k(list(simulation = .x)))) %>%
  group_by(scenario_id, scenario) %>%
  summarise(true_k = dplyr::first(true_k), .groups = "drop")

k_plot_data <- results %>%
  filter(!is.na(k)) %>%
  left_join(true_k_tbl, by = c("scenario_id", "scenario"))

# k: selected latent factor count; dashed line shows true k from the simulator.
plot_k <- ggplot(k_plot_data, aes(x = alpha_group, y = k, fill = alpha_group)) +
  geom_violin(width = 0.8, alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
  geom_hline(
    data = true_k_tbl,
    aes(yintercept = true_k),
    linetype = "dashed",
    color = "black"
  ) +
  facet_wrap(~ scenario, scales = "free_y") +
  labs(title = "Selected k by alpha", x = NULL, y = "Selected k", fill = NULL) +
  theme_bw() +
  theme(legend.position = "none")

sim_plots <- c(
  list(cindex = plot_cindex, precision = plot_precision),
  extra_metric_plots,
  list(k_hist = plot_k)
)
