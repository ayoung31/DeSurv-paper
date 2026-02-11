#!/usr/bin/env Rscript
# Regenerate simulation and BO heatmap figures with updated labels
# Run from project root: Rscript regenerate_figures.R

suppressPackageStartupMessages({
  library(targets)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(viridis)
  library(grid)
})

STORE <- "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"

read_target <- function(name, store = STORE) {
  tryCatch(
    tar_read_raw(name, store = store),
    error = function(e) {
      obj_path <- file.path(store, "objects", name)
      if (file.exists(obj_path)) {
        cat("  (reading", name, "from RDS directly)\n")
        readRDS(obj_path)
      } else {
        stop("Cannot find target: ", name, call. = FALSE)
      }
    }
  )
}

cat("=== Loading source files ===\n")
source("sim_figs.R")
source("R/figure_targets.R")

# -------------------------------------------------------
# 1. Regenerate simulation figures with new scenario labels
# -------------------------------------------------------
cat("\n=== Regenerating simulation figures ===\n")
sim_results_table <- read_target("sim_results_table")
cat("  Loaded sim_results_table:", nrow(sim_results_table), "rows\n")

sim_figs <- build_sim_figs_by_scenario(sim_results_table)

scenario_ids <- sapply(sim_figs, function(x) x$scenario_id)
analysis_ids <- sapply(sim_figs, function(x) x$analysis_id)
cat("  Scenarios:", paste(unique(scenario_ids), collapse = ", "), "\n")
cat("  Analyses:", paste(unique(analysis_ids), collapse = ", "), "\n")

# Save to figures/sim/ directory
dir.create("figures/sim", recursive = TRUE, showWarnings = FALSE)
paths <- save_sim_figs_by_scenario(sim_figs, "figures/sim", SIM_FIGURE_CONFIGS)
cat("  Saved simulation figures:\n")
for (p in paths) cat("    ", p, "\n")

# Also save the primary scenario plots used in the paper (Fig 3)
alt <- which(scenario_ids == "R0_easy" & analysis_ids == "bo_tune_ntop")
if (length(alt)) {
  alt_plots <- sim_figs[[alt]]

  fig_sim <- plot_grid(
    alt_plots$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(.5, 1)),
    alt_plots$precision_box + labs(title = NULL),
    ncol = 2, labels = c("A", "B")
  )
  ggsave("figures/fig_sim_primary.pdf", fig_sim, width = 6.5, height = 3, bg = "white")
  cat("  Saved figures/fig_sim_primary.pdf\n")
}

# -------------------------------------------------------
# 2. Regenerate BO heatmap with new axis labels
# -------------------------------------------------------
cat("\n=== Regenerating BO heatmap ===\n")

# Load the GP curve data
desurv_bo_results <- read_target("desurv_bo_results_tcgacptac")
tar_params_best <- read_target("tar_params_best_tcgacptac")

curve <- extract_gp_curve(desurv_bo_results, tar_params_best)
cat("  GP curve data:", nrow(curve), "rows\n")

fig_bo_heat <- ggplot(curve, aes(x = k, y = alpha, fill = mean)) +
  geom_tile(color = NA) +
  scale_x_continuous(
    breaks = function(x) seq(2, 12, by = 1)
  ) +
  scale_fill_viridis_c(
    name = "CV C-index",
    option = "D",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth  = unit(0.4, "cm")
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
  ) +
  labs(
    x = "Factorization rank (k)",
    y = "Supervision strength"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(color = "black")
  )

dir.create("figures/panels", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/panels/fig_bo_heat_tcgacptac.pdf", fig_bo_heat, width = 4, height = 3.5, bg = "white")
cat("  Saved figures/panels/fig_bo_heat_tcgacptac.pdf\n")

# Also do bladder
desurv_bo_results_bl <- read_target("desurv_bo_results_bladder")
tar_params_best_bl <- read_target("tar_params_best_bladder")

curve_bl <- extract_gp_curve(desurv_bo_results_bl, tar_params_best_bl)

fig_bo_heat_bl <- ggplot(curve_bl, aes(x = k, y = alpha, fill = mean)) +
  geom_tile(color = NA) +
  scale_x_continuous(
    breaks = function(x) seq(2, 12, by = 1)
  ) +
  scale_fill_viridis_c(
    name = "CV C-index",
    option = "D",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth  = unit(0.4, "cm")
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
  ) +
  labs(
    x = "Factorization rank (k)",
    y = "Supervision strength"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(color = "black")
  )

ggsave("figures/panels/fig_bo_heat_bladder.pdf", fig_bo_heat_bl, width = 4, height = 3.5, bg = "white")
cat("  Saved figures/panels/fig_bo_heat_bladder.pdf\n")

cat("\n=== Done ===\n")
