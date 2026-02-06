# targets_figure_configs.R - UNC Longleaf HPC mode
# Same as local version — figure configs don't need resource changes

targets_figure_configs <- function() {
  figures_dir <- "figures"
  list(
    figures_dir = figures_dir,
    panel_dir = file.path(figures_dir, "panels"),
    sim_dir = file.path(figures_dir, "sim"),
    paper_figure_keys = list(
      bo_key = "tcgacptac",
      run_key = "tcgacptac",
      val_key = "tcgacptac"
    ),
    sc_data_paths = list(
      all = "data/derv/Elyada_umap.Rds",
      caf = "data/derv/Elyada_caf_umap.Rds",
      tum = "data/derv/Elyada_PDAC_umap.Rds"
    ),
    sim_figures = list(
      k_hist = list(
        filename = "sim_selected_k_hist.pdf",
        width = 6.5,
        height = 4.5
      ),
      cindex_box = list(
        filename = "sim_cindex_boxplot.pdf",
        width = 6.0,
        height = 4.2
      ),
      precision_box = list(
        filename = "sim_precision_boxplot.pdf",
        width = 6.0,
        height = 4.2
      )
    )
  )
}
