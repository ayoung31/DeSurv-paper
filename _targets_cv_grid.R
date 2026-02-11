# _targets_cv_grid.R
#
# Cross-validation grid search pipeline for k x alpha x ntop parameter space.
# Performs exhaustive CV over:
#   - k: 2, 3, 4, ..., 12 (11 values)
#   - alpha: 0, 0.05, 0.10, ..., 0.95 (20 values)
#   - ntop: NULL (all genes), 300 (2 values)
#   - Total: 440 parameter combinations
#
# Fixed parameters:
#   - lambda = 0.3
#   - nu = 0.05
#   - lambdaW = 0
#   - lambdaH = 0
#   - ngene = 3000 (for preprocessing)

library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(DeSurv)
})

# ------ Configuration ------
CV_GRID_K_VALUES <- 2:12
CV_GRID_ALPHA_VALUES <- seq(0, 0.95, by = 0.05)
CV_GRID_NTOP_VALUES <- list(NULL, 300)  # NULL = all genes, 300 = top 300
CV_GRID_FIXED_PARAMS <- list(
  lambda = 0.3,
  nu = 0.05,
  lambdaW = 0,
  lambdaH = 0
)
CV_GRID_NGENE <- 3000
CV_GRID_NFOLDS <- 5
CV_GRID_NSTARTS <- 30
CV_GRID_SEED <- 123
CV_GRID_METHOD_TRANSFORM <- "rank"

# ------ Source helpers ------
source("R/cv_grid_helpers.R")
source("R/load_data.R")
source("R/load_data_internal.R")

# ------ Slurm controllers ------
default_controller <- crew_controller_sequential()

cv_grid_controller <- crew_controller_slurm(
  name = "cv_grid",
  workers = 440,
  seconds_idle = 300,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = CV_GRID_NSTARTS,
    memory_gigabytes_required = 16,
    time_minutes = 480,
    log_error = "logs/cv_grid_%A.err",
    log_output = "logs/cv_grid_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

active_controller <- crew_controller_group(
  default_controller,
  cv_grid_controller
)

# ------ Global options ------
tar_option_set(
  packages = c("DeSurv", "dplyr", "purrr", "tibble", "survival"),
  format = "rds",
  controller = active_controller,
  error = "continue"
)

# ------ Target definitions ------
list(
  # Target 1: Generate the 440 (k, alpha, ntop) parameter combinations
  tar_target(
    cv_grid_params,
    create_cv_grid(
      k_values = CV_GRID_K_VALUES,
      alpha_values = CV_GRID_ALPHA_VALUES,
      ntop_values = CV_GRID_NTOP_VALUES
    ),
    iteration = "list"
  ),

  # Target 2: Load and preprocess TCGA+CPTAC training data with ngene=3000
  tar_target(
    cv_grid_data,
    {
      # Load raw TCGA and CPTAC data
      raw <- load_data(c("TCGA_PAAD", "CPTAC"))
      keep_idx <- raw$samp_keeps

      # Preprocess using DeSurv::preprocess_data
      processed <- DeSurv::preprocess_data(
        X = raw$ex[, keep_idx, drop = FALSE],
        y = raw$sampInfo$time[keep_idx],
        d = raw$sampInfo$event[keep_idx],
        dataset = raw$sampInfo$dataset[keep_idx],
        samp_keeps = NULL,
        ngene = CV_GRID_NGENE,
        method_trans_train = CV_GRID_METHOD_TRANSFORM,
        verbose = TRUE
      )

      # Return in the format expected by run_cv_grid_point
      list(
        ex = processed$ex,
        sampInfo = processed$sampInfo,
        transform_target = processed$transform_target
      )
    }
  ),

  # Target 3: Run CV for each (k, alpha, ntop) combination (440 parallel jobs)
  tar_target(
    cv_grid_result,
    {
      params <- cv_grid_params
      fixed <- CV_GRID_FIXED_PARAMS
      fixed$ntop <- params$ntop
      run_cv_grid_point(
        data = cv_grid_data,
        k = params$k,
        alpha = params$alpha,
        fixed_params = fixed,
        nfolds = CV_GRID_NFOLDS,
        n_starts = CV_GRID_NSTARTS,
        seed = CV_GRID_SEED,
        verbose = TRUE
      )
    },
    pattern = map(cv_grid_params),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv_grid")
    )
  ),

  # Target 4: Aggregate all 220 results into final summary table
  tar_target(
    cv_grid_summary,
    aggregate_cv_grid_results(cv_grid_result)
  ),

  # Target 5: Save summary table to CSV
  tar_target(
    cv_grid_summary_file,
    {
      path <- "results/cv_grid/cv_grid_summary.csv"
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(cv_grid_summary, path, row.names = FALSE)
      path
    },
    format = "file"
  )
)
