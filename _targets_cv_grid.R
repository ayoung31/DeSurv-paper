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
CV_GRID_NSTARTS_FULL <- 100
CV_GRID_METHOD_TRANSFORM <- "rank"
CV_GRID_Z_CUTPOINTS <- seq(-2.0, 2.0, by = 0.2)  # 21 z-score cutpoints

# ------ Source helpers ------
source("R/cv_grid_helpers.R")
source("R/load_data.R")
source("R/load_data_internal.R")
source("R/targets_config.R")

# ------ Validation datasets ------
CV_GRID_VAL_DATASETS <- c("Dijk", "Moffitt_GEO_array",
                           "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

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

cv_grid_full_controller <- crew_controller_slurm(
  name = "cv_grid_full",
  workers = 440,
  seconds_idle = 300,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = CV_GRID_NSTARTS_FULL,
    memory_gigabytes_required = 32,
    time_minutes = 480,
    log_error = "logs/cv_grid_full_%A.err",
    log_output = "logs/cv_grid_full_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

cv_grid_val_controller <- crew_controller_slurm(
  name = "cv_grid_val",
  workers = 440,
  seconds_idle = 300,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_required = 4,
    time_minutes = 30,
    log_error = "logs/cv_grid_val_%A.err",
    log_output = "logs/cv_grid_val_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

active_controller <- crew_controller_group(
  default_controller,
  cv_grid_controller,
  cv_grid_full_controller,
  cv_grid_val_controller
)

# ------ Global options ------
tar_option_set(
  packages = c("DeSurv", "dplyr", "purrr", "tibble", "survival",
                "ggplot2", "survminer", "cowplot"),
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

  # Target N1: Evaluate z-score cutpoints for each grid point
  tar_target(
    cv_grid_cutpoint_eval,
    evaluate_cutpoint_zscores(
      cv_grid_result,
      z_grid = CV_GRID_Z_CUTPOINTS
    ),
    pattern = map(cv_grid_result),
    iteration = "list"
  ),

  # Target N2: Aggregate cutpoint evaluations across all grid points
  tar_target(
    cv_grid_cutpoint_summary,
    {
      all_evals <- dplyr::bind_rows(purrr::compact(cv_grid_cutpoint_eval))
      all_evals |>
        dplyr::mutate(abs_logrank_z = abs(logrank_z)) |>
        dplyr::group_by(k, alpha, ntop, z_cutpoint) |>
        dplyr::summarise(
          mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
          se_cindex_dichot = sd(cindex_dichot, na.rm = TRUE) /
            sqrt(sum(!is.na(cindex_dichot))),
          mean_abs_logrank_z = mean(abs_logrank_z, na.rm = TRUE),
          se_abs_logrank_z = sd(abs_logrank_z, na.rm = TRUE) /
            sqrt(sum(!is.na(abs_logrank_z))),
          .groups = "drop"
        )
    }
  ),

  # Target N3: Save cutpoint summary to CSV (with absolute cutpoints from full fits)
  tar_target(
    cv_grid_cutpoint_summary_file,
    {
      fit_stats <- aggregate_cv_grid_fit_results(cv_grid_fit) |>
        dplyr::select(k, alpha, ntop, lp_mean, lp_sd)
      out <- cv_grid_cutpoint_summary |>
        dplyr::left_join(fit_stats, by = c("k", "alpha", "ntop")) |>
        dplyr::mutate(cutpoint_abs = z_cutpoint * lp_sd + lp_mean) |>
        dplyr::select(k, alpha, ntop, z_cutpoint, cutpoint_abs,
                      mean_cindex_dichot, se_cindex_dichot,
                      mean_abs_logrank_z, se_abs_logrank_z)
      path <- "results/cv_grid/cv_grid_cutpoint_summary.csv"
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(out, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  # Target N4: Select optimal cutpoint per grid point
  tar_target(
    cv_grid_optimal_cutpoint,
    select_optimal_cutpoint(cv_grid_cutpoint_summary)
  ),

  # Target 4: Aggregate all CV results with optimal cutpoint info
  tar_target(
    cv_grid_summary,
    {
      cv_tbl <- aggregate_cv_grid_results(cv_grid_result)
      opt_tbl <- cv_grid_optimal_cutpoint
      dplyr::left_join(cv_tbl, opt_tbl, by = c("k", "alpha", "ntop"))
    }
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
  ),

  # Target 6: Full model fit per grid point (440 parallel jobs, 100 inits each)
  tar_target(
    cv_grid_fit,
    {
      params <- cv_grid_params
      fixed <- CV_GRID_FIXED_PARAMS
      fixed$ntop <- params$ntop
      ntop_match <- if (is.null(params$ntop)) NA_real_ else as.numeric(params$ntop)
      if (is.na(ntop_match)) {
        opt_row <- cv_grid_optimal_cutpoint |>
          dplyr::filter(k == params$k, alpha == params$alpha, is.na(ntop))
      } else {
        opt_row <- cv_grid_optimal_cutpoint |>
          dplyr::filter(k == params$k, alpha == params$alpha, ntop == ntop_match)
      }
      opt_z <- if (nrow(opt_row) == 1) opt_row$optimal_z_cutpoint else NA_real_
      fit_grid_point(
        data = cv_grid_data,
        k = params$k,
        alpha = params$alpha,
        fixed_params = fixed,
        n_starts = CV_GRID_NSTARTS_FULL,
        seed = CV_GRID_SEED,
        optimal_z_cutpoint = opt_z,
        verbose = TRUE
      )
    },
    pattern = map(cv_grid_params),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv_grid_full")
    )
  ),

  # Target 7: Aggregate fit metadata into summary table
  tar_target(
    cv_grid_fit_summary,
    aggregate_cv_grid_fit_results(cv_grid_fit)
  ),

  # Target 8: Save fit summary to CSV
  tar_target(
    cv_grid_fit_summary_file,
    {
      path <- "results/cv_grid/cv_grid_fit_summary.csv"
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(cv_grid_fit_summary, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  # Target 9: Load and preprocess external validation datasets
  tar_target(
    cv_grid_val_data,
    {
      train_genes <- rownames(cv_grid_data$ex)
      transform_target <- cv_grid_data$transform_target

      val_list <- lapply(CV_GRID_VAL_DATASETS, function(ds_name) {
        raw <- load_data(ds_name)
        preprocess_validation_data(
          dataset = raw,
          genes = train_genes,
          method_trans_train = CV_GRID_METHOD_TRANSFORM,
          dataname = ds_name,
          transform_target = transform_target,
          zero_fill_missing = TRUE
        )
      })
      names(val_list) <- CV_GRID_VAL_DATASETS
      val_list
    }
  ),

  # Target 9b: Merge PACA_AU_array + PACA_AU_seq for survival validation
  tar_target(
    cv_grid_val_data_surv,
    merge_paca_au_datasets(cv_grid_val_data)
  ),

  # Target 10: Validate each grid fit against external datasets (440 parallel jobs)
  tar_target(
    cv_grid_val_result,
    validate_grid_point(cv_grid_fit, cv_grid_val_data_surv),
    pattern = map(cv_grid_fit),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv_grid_val")
    )
  ),

  # Target 11: Aggregate all validation results
  tar_target(
    cv_grid_val_summary,
    aggregate_cv_grid_val_results(cv_grid_val_result)
  ),

  # Target 12: Save validation summary to CSV
  tar_target(
    cv_grid_val_summary_file,
    {
      path <- "results/cv_grid/cv_grid_val_summary.csv"
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(cv_grid_val_summary, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  # Target 13: Select best alpha per (k, ntop) using 3 criteria
  tar_target(
    cv_grid_best_alpha,
    select_best_alpha_per_k(cv_grid_summary)
  ),

  # Target 14: Save best alpha table to CSV
  tar_target(
    cv_grid_best_alpha_file,
    {
      path <- "results/cv_grid/cv_grid_best_alpha.csv"
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      write.csv(cv_grid_best_alpha, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  # Target 15a: Extract alpha=0 (standard NMF) combos for all k × ntop
  tar_target(
    cv_grid_alpha0,
    get_alpha0_combos(cv_grid_summary)
  ),

  # Target 15: Cutpoint evaluation plots for each selected (k, alpha, ntop, method)
  tar_target(
    cv_grid_cutpoint_plots,
    {
      out_dir <- "figures/cv_grid/cutpoint"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(nrow(cv_grid_best_alpha))

      for (i in seq_len(nrow(cv_grid_best_alpha))) {
        row <- cv_grid_best_alpha[i, ]
        k_i <- row$k
        alpha_i <- row$best_alpha
        ntop_i <- row$ntop
        method_i <- row$selection_method
        opt_z <- row$optimal_z_cutpoint

        # Filter cutpoint summary to this (k, alpha, ntop)
        if (is.na(ntop_i)) {
          cp_data <- cv_grid_cutpoint_summary |>
            dplyr::filter(k == k_i, alpha == alpha_i, is.na(ntop))
        } else {
          cp_data <- cv_grid_cutpoint_summary |>
            dplyr::filter(k == k_i, alpha == alpha_i, ntop == ntop_i)
        }

        if (nrow(cp_data) == 0) {
          paths[i] <- NA_character_
          next
        }

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)
        fname <- sprintf("cutpoint_k%d_alpha%.2f_ntop%s_%s.pdf",
                          k_i, alpha_i, ntop_label, method_i)
        fpath <- file.path(out_dir, fname)

        p <- plot_cutpoint_curves(cp_data, k_i, alpha_i, ntop_i,
                                   optimal_z = opt_z)
        ggplot2::ggsave(fpath, p, width = 10, height = 4)
        paths[i] <- fpath
      }

      paths[!is.na(paths)]
    },
    format = "file"
  ),

  # Target 16: Training KM plots for each selected combo
  tar_target(
    cv_grid_km_train_plots,
    {
      out_dir <- "figures/cv_grid/km_train"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(nrow(cv_grid_best_alpha))

      for (i in seq_len(nrow(cv_grid_best_alpha))) {
        row <- cv_grid_best_alpha[i, ]
        k_i <- row$k
        alpha_i <- row$best_alpha
        ntop_i <- row$ntop
        method_i <- row$selection_method

        # Find the matching fit entry
        fit_entry <- NULL
        for (fe in cv_grid_fit) {
          if (is.null(fe)) next
          fe_ntop <- if (is.null(fe$ntop)) NA_real_ else as.numeric(fe$ntop)
          if (fe$k == k_i && fe$alpha == alpha_i &&
              (is.na(ntop_i) && is.na(fe_ntop) ||
               !is.na(ntop_i) && !is.na(fe_ntop) && fe_ntop == ntop_i)) {
            fit_entry <- fe
            break
          }
        }

        if (is.null(fit_entry)) {
          paths[i] <- NA_character_
          next
        }

        # Override z_cutpoint with the selection-method-specific value
        fit_entry$z_cutpoint <- row$optimal_z_cutpoint

        p <- plot_km_training(fit_entry, cv_grid_data)
        if (is.null(p)) {
          paths[i] <- NA_character_
          next
        }

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)
        fname <- sprintf("km_train_k%d_alpha%.2f_ntop%s_%s.pdf",
                          k_i, alpha_i, ntop_label, method_i)
        fpath <- file.path(out_dir, fname)
        ggplot2::ggsave(fpath, p, width = 7, height = 5)
        paths[i] <- fpath
      }

      paths[!is.na(paths)]
    },
    format = "file"
  ),

  # Target 17: Validation KM plots for each combo × dataset (+ pooled)
  tar_target(
    cv_grid_km_val_plots,
    {
      out_dir <- "figures/cv_grid/km_val"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ds_names <- names(cv_grid_val_data_surv)
      paths <- character(0)

      for (i in seq_len(nrow(cv_grid_best_alpha))) {
        row <- cv_grid_best_alpha[i, ]
        k_i <- row$k
        alpha_i <- row$best_alpha
        ntop_i <- row$ntop
        method_i <- row$selection_method

        # Find the matching fit entry
        fit_entry <- NULL
        for (fe in cv_grid_fit) {
          if (is.null(fe)) next
          fe_ntop <- if (is.null(fe$ntop)) NA_real_ else as.numeric(fe$ntop)
          if (fe$k == k_i && fe$alpha == alpha_i &&
              (is.na(ntop_i) && is.na(fe_ntop) ||
               !is.na(ntop_i) && !is.na(fe_ntop) && fe_ntop == ntop_i)) {
            fit_entry <- fe
            break
          }
        }

        if (is.null(fit_entry)) next
        fit_entry$z_cutpoint <- row$optimal_z_cutpoint

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)

        # Per-dataset KM plots
        for (ds_name in ds_names) {
          p <- plot_km_validation(fit_entry, cv_grid_val_data_surv[[ds_name]],
                                   ds_name)
          if (is.null(p)) next

          fname <- sprintf("km_val_k%d_alpha%.2f_ntop%s_%s_%s.pdf",
                            k_i, alpha_i, ntop_label, method_i, ds_name)
          fpath <- file.path(out_dir, fname)
          ggplot2::ggsave(fpath, p, width = 7, height = 5)
          paths <- c(paths, fpath)
        }

        # Pooled KM plot
        p_pooled <- plot_km_validation_pooled(fit_entry, cv_grid_val_data_surv)
        if (!is.null(p_pooled)) {
          fname <- sprintf("km_val_k%d_alpha%.2f_ntop%s_%s_pooled.pdf",
                            k_i, alpha_i, ntop_label, method_i)
          fpath <- file.path(out_dir, fname)
          ggplot2::ggsave(fpath, p_pooled, width = 7, height = 5)
          paths <- c(paths, fpath)
        }
      }

      paths
    },
    format = "file"
  ),

  # Target 18: Cutpoint evaluation plots for alpha=0 (standard NMF) baseline
  tar_target(
    cv_grid_cutpoint_plots_alpha0,
    {
      out_dir <- "figures/cv_grid/cutpoint"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(nrow(cv_grid_alpha0))

      for (i in seq_len(nrow(cv_grid_alpha0))) {
        row <- cv_grid_alpha0[i, ]
        k_i <- row$k
        ntop_i <- row$ntop
        opt_z <- row$optimal_z_cutpoint

        if (is.na(ntop_i)) {
          cp_data <- cv_grid_cutpoint_summary |>
            dplyr::filter(k == k_i, alpha == 0, is.na(ntop))
        } else {
          cp_data <- cv_grid_cutpoint_summary |>
            dplyr::filter(k == k_i, alpha == 0, ntop == ntop_i)
        }

        if (nrow(cp_data) == 0) {
          paths[i] <- NA_character_
          next
        }

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)
        fname <- sprintf("cutpoint_k%d_alpha0.00_ntop%s_alpha0.pdf",
                          k_i, ntop_label)
        fpath <- file.path(out_dir, fname)

        p <- plot_cutpoint_curves(cp_data, k_i, 0, ntop_i,
                                   optimal_z = opt_z)
        ggplot2::ggsave(fpath, p, width = 10, height = 4)
        paths[i] <- fpath
      }

      paths[!is.na(paths)]
    },
    format = "file"
  ),

  # Target 19: Training KM plots for alpha=0 baseline
  tar_target(
    cv_grid_km_train_plots_alpha0,
    {
      out_dir <- "figures/cv_grid/km_train"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(nrow(cv_grid_alpha0))

      for (i in seq_len(nrow(cv_grid_alpha0))) {
        row <- cv_grid_alpha0[i, ]
        k_i <- row$k
        ntop_i <- row$ntop

        fit_entry <- NULL
        for (fe in cv_grid_fit) {
          if (is.null(fe)) next
          fe_ntop <- if (is.null(fe$ntop)) NA_real_ else as.numeric(fe$ntop)
          if (fe$k == k_i && fe$alpha == 0 &&
              (is.na(ntop_i) && is.na(fe_ntop) ||
               !is.na(ntop_i) && !is.na(fe_ntop) && fe_ntop == ntop_i)) {
            fit_entry <- fe
            break
          }
        }

        if (is.null(fit_entry)) {
          paths[i] <- NA_character_
          next
        }

        fit_entry$z_cutpoint <- row$optimal_z_cutpoint

        p <- plot_km_training(fit_entry, cv_grid_data)
        if (is.null(p)) {
          paths[i] <- NA_character_
          next
        }

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)
        fname <- sprintf("km_train_k%d_alpha0.00_ntop%s_alpha0.pdf",
                          k_i, ntop_label)
        fpath <- file.path(out_dir, fname)
        ggplot2::ggsave(fpath, p, width = 7, height = 5)
        paths[i] <- fpath
      }

      paths[!is.na(paths)]
    },
    format = "file"
  ),

  # Target 20: Validation KM plots for alpha=0 baseline (per-dataset + pooled)
  tar_target(
    cv_grid_km_val_plots_alpha0,
    {
      out_dir <- "figures/cv_grid/km_val"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ds_names <- names(cv_grid_val_data_surv)
      paths <- character(0)

      for (i in seq_len(nrow(cv_grid_alpha0))) {
        row <- cv_grid_alpha0[i, ]
        k_i <- row$k
        ntop_i <- row$ntop

        fit_entry <- NULL
        for (fe in cv_grid_fit) {
          if (is.null(fe)) next
          fe_ntop <- if (is.null(fe$ntop)) NA_real_ else as.numeric(fe$ntop)
          if (fe$k == k_i && fe$alpha == 0 &&
              (is.na(ntop_i) && is.na(fe_ntop) ||
               !is.na(ntop_i) && !is.na(fe_ntop) && fe_ntop == ntop_i)) {
            fit_entry <- fe
            break
          }
        }

        if (is.null(fit_entry)) next
        fit_entry$z_cutpoint <- row$optimal_z_cutpoint

        ntop_label <- if (is.na(ntop_i)) "ALL" else as.character(ntop_i)

        for (ds_name in ds_names) {
          p <- plot_km_validation(fit_entry, cv_grid_val_data_surv[[ds_name]],
                                   ds_name)
          if (is.null(p)) next

          fname <- sprintf("km_val_k%d_alpha0.00_ntop%s_alpha0_%s.pdf",
                            k_i, ntop_label, ds_name)
          fpath <- file.path(out_dir, fname)
          ggplot2::ggsave(fpath, p, width = 7, height = 5)
          paths <- c(paths, fpath)
        }

        p_pooled <- plot_km_validation_pooled(fit_entry, cv_grid_val_data_surv)
        if (!is.null(p_pooled)) {
          fname <- sprintf("km_val_k%d_alpha0.00_ntop%s_alpha0_pooled.pdf",
                            k_i, ntop_label)
          fpath <- file.path(out_dir, fname)
          ggplot2::ggsave(fpath, p_pooled, width = 7, height = 5)
          paths <- c(paths, fpath)
        }
      }

      paths
    },
    format = "file"
  )
)
