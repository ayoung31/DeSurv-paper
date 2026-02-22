bo_bundle_target <- tar_target(
  tar_bo_bundle,
  list(
    label = bo_config$label,
    config = bo_config,
    data = tar_data,
    data_split = tar_data_split,
    data_filtered = tar_data_filtered,
    data_filtered_alpha0 = tar_data_filtered_alpha0,
    params_best = tar_params_best,
    params_best_alpha0 = tar_params_best_alpha0,
    ngene_value = tar_ngene_value,
    ngene_value_alpha0 = tar_ngene_value_alpha0,
    ntop_value = tar_ntop_value,
    ntop_value_alpha0 = tar_ntop_value_alpha0,
    lambdaW_value = tar_lambdaW_value,
    lambdaH_value = tar_lambdaH_value,
    lambdaW_value_alpha0 = tar_lambdaW_value_alpha0,
    lambdaH_value_alpha0 = tar_lambdaH_value_alpha0
  )
)

COMMON_DESURV_BO_TARGETS <- list(
  tar_target(
    desurv_bo_results,
    {
      bounds <- bo_config$desurv_bo_bounds
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ngene_config, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ntop_config, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdaw_config, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdah_config, "lambdaH_grid", log_scale = TRUE)

      bo_fixed <- list(n_starts = bo_config$ninit)
      if (!bo_config$tune_ngene) {
        bo_fixed$ngene <- bo_config$ngene_default
      }
      if (!bo_config$tune_ntop) {
        bo_fixed$ntop <- bo_config$ntop_default
      }
      if (!bo_config$tune_lambdaw) {
        bo_fixed$lambdaW_grid <- bo_config$lambdaw_default
      }
      if (!bo_config$tune_lambdah) {
        bo_fixed$lambdaH_grid <- bo_config$lambdah_default
      }

      DeSurv::desurv_cv_bayesopt_refine(
        X = tar_data$ex,
        y = tar_data$sampInfo$time,
        d = tar_data$sampInfo$event,
        dataset = tar_data$sampInfo$dataset,
        samp_keeps = tar_data$samp_keeps,
        preprocess = TRUE,
        method_trans_train = bo_config$method_trans_train,
        engine = "warmstart",
        nfolds = bo_config$nfold,
        tol = bo_config$bo_tol,
        maxit = bo_config$bo_maxit,
        coarse_bounds = bounds,
        bo_fixed = bo_fixed,
        max_refinements = bo_config$bo_max_refinements,
        tol_gain = bo_config$bo_tol_gain,
        plateau = bo_config$bo_plateau,
        top_k = bo_config$bo_top_k,
        shrink_base = bo_config$bo_shrink_base,
        importance_gain = bo_config$bo_importance_gain,
        coarse_control = bo_config$bo_coarse_control,
        refine_control = bo_config$bo_refine_control,
        verbose = TRUE,
        parallel_grid = bo_config$desurv_parallel_grid,
        ncores_grid = bo_config$desurv_ncores_grid
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    ),
    cue = tar_cue(mode = "never")
  ),
  
  tar_target(
    tar_k_selection,
    select_bo_k_by_cv_se(desurv_bo_results)
  ),

  tar_target(
    tar_params_best,
    {
      params <- standardize_bo_params(desurv_bo_results$overall_best$params)
      if (!is.null(tar_k_selection$k_selected)) {
        params$k <- tar_k_selection$k_selected
      }
      params
    }
  ),
  
  tar_target(
    tar_ngene_value,
    {
      value <- tar_params_best$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ngene_default)
      } else {
        as.integer(round(value))
      }
    }
  ),
  tar_target(
    tar_ntop_value,
    {
      value <- tar_params_best$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ntop_default)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    tar_lambdaW_value,
    {
      value <- tar_params_best$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdaw_default)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    tar_lambdaH_value,
    {
      value <- tar_params_best$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdah_default)
      } else {
        as.numeric(value)
      }
    }
  ),
  
  tar_target(
    bo_results_dir,
    get_bo_results_dir(
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = bo_config$train_prefix,
      method_trans_train = bo_config$method_trans_train,
      bo_config_tag = bo_config$path_tag,
      bo_config_id = bo_config$config_id
    )
  ),
  
  tar_target(
    desurv_bo_history,
    {
      path <- file.path(bo_results_dir, "desurv_bo_history.csv")
      utils::write.csv(desurv_bo_results$history, path, row.names = FALSE)
      path
    },
    format = "file"
  ),
  
  tar_target(
    tar_data_filtered,
    preprocess_training_data(
      data = tar_data,
      ngene = tar_ngene_value,
      method_trans_train = bo_config$method_trans_train
    )
  ),

  tar_target(
    desurv_bo_results_alpha0,
    {
      bounds <- bo_config$desurv_bo_bounds
      bounds$alpha_grid <- NULL
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ngene_config, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ntop_config, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdaw_config, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdah_config, "lambdaH_grid", log_scale = TRUE)

      bo_fixed <- list(
        n_starts = bo_config$ninit,
        alpha_grid = 0
      )
      if (!bo_config$tune_ngene) {
        bo_fixed$ngene <- bo_config$ngene_default
      }
      if (!bo_config$tune_ntop) {
        bo_fixed$ntop <- bo_config$ntop_default
      }
      if (!bo_config$tune_lambdaw) {
        bo_fixed$lambdaW_grid <- bo_config$lambdaw_default
      }
      if (!bo_config$tune_lambdah) {
        bo_fixed$lambdaH_grid <- bo_config$lambdah_default
      }

      DeSurv::desurv_cv_bayesopt_refine(
        X = tar_data$ex,
        y = tar_data$sampInfo$time,
        d = tar_data$sampInfo$event,
        dataset = tar_data$sampInfo$dataset,
        samp_keeps = tar_data$samp_keeps,
        preprocess = TRUE,
        method_trans_train = bo_config$method_trans_train,
        engine = "warmstart",
        nfolds = bo_config$nfold,
        tol = bo_config$bo_tol,
        maxit = bo_config$bo_maxit,
        coarse_bounds = bounds,
        bo_fixed = bo_fixed,
        max_refinements = bo_config$bo_max_refinements,
        tol_gain = bo_config$bo_tol_gain,
        plateau = bo_config$bo_plateau,
        top_k = bo_config$bo_top_k,
        shrink_base = bo_config$bo_shrink_base,
        importance_gain = bo_config$bo_importance_gain,
        coarse_control = bo_config$bo_coarse_control,
        refine_control = bo_config$bo_refine_control,
        verbose = TRUE,
        parallel_grid = bo_config$desurv_parallel_grid,
        ncores_grid = bo_config$desurv_ncores_grid
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    tar_k_selection_alpha0,
    select_bo_k_by_cv_se(desurv_bo_results_alpha0)
  ),

  tar_target(
    tar_params_best_alpha0,
    {
      params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
      if (!is.null(tar_k_selection_alpha0$k_selected)) {
        params$k <- tar_k_selection_alpha0$k_selected
      }
      params$alpha <- 0
      params
    }
  ),

  tar_target(
    tar_ngene_value_alpha0,
    {
      value <- tar_params_best_alpha0$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ngene_default)
      } else {
        as.integer(round(value))
      }
    }
  ),
  tar_target(
    tar_ntop_value_alpha0,
    {
      value <- tar_params_best_alpha0$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ntop_default)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    tar_lambdaW_value_alpha0,
    {
      value <- tar_params_best_alpha0$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdaw_default)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    tar_lambdaH_value_alpha0,
    {
      value <- tar_params_best_alpha0$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdah_default)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    desurv_bo_history_alpha0,
    {
      path <- file.path(bo_results_dir, "desurv_bo_history_alpha0.csv")
      utils::write.csv(desurv_bo_results_alpha0$history, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  tar_target(
    tar_data_filtered_alpha0,
    {
      preprocess_training_data(
        data = tar_data,
        ngene = tar_ngene_value_alpha0,
        method_trans_train = bo_config$method_trans_train
      )
    }
    
  ),
  
  tar_target(
    desurv_bo_results_elbowk,
    {
      bounds <- bo_config$desurv_bo_bounds
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ngene_config, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$ntop_config, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdaw_config, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, bo_config$lambdah_config, "lambdaH_grid", log_scale = TRUE)
      bounds$k_grid=NULL
      bo_fixed <- listbo_fixed <- listbo_fixed <- list(n_starts = bo_config$ninit)
      if (!bo_config$tune_ngene) {
        bo_fixed$ngene <- bo_config$ngene_default
      }
      if (!bo_config$tune_ntop) {
        bo_fixed$ntop <- bo_config$ntop_default
      }
      if (!bo_config$tune_lambdaw) {
        bo_fixed$lambdaW_grid <- bo_config$lambdaw_default
      }
      if (!bo_config$tune_lambdah) {
        bo_fixed$lambdaH_grid <- bo_config$lambdah_default
      }
      
      bo_fixed$k_grid = std_nmf_selected_k

      DeSurv::desurv_cv_bayesopt_refine(
        X = tar_data$ex,
        y = tar_data$sampInfo$time,
        d = tar_data$sampInfo$event,
        dataset = tar_data$sampInfo$dataset,
        samp_keeps = tar_data$samp_keeps,
        preprocess = TRUE,
        method_trans_train = bo_config$method_trans_train,
        engine = "warmstart",
        nfolds = bo_config$nfold,
        tol = bo_config$bo_tol,
        maxit = bo_config$bo_maxit,
        coarse_bounds = bounds,
        bo_fixed = bo_fixed,
        max_refinements = bo_config$bo_max_refinements,
        tol_gain = bo_config$bo_tol_gain,
        plateau = bo_config$bo_plateau,
        top_k = bo_config$bo_top_k,
        shrink_base = bo_config$bo_shrink_base,
        importance_gain = bo_config$bo_importance_gain,
        coarse_control = bo_config$bo_coarse_control,
        refine_control = bo_config$bo_refine_control,
        verbose = TRUE,
        parallel_grid = bo_config$desurv_parallel_grid,
        ncores_grid = bo_config$desurv_ncores_grid
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    tar_params_best_elbowk,
    {
      params <- standardize_bo_params(desurv_bo_results_elbowk$overall_best$params)
      params
    }
  ),
  
  tar_target(
    tar_ngene_value_elbowk,
    {
      value <- tar_params_best_elbowk$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ngene_default)
      } else {
        as.integer(round(value))
      }
    }
  ),
  tar_target(
    tar_ntop_value_elbowk,
    {
      value <- tar_params_best_elbowk$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(bo_config$ntop_default)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    tar_lambdaW_value_elbowk,
    {
      value <- tar_params_best_elbowk$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdaw_default)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    tar_lambdaH_value_elbowk,
    {
      value <- tar_params_best_elbowk$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(bo_config$lambdah_default)
      } else {
        as.numeric(value)
      }
    }
  ),
  
  
  tar_target(
    tar_data_filtered_elbowk,
    {
      preprocess_training_data(
        data = tar_data,
        ngene = tar_ngene_value_elbowk,
        method_trans_train = bo_config$method_trans_train
      )
    }
    
  ),
  

  bo_bundle_target
)

run_bundle_target <- tar_target(
  run_bundle,
  list(
    label = run_config$label,
    config = run_config,
    bo_bundle = bo_bundle_selected,
    ntop_value = tar_ntop_value,
    ntop_value_alpha0 = tar_ntop_value_alpha0,
    training_results_dir = tar_training_results_dir,
    training_results_dir_alpha0 = tar_training_results_dir_alpha0,
    fit_desurv = tar_fit_desurv,
    fit_desurv_alpha0 = tar_fit_desurv_alpha0,
    tops_desurv = tar_tops_desurv,
    tops_desurv_alpha0 = tar_tops_desurv_alpha0
  )
)

run_bundles_target <- tar_target(
  run_bundles,
  list(run_bundle)
)

COMMON_DESURV_RUN_TARGETS <- list(
  tar_target(
    bo_bundle_selected,
    tar_bo_bundle
  ),
  tar_target(
    tar_training_results_dir,
    results_root_dir(
      ngene = bo_bundle_selected$ngene_value,
      tol = run_config$run_tol,
      maxit = run_config$run_maxit,
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = bo_bundle_selected$config$train_prefix,
      method_trans_train = bo_bundle_selected$config$method_trans_train,
      config_tag = run_config$path_tag
    )
  ),
  tar_target(
    tar_training_results_dir_alpha0,
    results_root_dir(
      ngene = bo_bundle_selected$ngene_value_alpha0,
      tol = run_config$run_tol,
      maxit = run_config$run_maxit,
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = paste0(bo_bundle_selected$config$train_prefix, "_alpha0"),
      method_trans_train = bo_bundle_selected$config$method_trans_train,
      config_tag = run_config$path_tag
    )
  ),
  tar_target(
    desurv_seed_fits,
    {
      seeds <- seq_len(run_config$ninit_full)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fit_i <- try(
          desurv_fit(
            X = bo_bundle_selected$data_filtered$ex,
            y = bo_bundle_selected$data_filtered$sampInfo$time,
            d = bo_bundle_selected$data_filtered$sampInfo$event,
            k = bo_bundle_selected$params_best$k,
            alpha = bo_bundle_selected$params_best$alpha,
            lambda = bo_bundle_selected$params_best$lambda,
            nu = bo_bundle_selected$params_best$nu,
            lambdaW = bo_bundle_selected$lambdaW_value,
            lambdaH = bo_bundle_selected$lambdaH_value,
            seed = seeds[i],
            tol = run_config$run_tol / 100,
            tol_init = run_config$run_tol,
            maxit = run_config$run_maxit,
            imaxit = run_config$run_maxit,
            ninit = 1,
            parallel_init = FALSE,
            verbose = FALSE
          ),
          silent = TRUE
        )
        if (!inherits(fit_i, "try-error") && inherits(fit_i, "desurv_fit")) {
          fits[[i]] <- fit_i
          scores[i] <- if (!is.null(fit_i$cindex)) fit_i$cindex else NA_real_
        }
      }
      keep <- !vapply(fits, is.null, logical(1))
      if (!any(keep)) {
        stop("No successful full-model fits were obtained.")
      }
      list(
        fits = fits[keep],
        seeds = seeds[keep],
        cindex = scores[keep]
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    ),
    cue = tar_cue(mode = "never")
  ),
  tar_target(
    desurv_consensus_init,
    {
      if (is.null(desurv_seed_fits$fits) || !length(desurv_seed_fits$fits)) {
        stop("Consensus initialization requires at least one successful seed fit.")
      }
      DeSurv::desurv_consensus_seed(
        fits = desurv_seed_fits$fits,
        X = bo_bundle_selected$data_filtered$ex,
        ntop = tar_ntop_value,
        k = bo_bundle_selected$params_best$k,
        min_frequency = 0.3 * run_config$ninit_full
      )
    }
  ),
  tar_target(
    tar_fit_desurv,
    {
      init_vals <- desurv_consensus_init
      desurv_fit(
        X = bo_bundle_selected$data_filtered$ex,
        y = bo_bundle_selected$data_filtered$sampInfo$time,
        d = bo_bundle_selected$data_filtered$sampInfo$event,
        k = bo_bundle_selected$params_best$k,
        alpha = bo_bundle_selected$params_best$alpha,
        lambda = bo_bundle_selected$params_best$lambda,
        nu = bo_bundle_selected$params_best$nu,
        lambdaW = bo_bundle_selected$lambdaW_value,
        lambdaH = bo_bundle_selected$lambdaH_value,
        W0 = init_vals$W0,
        H0 = init_vals$H0,
        beta0 = init_vals$beta0,
        seed = NULL,
        tol = run_config$run_tol / 100,
        maxit = run_config$run_maxit,
        verbose = FALSE
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),

  # --- CV-based cutpoint selection for DeSurv model ---
  tar_target(
    desurv_cv_cutpoint_result,
    run_cv_grid_point(
      data = bo_bundle_selected$data_filtered,
      k = bo_bundle_selected$params_best$k,
      alpha = bo_bundle_selected$params_best$alpha,
      fixed_params = list(
        lambda = bo_bundle_selected$params_best$lambda,
        nu = bo_bundle_selected$params_best$nu,
        lambdaW = bo_bundle_selected$lambdaW_value,
        lambdaH = bo_bundle_selected$lambdaH_value,
        ntop = tar_ntop_value
      ),
      nfolds = 5,
      n_starts = 30,
      seed = 123
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    desurv_cutpoint_eval,
    evaluate_cutpoint_zscores(
      desurv_cv_cutpoint_result,
      z_grid = seq(-2.0, 2.0, by = 0.2)
    )
  ),

  tar_target(
    desurv_cutpoint_summary,
    {
      desurv_cutpoint_eval |>
        dplyr::group_by(z_cutpoint) |>
        dplyr::summarise(
          mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
          se_cindex_dichot = sd(cindex_dichot, na.rm = TRUE) / sqrt(sum(!is.na(cindex_dichot))),
          mean_abs_logrank_z = mean(abs(logrank_z), na.rm = TRUE),
          se_abs_logrank_z = sd(abs(logrank_z), na.rm = TRUE) / sqrt(sum(!is.na(logrank_z))),
          .groups = "drop"
        )
    }
  ),

  tar_target(
    desurv_optimal_z_cutpoint,
    {
      desurv_cutpoint_summary |>
        dplyr::slice_max(mean_abs_logrank_z, n = 1, with_ties = FALSE) |>
        dplyr::pull(z_cutpoint)
    }
  ),

  tar_target(
    desurv_optimal_z_cutpoint_cindex,
    {
      desurv_cutpoint_summary |>
        dplyr::slice_max(mean_cindex_dichot, n = 1, with_ties = FALSE) |>
        dplyr::pull(z_cutpoint)
    }
  ),

  tar_target(
    fig_cutpoint_curves,
    {
      k_val <- bo_bundle_selected$params_best$k
      alpha_val <- bo_bundle_selected$params_best$alpha
      ntop_val <- tar_ntop_value
      ntop_na <- if (is.null(ntop_val) || is.na(ntop_val)) NA_real_ else as.numeric(ntop_val)

      out_dir <- "figures"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      p_cindex <- plot_cutpoint_curve_cindex(
        desurv_cutpoint_summary,
        k = k_val, alpha = alpha_val, ntop = ntop_na,
        optimal_z = desurv_optimal_z_cutpoint_cindex
      )
      fpath_cindex <- file.path(out_dir, sprintf("cutpoint_curve_cindex_%s.pdf", bo_label))
      ggplot2::ggsave(fpath_cindex, p_cindex, width = 5, height = 4)

      p_logrank <- plot_cutpoint_curve_logrank(
        desurv_cutpoint_summary,
        k = k_val, alpha = alpha_val, ntop = ntop_na,
        optimal_z = desurv_optimal_z_cutpoint
      )
      fpath_logrank <- file.path(out_dir, sprintf("cutpoint_curve_logrank_%s.pdf", bo_label))
      ggplot2::ggsave(fpath_logrank, p_logrank, width = 5, height = 4)

      c(fpath_cindex, fpath_logrank)
    },
    format = "file"
  ),

  # --- Per-factor cutpoint evaluation ---
  tar_target(
    desurv_cutpoint_eval_per_factor,
    evaluate_cutpoint_zscores_per_factor(
      desurv_cv_cutpoint_result,
      z_grid = seq(-2.0, 2.0, by = 0.2),
      k = bo_bundle_selected$params_best$k
    )
  ),

  tar_target(
    desurv_cutpoint_summary_per_factor,
    {
      desurv_cutpoint_eval_per_factor |>
        dplyr::group_by(factor_id, z_cutpoint) |>
        dplyr::summarise(
          mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
          se_cindex_dichot = sd(cindex_dichot, na.rm = TRUE) / sqrt(sum(!is.na(cindex_dichot))),
          mean_abs_logrank_z = mean(abs(logrank_z), na.rm = TRUE),
          se_abs_logrank_z = sd(abs(logrank_z), na.rm = TRUE) / sqrt(sum(!is.na(logrank_z))),
          .groups = "drop"
        )
    }
  ),

  tar_target(
    desurv_optimal_z_cutpoint_per_factor,
    {
      desurv_cutpoint_summary_per_factor |>
        dplyr::group_by(factor_id) |>
        dplyr::slice_max(mean_abs_logrank_z, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(factor_id, z_cutpoint)
    }
  ),

  tar_target(
    desurv_optimal_z_cutpoint_cindex_per_factor,
    {
      desurv_cutpoint_summary_per_factor |>
        dplyr::group_by(factor_id) |>
        dplyr::slice_max(mean_cindex_dichot, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(factor_id, z_cutpoint)
    }
  ),

  tar_target(
    desurv_factor_stats,
    {
      W <- tar_fit_desurv$W
      X <- bo_bundle_selected$data_filtered$ex
      Z <- t(X) %*% W
      k <- ncol(Z)
      logrank_cuts <- desurv_optimal_z_cutpoint_per_factor
      cindex_cuts <- desurv_optimal_z_cutpoint_cindex_per_factor
      lapply(seq_len(k), function(j) {
        list(
          factor_id = j,
          z_mean = mean(Z[, j]),
          z_sd = sd(Z[, j]),
          optimal_z_logrank = logrank_cuts$z_cutpoint[logrank_cuts$factor_id == j],
          optimal_z_cindex = cindex_cuts$z_cutpoint[cindex_cuts$factor_id == j]
        )
      })
    }
  ),

  tar_target(
    fig_cutpoint_curves_per_factor,
    {
      k_val <- bo_bundle_selected$params_best$k
      alpha_val <- bo_bundle_selected$params_best$alpha
      out_dir <- "figures"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(0)

      for (j in seq_len(k_val)) {
        fac_data <- desurv_cutpoint_summary_per_factor |>
          dplyr::filter(factor_id == j)
        if (nrow(fac_data) == 0) next

        opt_logrank <- desurv_optimal_z_cutpoint_per_factor$z_cutpoint[
          desurv_optimal_z_cutpoint_per_factor$factor_id == j]
        opt_cindex <- desurv_optimal_z_cutpoint_cindex_per_factor$z_cutpoint[
          desurv_optimal_z_cutpoint_cindex_per_factor$factor_id == j]

        # Logrank curve
        p_lr <- ggplot2::ggplot(
          fac_data, ggplot2::aes(x = z_cutpoint, y = mean_abs_logrank_z)
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 1.5) +
          ggplot2::geom_errorbar(
            ggplot2::aes(
              ymin = mean_abs_logrank_z - se_abs_logrank_z,
              ymax = mean_abs_logrank_z + se_abs_logrank_z
            ), width = 0.05
          ) +
          ggplot2::labs(
            x = "z-score cutpoint", y = "Mean |log-rank z|",
            title = sprintf("Factor %d |Log-rank z|", j)
          ) +
          ggplot2::theme_bw(base_size = 10)
        if (length(opt_logrank) && is.finite(opt_logrank)) {
          p_lr <- p_lr + ggplot2::geom_vline(
            xintercept = opt_logrank, linetype = "dashed", color = "red", linewidth = 0.5)
        }
        fpath_lr <- file.path(out_dir, sprintf("cutpoint_curve_factor%d_logrank_%s.pdf", j, bo_label))
        ggplot2::ggsave(fpath_lr, p_lr, width = 5, height = 4)
        paths <- c(paths, fpath_lr)

        # Cindex curve
        p_ci <- ggplot2::ggplot(
          fac_data, ggplot2::aes(x = z_cutpoint, y = mean_cindex_dichot)
        ) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 1.5) +
          ggplot2::geom_errorbar(
            ggplot2::aes(
              ymin = mean_cindex_dichot - se_cindex_dichot,
              ymax = mean_cindex_dichot + se_cindex_dichot
            ), width = 0.05
          ) +
          ggplot2::labs(
            x = "z-score cutpoint", y = "Mean C-index (dichotomized)",
            title = sprintf("Factor %d Dichot. C-index", j)
          ) +
          ggplot2::theme_bw(base_size = 10)
        if (length(opt_cindex) && is.finite(opt_cindex)) {
          p_ci <- p_ci + ggplot2::geom_vline(
            xintercept = opt_cindex, linetype = "dashed", color = "red", linewidth = 0.5)
        }
        fpath_ci <- file.path(out_dir, sprintf("cutpoint_curve_factor%d_cindex_%s.pdf", j, bo_label))
        ggplot2::ggsave(fpath_ci, p_ci, width = 5, height = 4)
        paths <- c(paths, fpath_ci)
      }

      paths
    },
    format = "file"
  ),

  tar_target(
    desurv_lp_stats,
    {
      lp <- compute_lp(
        tar_fit_desurv$W,
        tar_fit_desurv$beta,
        bo_bundle_selected$data_filtered$ex,
        tar_ntop_value
      )
      lp_mean <- mean(lp, na.rm = TRUE)
      lp_sd <- sd(lp, na.rm = TRUE)
      list(
        lp_mean = lp_mean,
        lp_sd = lp_sd,
        optimal_z_cutpoint = desurv_optimal_z_cutpoint,
        cutpoint_abs = desurv_optimal_z_cutpoint * lp_sd + lp_mean,
        optimal_z_cutpoint_cindex = desurv_optimal_z_cutpoint_cindex,
        cutpoint_abs_cindex = desurv_optimal_z_cutpoint_cindex * lp_sd + lp_mean
      )
    }
  ),

  tar_target(
    desurv_seed_fits_alpha0,
    {
      seeds <- seq_len(run_config$ninit_full)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fit_i <- try(
          desurv_fit(
            X = bo_bundle_selected$data_filtered_alpha0$ex,
            y = bo_bundle_selected$data_filtered_alpha0$sampInfo$time,
            d = bo_bundle_selected$data_filtered_alpha0$sampInfo$event,
            k = bo_bundle_selected$params_best_alpha0$k,
            alpha = bo_bundle_selected$params_best_alpha0$alpha,
            lambda = bo_bundle_selected$params_best_alpha0$lambda,
            nu = bo_bundle_selected$params_best_alpha0$nu,
            lambdaW = bo_bundle_selected$lambdaW_value_alpha0,
            lambdaH = bo_bundle_selected$lambdaH_value_alpha0,
            seed = seeds[i],
            tol = run_config$run_tol / 100,
            tol_init = run_config$run_tol,
            maxit = run_config$run_maxit,
            imaxit = run_config$run_maxit,
            ninit = 1,
            parallel_init = FALSE,
            verbose = FALSE
          ),
          silent = TRUE
        )
        if (!inherits(fit_i, "try-error") && inherits(fit_i, "desurv_fit")) {
          fits[[i]] <- fit_i
          scores[i] <- if (!is.null(fit_i$cindex)) fit_i$cindex else NA_real_
        }
      }
      keep <- !vapply(fits, is.null, logical(1))
      if (!any(keep)) {
        stop("No successful full-model fits were obtained for alpha=0.")
      }
      list(
        fits = fits[keep],
        seeds = seeds[keep],
        cindex = scores[keep]
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    desurv_consensus_init_alpha0,
    {
      if (is.null(desurv_seed_fits_alpha0$fits) || !length(desurv_seed_fits_alpha0$fits)) {
        stop("Consensus initialization requires at least one successful alpha=0 seed fit.")
      }
      DeSurv::desurv_consensus_seed(
        fits = desurv_seed_fits_alpha0$fits,
        X = bo_bundle_selected$data_filtered_alpha0$ex,
        ntop = tar_ntop_value_alpha0,
        k = bo_bundle_selected$params_best_alpha0$k,
        min_frequency = 0.3 * run_config$ninit_full
      )
    }
  ),

  tar_target(
    tar_fit_desurv_alpha0,
    {
      init_vals <- desurv_consensus_init_alpha0
      desurv_fit(
        X = bo_bundle_selected$data_filtered_alpha0$ex,
        y = bo_bundle_selected$data_filtered_alpha0$sampInfo$time,
        d = bo_bundle_selected$data_filtered_alpha0$sampInfo$event,
        k = bo_bundle_selected$params_best_alpha0$k,
        alpha = bo_bundle_selected$params_best_alpha0$alpha,
        lambda = bo_bundle_selected$params_best_alpha0$lambda,
        nu = bo_bundle_selected$params_best_alpha0$nu,
        lambdaW = bo_bundle_selected$lambdaW_value_alpha0,
        lambdaH = bo_bundle_selected$lambdaH_value_alpha0,
        W0 = init_vals$W0,
        H0 = init_vals$H0,
        beta0 = init_vals$beta0,
        seed = NULL,
        tol = run_config$run_tol / 100,
        maxit = run_config$run_maxit,
        verbose = FALSE
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),
  
  tar_target(
    desurv_seed_fits_elbowk,
    {
      seeds <- seq_len(run_config$ninit_full)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fit_i <- try(
          desurv_fit(
            X = tar_data_filtered_elbowk$ex,
            y = tar_data_filtered_elbowk$sampInfo$time,
            d = tar_data_filtered_elbowk$sampInfo$event,
            k = std_nmf_selected_k,
            alpha = tar_params_best_elbowk$alpha,
            lambda = tar_params_best_elbowk$lambda,
            nu = tar_params_best_elbowk$nu,
            lambdaW = tar_lambdaW_value_elbowk,
            lambdaH = tar_lambdaH_value_elbowk,
            seed = seeds[i],
            tol = run_config$run_tol / 100,
            tol_init = run_config$run_tol,
            maxit = run_config$run_maxit,
            imaxit = run_config$run_maxit,
            ninit = 1,
            parallel_init = FALSE,
            verbose = FALSE
          ),
          silent = TRUE
        )
        if (!inherits(fit_i, "try-error") && inherits(fit_i, "desurv_fit")) {
          fits[[i]] <- fit_i
          scores[i] <- if (!is.null(fit_i$cindex)) fit_i$cindex else NA_real_
        }
      }
      keep <- !vapply(fits, is.null, logical(1))
      if (!any(keep)) {
        stop("No successful full-model fits were obtained.")
      }
      list(
        fits = fits[keep],
        seeds = seeds[keep],
        cindex = scores[keep]
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    ),
    cue = tar_cue(mode = "never")
  ),
  tar_target(
    desurv_consensus_init_elbowk,
    {
      if (is.null(desurv_seed_fits_elbowk$fits) || !length(desurv_seed_fits_elbowk$fits)) {
        stop("Consensus initialization requires at least one successful seed fit.")
      }
      DeSurv::desurv_consensus_seed(
        fits = desurv_seed_fits_elbowk$fits,
        X = tar_data_filtered_elbowk$ex,
        ntop = tar_ntop_value_elbowk,
        k = tar_params_best_elbowk$k,
        min_frequency = 0.3 * run_config$ninit_full
      )
    }
  ),
  tar_target(
    tar_fit_desurv_elbowk,
    {
      init_vals <- desurv_consensus_init_elbowk
      desurv_fit(
        X = tar_data_filtered_elbowk$ex,
        y = tar_data_filtered_elbowk$sampInfo$time,
        d = tar_data_filtered_elbowk$sampInfo$event,
        k = std_nmf_selected_k,
        alpha = tar_params_best_elbowk$alpha,
        lambda = tar_params_best_elbowk$lambda,
        nu = tar_params_best_elbowk$nu,
        lambdaW = tar_lambdaW_value_elbowk,
        lambdaH = tar_lambdaH_value_elbowk,
        W0 = init_vals$W0,
        H0 = init_vals$H0,
        beta0 = init_vals$beta0,
        seed = NULL,
        tol = run_config$run_tol / 100,
        maxit = run_config$run_maxit,
        verbose = FALSE
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),
  

  tar_target(
    fit_std,
    nmf(
      bo_bundle_selected$data_filtered$ex,
      run_config$std_nmf_k_grid,
      nrun = bo_bundle_selected$config$ninit,
      method = "lee",
      .options = paste0("p", bo_bundle_selected$config$ninit)
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    ),
    cue = tar_cue(mode = "never")
  ),
  
  tar_target(
    std_nmf_k_selection_plots,
    {
      save_dir=file.path(tar_training_results_dir,"std_nmf_k_selection")
      dir.create(save_dir,showWarnings = FALSE)
      path = file.path(save_dir,paste0("std_nmf_k_selection_plots.png"))
      if (file.exists(path)) {
        # 1. make a timestamped backup of the old file
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        backup_path <- file.path(dirname(path),
                                 paste0("std_nmf_k_selection_plots_", timestamp, ".png")
        )
        file.copy(from = path, to = backup_path)
      }
      png(filename=path)
      print(NMF::plot(fit_std))
      dev.off()
      path
    },
    format="file"
  ),
  
  tar_target(
    std_nmf_selected_k,
    {
      #select k
      ranks = fit_std$measures$rank
      resids = fit_std$measures$residuals
      pick_k_elbow(ranks,resids)
    }
  ),
  
 

  tar_target(
    fit_std_elbowk,
    {
      # selected k model
      selected_k = std_nmf_selected_k
      fit_nmf = fit_std$fit[[as.character(selected_k)]]
      
      # WTX
      W = fit_nmf@fit@W
      X = tar_data_filtered_elbowk$ex
      keep_genes = intersect(rownames(W),rownames(X))
      W = W[keep_genes,]
      X = X[keep_genes,]
      XtW = t(X)%*%W
      colnames(XtW) = paste0("XtW",1:ncol(XtW))
      
      # dataframe
      df = cbind(tar_data_filtered_elbowk$sampInfo,XtW)
      
      beta = fit_cox_model(XtW,df,bo_config$nfold)
      
      fit = list()
      fit$W = fit_nmf@fit@W
      fit$H = fit_nmf@fit@H
      fit$beta = beta
      
      fit
    }
  ),
  
  
  tar_target(
    fit_std_desurvk,
    {
      # selected k model
      
      selected_k = bo_bundle_selected$params_best$k
      fit_nmf = fit_std$fit[[as.character(selected_k)]]
      
      # WTX
      W = fit_nmf@fit@W
      X = tar_data_filtered$ex
      keep_genes = intersect(rownames(W),rownames(X))
      W = W[keep_genes,]
      X = X[keep_genes,]
      XtW = t(X)%*%W
      colnames(XtW) = paste0("XtW",1:ncol(XtW))
      
      # dataframe
      df = cbind(tar_data_filtered$sampInfo,XtW)
      
      beta = fit_cox_model(XtW,df,bo_config$nfold)
      
      fit = list()
      fit$W = fit_nmf@fit@W
      fit$H = fit_nmf@fit@H
      fit$beta = beta
      
      fit
    }
  ),
  
  
  
  run_bundle_target,
  run_bundles_target,
  tar_target(
    tar_tops_desurv,
    get_top_genes(W = tar_fit_desurv$W, ntop = tar_ntop_value)
  ),
  tar_target(
    gene_overlap_desurv,
    {
      create_table(tops = tar_tops_desurv$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),

  tar_target(
    tar_tops_desurv_alpha0,
    get_top_genes(W = tar_fit_desurv_alpha0$W, ntop = tar_ntop_value_alpha0)
  ),
  tar_target(
    gene_overlap_desurv_alpha0,
    {
      create_table(tops = tar_tops_desurv_alpha0$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),
  tar_target(
    tar_tops_std_elbowk,
    {
      get_top_genes(W = fit_std_elbowk$W, ntop = tar_ntop_value_elbowk)
    }
  ),
  tar_target(
    gene_overlap_std_elbowk,
    {
      create_table(tops = tar_tops_std_elbowk$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),
  
  tar_target(
    tar_tops_std_desurvk,
    {
      get_top_genes(W = fit_std_desurvk$W, ntop = tar_ntop_value)
    }
  ),
  tar_target(
    gene_overlap_std_desurvk,
    {
      create_table(tops = tar_tops_std_desurvk$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),
  
  tar_target(
    tar_tops_desurv_elbowk,
    {
      get_top_genes(W = tar_fit_desurv_elbowk$W, ntop = bo_config$ntop_default)
    }
  ),
  tar_target(
    gene_overlap_desurv_elbowk,
    {
      create_table(tops = tar_tops_desurv_elbowk$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),
  # tar_target(
  #   tops_std,
  #   {
  #     W = fit_std_selected_k@fit@W
  #     get_top_genes(W = W, ntop = ntop_value)
  #   }
  #   
  # ),
  # tar_target(
  #   gene_overlap_std,
  #   {
  #     create_table(tops = tops_std$top_genes, gene_lists = top_genes,
  #                  which.lists = "DECODER", color.lists = colors)
  #   }
  # ),
  # 
  # 
  ### characterization
  tar_target(
    ora_analysis_desurv,
    {
      universe = rownames(bo_bundle_selected$data_filtered$ex)
      organism <- org.Hs.eg.db
      ora(tar_tops_desurv$top_genes,universe,organism)
    }
  ),

  tar_target(
    ora_analysis_desurv_alpha0,
    {
      universe = rownames(bo_bundle_selected$data_filtered_alpha0$ex)
      organism <- org.Hs.eg.db
      ora(tar_tops_desurv_alpha0$top_genes,universe,organism)
    }
  ),
  
  tar_target(
    ora_analysis_desurv_elbowk,
    {
      universe = rownames(tar_data_filtered_elbowk$ex)
      organism <- org.Hs.eg.db
      ora(tar_tops_desurv_elbowk$top_genes,universe,organism)
    }
  ),
  
  tar_target(
    ora_analysis_std_elbowk,
    {
      universe = rownames(tar_data_filtered$ex)
      organism <- org.Hs.eg.db
      ora(tar_tops_std_elbowk$top_genes,universe,organism)
    }
  ),
  
  tar_target(
    ora_analysis_std_desurvk,
    {
      universe = rownames(tar_data_filtered$ex)
      organism <- org.Hs.eg.db
      ora(tar_tops_std_desurvk$top_genes,universe,organism)
    }
  )
)

COMMON_DESURV_VAL_TARGETS <- list(
  tar_target(
    val_run_bundle,
    run_bundle
  ),
  tar_target(
    data_val_filtered,
    {
      datasets_named <- ensure_named_datasets(data_val)
      setNames(
        lapply(
          seq_along(datasets_named),
          function(idx) {
            dataset <- datasets_named[[idx]]
            dataname <- names(datasets_named)[idx]
            genes_train <- rownames(val_run_bundle$bo_bundle$data_filtered$ex)
            preprocess_validation_data(
              dataset = dataset,
              genes = genes_train,
              ngene = val_run_bundle$bo_bundle$ngene_value,
              method_trans_train = val_run_bundle$bo_bundle$config$method_trans_train,
              dataname = dataname,
              transform_target = val_run_bundle$bo_bundle$data_filtered$transform_target,
              zero_fill_missing = TRUE
            )
          }
        ),
        names(datasets_named)
      )
    },
    iteration="list"
  ),

  # Merged PACA_AU for survival validation (deduplicates overlapping subjects)
  # Reconstruct names lost by iteration="list" aggregation, then merge
  tar_target(
    data_val_filtered_surv,
    {
      val_named <- data_val_filtered
      # Restore dataset names from each element's dataname field
      nms <- vapply(val_named, function(x) {
        if (!is.null(x$dataname) && nzchar(x$dataname)) x$dataname
        else infer_dataset_name(x, fallback = "unknown")
      }, character(1))
      names(val_named) <- nms
      merge_paca_au_datasets(val_named)
    }
  ),

  tar_target(
    data_val_filtered_elbowk,
    {
      datasets_named <- ensure_named_datasets(data_val)
      setNames(
        lapply(
          seq_along(datasets_named),
          function(idx) {
            dataset <- datasets_named[[idx]]
            dataname <- names(datasets_named)[idx]
            genes_train <- rownames(tar_data_filtered_elbowk$ex)
            preprocess_validation_data(
              dataset = dataset,
              genes = genes_train,
              ngene = tar_ngene_value_elbowk,
              method_trans_train = val_run_bundle$bo_bundle$config$method_trans_train,
              dataname = dataname,
              transform_target = val_run_bundle$bo_bundle$data_filtered$transform_target,
              zero_fill_missing = TRUE
            )
          }
        ),
        names(datasets_named)
      )
    },
    iteration = "list"
  ),
  
  tar_target(
    val_predictions_desurv,
    desurv_predict_validation(
      fit = val_run_bundle$fit_desurv,
      data_list = data_val_filtered,
      top_genes = val_run_bundle$tops_desurv$top_genes
    )
  ),
  tar_target(
    val_predictions_desurv_alpha0,
    desurv_predict_validation(
      fit = val_run_bundle$fit_desurv_alpha0,
      data_list = data_val_filtered,
      top_genes = val_run_bundle$tops_desurv_alpha0$top_genes
    )
  ),
  tar_target(
    val_predictions_desurv_elbowk,
    desurv_predict_validation(
      fit = tar_fit_desurv_elbowk,
      data_list = data_val_filtered_elbowk,
      top_genes = tar_tops_desurv_elbowk$top_genes
    )
  ),
  tar_target(
    val_predictions_std_elbowk,
    desurv_predict_validation(
      fit = fit_std_elbowk,
      data_list = data_val_filtered,
      top_genes = tar_tops_std_elbowk$top_genes
    )
  ),
  tar_target(
    val_predictions_std_desurvk,
    desurv_predict_validation(
      fit = fit_std_desurvk,
      data_list = data_val_filtered,
      top_genes = tar_tops_std_desurvk$top_genes
    )
  ),
  
  tar_target(
    val_latent_desurv,
    {
      latent <- desurv_collect_validation_latent(
        fit = val_run_bundle$fit_desurv,
        data_list = data_val_filtered,
        top_genes = val_run_bundle$tops_desurv$top_genes
      )
      write_validation_latent_outputs(
        latent_list = latent,
        base_dir = file.path(
          val_run_bundle$training_results_dir,
          "validation",
          val_config_effective$config_id,
          "desurv"
        )
      )
      latent
    },
    iteration = "list"
  ),
  tar_target(
    val_latent_desurv_alpha0,
    {
      latent <- desurv_collect_validation_latent(
        fit = val_run_bundle$fit_desurv_alpha0,
        data_list = data_val_filtered,
        top_genes = val_run_bundle$tops_desurv_alpha0$top_genes
      )
      write_validation_latent_outputs(
        latent_list = latent,
        base_dir = file.path(
          val_run_bundle$training_results_dir_alpha0,
          "validation",
          val_config_effective$config_id,
          "desurv_alpha0"
        )
      )
      latent
    }
  ),
  tar_target(
    val_latent_desurv_elbowk,
    {
      latent <- desurv_collect_validation_latent(
        fit = tar_fit_desurv_elbowk,
        data_list = data_val_filtered_elbowk,
        top_genes = tar_tops_desurv_elbowk$top_genes
      )
      latent
    },
    iteration = "list"
  ),
  
  tar_target(
    val_latent_std_elbowk,
    {
      latent <- desurv_collect_validation_latent(
        fit = fit_std_elbowk,
        data_list = data_val_filtered_elbowk,
        top_genes = tar_tops_std_elbowk$top_genes
      )
      latent
    },
    iteration = "list"
  ),
  
  tar_target(
    val_latent_std_desurvk,
    {
      latent <- desurv_collect_validation_latent(
        fit = fit_std_desurvk,
        data_list = data_val_filtered,
        top_genes = tar_tops_std_desurvk$top_genes
      )
      latent
    },
    iteration = "list"
  ),
  
  
  tar_target(
    val_cindex_desurv,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_desurv)
      if (nrow(summary_tbl)) {
        dir.create(
          file.path(val_run_bundle$training_results_dir, "validation", val_config_effective$config_id),
          recursive = TRUE,
          showWarnings = FALSE
        )
        utils::write.csv(
          summary_tbl,
          file = file.path(
            val_run_bundle$training_results_dir,
            "validation",
            val_config_effective$config_id,
            "val_cindex_desurv.csv"
          ),
          row.names = FALSE
        )
      }
      summary_tbl
    }
  ),

  # --- Dichotomized risk group evaluation on validation datasets ---
  tar_target(
    val_dichot_desurv,
    {
      W <- val_run_bundle$fit_desurv$W
      beta <- val_run_bundle$fit_desurv$beta
      ntop <- val_run_bundle$ntop_value
      train_lp_mean <- desurv_lp_stats$lp_mean
      train_lp_sd <- desurv_lp_stats$lp_sd

      cutpoints <- list(
        logrank = desurv_optimal_z_cutpoint,
        cindex = desurv_optimal_z_cutpoint_cindex
      )

      ds_names <- names(data_val_filtered_surv)
      all_rows <- list()

      for (method in names(cutpoints)) {
        z_cut <- cutpoints[[method]]
        rows <- lapply(ds_names, function(ds_name) {
          val_ds <- data_val_filtered_surv[[ds_name]]
          val_genes <- rownames(val_ds$ex)
          common_genes <- intersect(rownames(W), val_genes)
          if (length(common_genes) < 2) {
            return(tibble::tibble(
              dataset = ds_name, n_samples = 0L, n_events = 0L,
              cindex_dichot = NA_real_, logrank_z = NA_real_,
              z_cutpoint = z_cut, method = method
            ))
          }

          W_common <- W[common_genes, , drop = FALSE]
          X_val <- val_ds$ex[common_genes, , drop = FALSE]
          lp <- compute_lp(W_common, beta, X_val, ntop)

          time_val <- val_ds$sampInfo$time
          event_val <- as.integer(val_ds$sampInfo$event)
          valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
          if (length(valid_idx) < 2) {
            return(tibble::tibble(
              dataset = ds_name, n_samples = length(valid_idx),
              n_events = sum(event_val[valid_idx]),
              cindex_dichot = NA_real_, logrank_z = NA_real_,
              z_cutpoint = z_cut, method = method
            ))
          }

          lp <- lp[valid_idx]
          time_val <- time_val[valid_idx]
          event_val <- event_val[valid_idx]

          if (!is.finite(train_lp_sd) || train_lp_sd <= 0) {
            return(tibble::tibble(
              dataset = ds_name, n_samples = length(valid_idx),
              n_events = sum(event_val),
              cindex_dichot = NA_real_, logrank_z = NA_real_,
              z_cutpoint = z_cut, method = method
            ))
          }

          z_val <- (lp - train_lp_mean) / train_lp_sd
          group <- as.integer(z_val > z_cut)

          ci <- NA_real_
          lr_z <- NA_real_
          if (length(unique(group)) == 2 && sum(event_val) > 0) {
            ci <- tryCatch({
              cc <- survival::concordance(
                survival::Surv(time_val, event_val) ~ group,
                reverse = TRUE
              )
              cc$concordance
            }, error = function(e) NA_real_)
            ds_col <- val_ds$sampInfo$dataset[valid_idx]
            has_strata <- !is.null(ds_col) && length(unique(ds_col)) > 1
            lr_z <- compute_logrank_z(
              time_val, event_val, group,
              strata = if (has_strata) ds_col else NULL
            )
          }

          tibble::tibble(
            dataset = ds_name,
            n_samples = length(valid_idx),
            n_events = sum(event_val),
            cindex_dichot = ci,
            logrank_z = lr_z,
            z_cutpoint = z_cut,
            method = method
          )
        })
        all_rows <- c(all_rows, rows)
      }

      dplyr::bind_rows(all_rows)
    }
  ),

  # --- Validation KM plots for dichotomized risk groups ---
  tar_target(
    fig_val_km_dichot,
    {
      W <- val_run_bundle$fit_desurv$W
      beta <- val_run_bundle$fit_desurv$beta
      ntop <- val_run_bundle$ntop_value
      lp_mean <- desurv_lp_stats$lp_mean
      lp_sd <- desurv_lp_stats$lp_sd

      cutpoints <- list(
        logrank = desurv_optimal_z_cutpoint,
        cindex = desurv_optimal_z_cutpoint_cindex
      )

      out_dir <- file.path("figures", "km_dichot")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(0)

      for (method in names(cutpoints)) {
        z_cut <- cutpoints[[method]]

        fit_entry <- list(
          fit = list(W = W, beta = beta),
          ntop = ntop,
          z_cutpoint = z_cut,
          lp_mean = lp_mean,
          lp_sd = lp_sd,
          k = val_run_bundle$bo_bundle$params_best$k,
          alpha = val_run_bundle$bo_bundle$params_best$alpha %||% NA_real_
        )

        # Per-dataset KM plots
        for (ds_name in names(data_val_filtered_surv)) {
          p <- plot_km_validation(fit_entry, data_val_filtered_surv[[ds_name]],
                                   ds_name)
          if (is.null(p)) next
          fpath <- file.path(out_dir, sprintf("km_val_%s_%s_%s.pdf",
                                               ds_name, method, bo_label))
          save_ggsurvplot(p, fpath, width = 7, height = 5)
          paths <- c(paths, fpath)
        }

        # Pooled KM plot
        p_pooled <- plot_km_validation_pooled(fit_entry, data_val_filtered_surv)
        if (!is.null(p_pooled)) {
          fpath <- file.path(out_dir, sprintf("km_val_pooled_%s_%s.pdf",
                                               method, bo_label))
          save_ggsurvplot(p_pooled, fpath, width = 7, height = 5)
          paths <- c(paths, fpath)
        }
      }

      paths
    },
    format = "file"
  ),

  # --- Subtype overlap with dichotomized risk groups ---
  tar_target(
    fig_subtype_overlap,
    {
      W <- val_run_bundle$fit_desurv$W
      beta <- val_run_bundle$fit_desurv$beta
      ntop <- val_run_bundle$ntop_value
      lp_mean <- desurv_lp_stats$lp_mean
      lp_sd <- desurv_lp_stats$lp_sd

      cutpoints <- list(
        logrank = desurv_optimal_z_cutpoint,
        cindex = desurv_optimal_z_cutpoint_cindex
      )

      out_dir <- file.path("figures", "km_dichot")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(0)

      for (method in names(cutpoints)) {
        z_cut <- cutpoints[[method]]

        # --- Training data ---
        train_ids <- colnames(bo_bundle_selected$data_filtered$ex)
        train_si <- tar_data$sampInfo
        # Match by rownames (set during preprocessing) or ID column
        train_match <- match(train_ids, rownames(train_si))
        if (any(is.na(train_match)) && "ID" %in% names(train_si)) {
          train_match <- match(train_ids, train_si$ID)
        }
        if (!all(is.na(train_match))) {
          train_si_matched <- train_si[train_match[!is.na(train_match)], ]
          train_X <- bo_bundle_selected$data_filtered$ex[, !is.na(train_match), drop = FALSE]
          if ("PurIST" %in% names(train_si_matched) && "DeCAF" %in% names(train_si_matched)) {
            train_group <- compute_risk_group(W, beta, train_X, ntop, lp_mean, lp_sd, z_cut)
            train_df <- data.frame(
              group = train_group,
              PurIST = train_si_matched$PurIST,
              DeCAF = train_si_matched$DeCAF,
              stringsAsFactors = FALSE
            )
            train_df <- train_df[complete.cases(train_df), ]
            if (nrow(train_df) > 0) {
              p <- plot_subtype_overlap(train_df,
                dataset_label = sprintf("Training (%s cutpoint)", method))
              fpath <- file.path(out_dir, sprintf("subtype_overlap_train_%s_%s.pdf",
                                                   method, bo_label))
              ggplot2::ggsave(fpath, p, width = 8, height = 5)
              paths <- c(paths, fpath)

              # Enrichment plot (confounding-aware)
              pe <- plot_subtype_enrichment(train_df,
                dataset_label = sprintf("Training (%s cutpoint)", method))
              fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_train_%s_%s.pdf",
                                                     method, bo_label))
              ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
              paths <- c(paths, fpath_e)
            }
          }
        }

        # --- Per validation dataset ---
        pooled_rows <- list()
        for (ds_name in names(data_val_filtered_surv)) {
          val_ds <- data_val_filtered_surv[[ds_name]]
          val_genes <- rownames(val_ds$ex)
          common_genes <- intersect(rownames(W), val_genes)
          if (length(common_genes) < 2) next

          W_common <- W[common_genes, , drop = FALSE]
          X_val <- val_ds$ex[common_genes, , drop = FALSE]

          time_val <- val_ds$sampInfo$time
          event_val <- val_ds$sampInfo$event
          valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
          if (length(valid_idx) < 2) next

          val_group <- compute_risk_group(W_common, beta,
            X_val[, valid_idx, drop = FALSE], ntop, lp_mean, lp_sd, z_cut)
          val_si <- val_ds$sampInfo[valid_idx, ]

          if ("PurIST" %in% names(val_si) && "DeCAF" %in% names(val_si)) {
            ds_df <- data.frame(
              group = val_group,
              PurIST = val_si$PurIST,
              DeCAF = val_si$DeCAF,
              stringsAsFactors = FALSE
            )
            ds_df <- ds_df[complete.cases(ds_df), ]
            pooled_rows[[ds_name]] <- ds_df
            if (nrow(ds_df) > 0) {
              p <- plot_subtype_overlap(ds_df,
                dataset_label = sprintf("%s (%s cutpoint)", ds_name, method))
              fpath <- file.path(out_dir, sprintf("subtype_overlap_%s_%s_%s.pdf",
                                                   ds_name, method, bo_label))
              ggplot2::ggsave(fpath, p, width = 8, height = 5)
              paths <- c(paths, fpath)

              # Enrichment plot (confounding-aware)
              pe <- plot_subtype_enrichment(ds_df,
                dataset_label = sprintf("%s (%s cutpoint)", ds_name, method))
              fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_%s_%s_%s.pdf",
                                                     ds_name, method, bo_label))
              ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
              paths <- c(paths, fpath_e)
            }
          }
        }

        # --- Pooled validation ---
        if (length(pooled_rows) > 0) {
          pooled_df <- do.call(rbind, pooled_rows)
          if (nrow(pooled_df) > 0) {
            p <- plot_subtype_overlap(pooled_df,
              dataset_label = sprintf("Pooled validation (%s cutpoint)", method))
            fpath <- file.path(out_dir, sprintf("subtype_overlap_pooled_%s_%s.pdf",
                                                 method, bo_label))
            ggplot2::ggsave(fpath, p, width = 8, height = 5)
            paths <- c(paths, fpath)

            # Enrichment plot (confounding-aware)
            pe <- plot_subtype_enrichment(pooled_df,
              dataset_label = sprintf("Pooled validation (%s cutpoint)", method))
            fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_pooled_%s_%s.pdf",
                                                   method, bo_label))
            ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
            paths <- c(paths, fpath_e)
          }
        }
      }

      paths
    },
    format = "file"
  ),

  # --- Per-factor KM figures (validation + training) ---
  tar_target(
    fig_val_km_dichot_per_factor,
    {
      fit <- val_run_bundle$fit_desurv
      k_val <- val_run_bundle$bo_bundle$params_best$k
      alpha_val <- val_run_bundle$bo_bundle$params_best$alpha %||% NA_real_
      factor_stats <- desurv_factor_stats

      out_dir <- file.path("figures", "km_dichot")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(0)

      for (fstat in factor_stats) {
        j <- fstat$factor_id
        z_mean <- fstat$z_mean
        z_sd <- fstat$z_sd

        cutpoints <- list(
          logrank = fstat$optimal_z_logrank,
          cindex = fstat$optimal_z_cindex
        )

        for (method in names(cutpoints)) {
          z_cut <- cutpoints[[method]]
          if (!length(z_cut) || is.na(z_cut)) next

          # Training KM
          p_train <- plot_km_training_factor(
            fit, bo_bundle_selected$data_filtered,
            j, z_mean, z_sd, z_cut, k_val, alpha_val
          )
          if (!is.null(p_train)) {
            fpath <- file.path(out_dir, sprintf("km_val_train_factor%d_%s_%s.pdf",
                                                 j, method, bo_label))
            save_ggsurvplot(p_train, fpath, width = 7, height = 5)
            paths <- c(paths, fpath)
          }

          # Per-dataset validation KM
          for (ds_name in names(data_val_filtered_surv)) {
            p <- plot_km_validation_factor(
              fit, data_val_filtered_surv[[ds_name]], ds_name,
              j, z_mean, z_sd, z_cut, k_val, alpha_val
            )
            if (is.null(p)) next
            fpath <- file.path(out_dir, sprintf("km_val_%s_factor%d_%s_%s.pdf",
                                                 ds_name, j, method, bo_label))
            save_ggsurvplot(p, fpath, width = 7, height = 5)
            paths <- c(paths, fpath)
          }

          # Pooled validation KM
          p_pooled <- plot_km_validation_pooled_factor(
            fit, data_val_filtered_surv,
            j, z_mean, z_sd, z_cut, k_val, alpha_val
          )
          if (!is.null(p_pooled)) {
            fpath <- file.path(out_dir, sprintf("km_val_pooled_factor%d_%s_%s.pdf",
                                                 j, method, bo_label))
            save_ggsurvplot(p_pooled, fpath, width = 7, height = 5)
            paths <- c(paths, fpath)
          }
        }
      }

      paths
    },
    format = "file"
  ),

  # --- Per-factor subtype overlap figures ---
  tar_target(
    fig_subtype_overlap_per_factor,
    {
      W <- val_run_bundle$fit_desurv$W
      factor_stats <- desurv_factor_stats

      out_dir <- file.path("figures", "km_dichot")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      paths <- character(0)

      for (fstat in factor_stats) {
        j <- fstat$factor_id
        z_mean <- fstat$z_mean
        z_sd <- fstat$z_sd
        if (!is.finite(z_sd) || z_sd <= 0) next

        cutpoints <- list(
          logrank = fstat$optimal_z_logrank,
          cindex = fstat$optimal_z_cindex
        )

        for (method in names(cutpoints)) {
          z_cut <- cutpoints[[method]]
          if (!length(z_cut) || is.na(z_cut)) next

          # --- Training data ---
          train_ids <- colnames(bo_bundle_selected$data_filtered$ex)
          train_si <- tar_data$sampInfo
          train_match <- match(train_ids, rownames(train_si))
          if (any(is.na(train_match)) && "ID" %in% names(train_si)) {
            train_match <- match(train_ids, train_si$ID)
          }
          if (!all(is.na(train_match))) {
            train_si_matched <- train_si[train_match[!is.na(train_match)], ]
            train_X <- bo_bundle_selected$data_filtered$ex[, !is.na(train_match), drop = FALSE]
            if ("PurIST" %in% names(train_si_matched) && "DeCAF" %in% names(train_si_matched)) {
              train_group <- compute_factor_group(W, train_X, j, z_mean, z_sd, z_cut)
              train_df <- data.frame(
                group = train_group,
                PurIST = train_si_matched$PurIST,
                DeCAF = train_si_matched$DeCAF,
                stringsAsFactors = FALSE
              )
              train_df <- train_df[complete.cases(train_df), ]
              if (nrow(train_df) > 0) {
                p <- plot_subtype_overlap(train_df,
                  dataset_label = sprintf("Training Factor %d (%s cutpoint)", j, method))
                fpath <- file.path(out_dir, sprintf("subtype_overlap_train_factor%d_%s_%s.pdf",
                                                     j, method, bo_label))
                ggplot2::ggsave(fpath, p, width = 8, height = 5)
                paths <- c(paths, fpath)

                # Enrichment plot (confounding-aware)
                pe <- plot_subtype_enrichment(train_df,
                  dataset_label = sprintf("Training Factor %d (%s cutpoint)", j, method))
                fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_train_factor%d_%s_%s.pdf",
                                                       j, method, bo_label))
                ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
                paths <- c(paths, fpath_e)
              }
            }
          }

          # --- Per validation dataset ---
          pooled_rows <- list()
          for (ds_name in names(data_val_filtered_surv)) {
            val_ds <- data_val_filtered_surv[[ds_name]]
            common_genes <- intersect(rownames(W), rownames(val_ds$ex))
            if (length(common_genes) < 2) next

            W_common <- W[common_genes, , drop = FALSE]
            X_val <- val_ds$ex[common_genes, , drop = FALSE]

            time_val <- val_ds$sampInfo$time
            event_val <- val_ds$sampInfo$event
            valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
            if (length(valid_idx) < 2) next

            val_group <- compute_factor_group(
              W_common, X_val[, valid_idx, drop = FALSE], j, z_mean, z_sd, z_cut)
            val_si <- val_ds$sampInfo[valid_idx, ]

            if ("PurIST" %in% names(val_si) && "DeCAF" %in% names(val_si)) {
              ds_df <- data.frame(
                group = val_group,
                PurIST = val_si$PurIST,
                DeCAF = val_si$DeCAF,
                stringsAsFactors = FALSE
              )
              ds_df <- ds_df[complete.cases(ds_df), ]
              pooled_rows[[ds_name]] <- ds_df
              if (nrow(ds_df) > 0) {
                p <- plot_subtype_overlap(ds_df,
                  dataset_label = sprintf("%s Factor %d (%s cutpoint)", ds_name, j, method))
                fpath <- file.path(out_dir, sprintf("subtype_overlap_%s_factor%d_%s_%s.pdf",
                                                     ds_name, j, method, bo_label))
                ggplot2::ggsave(fpath, p, width = 8, height = 5)
                paths <- c(paths, fpath)

                # Enrichment plot (confounding-aware)
                pe <- plot_subtype_enrichment(ds_df,
                  dataset_label = sprintf("%s Factor %d (%s cutpoint)", ds_name, j, method))
                fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_%s_factor%d_%s_%s.pdf",
                                                       ds_name, j, method, bo_label))
                ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
                paths <- c(paths, fpath_e)
              }
            }
          }

          # --- Pooled validation ---
          if (length(pooled_rows) > 0) {
            pooled_df <- do.call(rbind, pooled_rows)
            if (nrow(pooled_df) > 0) {
              p <- plot_subtype_overlap(pooled_df,
                dataset_label = sprintf("Pooled Factor %d (%s cutpoint)", j, method))
              fpath <- file.path(out_dir, sprintf("subtype_overlap_pooled_factor%d_%s_%s.pdf",
                                                   j, method, bo_label))
              ggplot2::ggsave(fpath, p, width = 8, height = 5)
              paths <- c(paths, fpath)

              # Enrichment plot (confounding-aware)
              pe <- plot_subtype_enrichment(pooled_df,
                dataset_label = sprintf("Pooled Factor %d (%s cutpoint)", j, method))
              fpath_e <- file.path(out_dir, sprintf("subtype_enrichment_pooled_factor%d_%s_%s.pdf",
                                                     j, method, bo_label))
              ggplot2::ggsave(fpath_e, pe, width = 12, height = 8)
              paths <- c(paths, fpath_e)
            }
          }
        }
      }

      paths
    },
    format = "file"
  ),

  tar_target(
    val_cindex_desurv_alpha0,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_desurv_alpha0)
      if (nrow(summary_tbl)) {
        dir.create(
          file.path(val_run_bundle$training_results_dir_alpha0, "validation", val_config_effective$config_id),
          recursive = TRUE,
          showWarnings = FALSE
        )
        utils::write.csv(
          summary_tbl,
          file = file.path(
            val_run_bundle$training_results_dir_alpha0,
            "validation",
            val_config_effective$config_id,
            "val_cindex_desurv_alpha0.csv"
          ),
          row.names = FALSE
        )
      }
      summary_tbl
    }
  ),
  tar_target(
    val_cindex_desurv_elbowk,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_desurv_elbowk)
      summary_tbl
    }
  ),
  tar_target(
    val_cindex_std_elbowk,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_std_elbowk)
      summary_tbl
    }
  ),
  tar_target(
    val_cindex_std_desurvk,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_std_desurvk)
      summary_tbl
    }
  )
)

FIGURE_TARGETS <- list(

  tar_target(
    fig_bo_cvk,
    {
      p = make_bo_best_observed_plot(bo_history_path = desurv_bo_history,
                                 bo_results = desurv_bo_results,
                                 method_label = "DeSurv")
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_bo_cvk_%s.pdf", bo_label)
        )
      )
      p
    }
    
  ),
  
  tar_target(
    fig_bo_cvalpha,
    {
      browser()
      p = make_bo_best_observed_alpha_plot(bo_history_path = desurv_bo_history,
                                     bo_results = desurv_bo_results,
                                     method_label = "DeSurv")
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_bo_cvalpha_%s.pdf", bo_label)
        )
      )
      p
    }
    
  ),
  
  tar_target(
    fig_bo_heat,
    {
      curve = extract_gp_curve(desurv_bo_results,tar_params_best)
      ggplot(curve, aes(x = k, y = alpha, fill = mean)) +
        geom_tile(color = NA) +
        scale_x_continuous(
          breaks = function(x) seq(2,12, by = 1)
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
    }
  ),
  
  tar_target(
    fig_residuals,
    {
      p=make_nmf_metric_plot(fit_std, "residuals")
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_residuals_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  tar_target(
    fig_cophenetic,
    {
      p = make_nmf_metric_plot(fit_std, "cophenetic")
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_cophenetic_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  tar_target(
    fig_silhouette,
    {
      p = make_nmf_metric_plot(fit_std, "silhouette")
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_cophenetic_%s.pdf", bo_label)
        )
      )
      p
    }
  ),

  tar_target(
    fig_dotplots_desurv,
    {
      p = make_ora_dotplots(ora_analysis_desurv)
      for(i in 1:length(p)){
        if(!is.null(p[[i]])){
          save_plot_pdf(
            p[[i]],
            file.path(
              FIGURE_CONFIGS$panel_dir,
              sprintf("fig_dotplot_desurv_factor%d_%s.pdf", i, bo_label)
            )
          )
        }
      }
      p
    }
  ),
  tar_target(
    fig_gene_overlap_heatmap_desurv,
    {
      browser()
      p = make_gene_overlap_heatmap(tar_fit_desurv,tar_tops_desurv$top_genes,top_genes,
                                    factor_labels = FIGURE_CONFIGS$heatmap_factor_labels[[bo_label]])
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_gene_overlap_heatmap_desurv_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  
  tar_target(
    fig_dotplots_desurv_elbowk,
    {
      p = make_ora_dotplots(ora_analysis_desurv_elbowk)
      for(i in 1:length(p)){
        if(!is.null(p[[i]])){
          save_plot_pdf(
            p[[i]],
            file.path(
              FIGURE_CONFIGS$panel_dir,
              sprintf("fig_dotplot_desurv_elbowk_factor%d_%s.pdf", i, bo_label)
            )
          )
        }
      }
      p
    }
  ),
  tar_target(
    fig_gene_overlap_heatmap_desurv_elbowk,
    {
      p = make_gene_overlap_heatmap(tar_fit_desurv_elbowk,
                                    tar_tops_desurv_elbowk$top_genes,
                                    top_genes)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_gene_overlap_heatmap_desurv_elbowk_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  
  tar_target(
    fig_dotplots_std_elbowk,
    {
      p = make_ora_dotplots(ora_analysis_std_elbowk)
      for(i in 1:length(p)){
        if(!is.null(p[[i]])){
          save_plot_pdf(
            p[[i]],
            file.path(
              FIGURE_CONFIGS$panel_dir,
              sprintf("fig_dotplot_std_elbowk_factor%d_%s.pdf", i, bo_label)
            )
          )
        }
      }
      p
    }
  ),
  tar_target(
    fig_gene_overlap_heatmap_std_elbowk,
    {
      p = make_gene_overlap_heatmap(fit_std_elbowk,
                                    tar_tops_std_elbowk$top_genes,
                                    top_genes)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_gene_overlap_heatmap_std_elbowk_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  
  tar_target(
    fig_dotplots_std_desurvk,
    {
      p = make_ora_dotplots(ora_analysis_std_desurvk)
      for(i in 1:length(p)){
        if(!is.null(p[[i]])){
          save_plot_pdf(
            p[[i]],
            file.path(
              FIGURE_CONFIGS$panel_dir,
              sprintf("fig_dotplot_std_desurvk_factor%d_%s.pdf", i, bo_label)
            )
          )
        }
      }
      p
    }
  ),
  tar_target(
    fig_gene_overlap_heatmap_std_desurvk,
    {
      p = make_gene_overlap_heatmap(fit_std_desurvk,
                                    tar_tops_std_desurvk$top_genes,
                                    top_genes,
                                    show_legend = FALSE)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_gene_overlap_heatmap_std_desurvk_%s.pdf", bo_label)
        )
      )
      p
    }
  ),
  
  tar_target(
    fig_variation_explained,
    {
      df_nmf <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = fit_std_desurvk$W,
        loadings = fit_std_desurvk$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "NMF"
      )
      
      df_desurv <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = tar_fit_desurv$W,
        loadings = tar_fit_desurv$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "DeSurv"
      )
      
      df_plot <- dplyr::bind_rows(df_nmf, df_desurv) %>%
        dplyr::mutate(
          factor_label = dplyr::case_when(
            method == "NMF" ~ paste0("N", factor),
            method == "DeSurv" ~ paste0("D", factor),
            TRUE ~ paste0("F", factor)
          )
        )
      
      ggplot(df_plot,
             aes(x = variance_explained,
                 y = delta_loglik,
                 label = factor_label,
                 color = method)) +
        geom_point(size = 4) +
        geom_text_repel(
          size = 4,
          max.overlaps = Inf,
          box.padding = 0.4,
          point.padding = 0.3,
          segment.size = 0.3
        ) +
        scale_color_manual(
          values = c(
            "NMF" = "red",
            "DeSurv" = "blue"
          )
        ) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(
          x = "Conditional variance explained\n(semi-partial R\u00b2)",
          y = "\u0394 partial log-likelihood\n(full vs. k\u22121 factor model)",
          color = "Method"
        ) +
        theme_classic(base_size = 10)
      
    }
  ),
  
  tar_target(
    fig_desurv_std_correlation,
    {
      c = cor(fit_std_desurvk$W,tar_fit_desurv$W)
      rownames(c) = paste0("NMF F",1:ncol(c))
      colnames(c) = paste0("DeSurv F",1:ncol(c))
      ph = pheatmap::pheatmap(c,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         show_colnames = TRUE,
                         show_rownames = TRUE,
                         silent = TRUE)
      ph_grob <- ph$gtable
      pheat <- cowplot::plot_grid(NULL, cowplot::ggdraw(ph_grob), nrow = 2, rel_heights = c(0.25, 4))
      pheat
    }
  ),
  
  tar_target(
    fig_sc_panels,
    build_fig_sc_panels(
      tops_desurv = tar_tops_desurv,
      sc_all_path = FIGURE_CONFIGS$sc_data_paths$all,
      sc_caf_path = FIGURE_CONFIGS$sc_data_paths$caf,
      sc_tum_path = FIGURE_CONFIGS$sc_data_paths$tum
    ),
    packages = c("ggplot2", "dplyr", "viridis", "cowplot", "pheatmap", "ggplotify", "Seurat", "VAM")
  ),
  tar_target(fig_sc_panel_a, fig_sc_panels$A),
  tar_target(fig_sc_panel_b, fig_sc_panels$B),
  tar_target(fig_sc_panel_c, fig_sc_panels$C),
  tar_target(fig_sc_panel_d, fig_sc_panels$D),
  tar_target(fig_sc_panel_e, fig_sc_panels$E),
  tar_target(
    fig_sc_panel_a_file,
    save_plot_pdf(
      fig_sc_panel_a,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_sc_%s_panel_a.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_sc_panel_b_file,
    save_plot_pdf(
      fig_sc_panel_b,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_sc_%s_panel_b.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_sc_panel_c_file,
    save_plot_pdf(
      fig_sc_panel_c,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_sc_%s_panel_c.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_sc_panel_d_file,
    save_plot_pdf(
      fig_sc_panel_d,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_sc_%s_panel_d.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_sc_panel_e_file,
    save_plot_pdf(
      fig_sc_panel_e,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_sc_%s_panel_e.pdf", bo_label)
      )
    ),
    format = "file"
  )
)

FIGURE_VAL_TARGETS <- list(
  tar_target(
    fig_extval_panels,
    build_fig_extval_panels(
      data_val_filtered = data_val_filtered,
      fit_desurv = tar_fit_desurv,
      tops_desurv = tar_tops_desurv
    ),
    packages = c("ggplot2", "dplyr", "pheatmap", "cowplot", "survival", "survminer")
  ),
  tar_target(fig_extval_panel_a, fig_extval_panels$A),
  tar_target(fig_extval_panel_b, fig_extval_panels$B),
  tar_target(fig_extval_panel_c, fig_extval_panels$C),
  tar_target(
    fig_extval_panel_a_file,
    save_plot_pdf(
      fig_extval_panel_a,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_extval_%s_panel_a.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_extval_panel_b_file,
    save_plot_pdf(
      fig_extval_panel_b,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_extval_%s_panel_b.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_extval_panel_c_file,
    save_plot_pdf(
      fig_extval_panel_c,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_extval_%s_panel_c.pdf", bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_hr_forest,
    {
      desurv_df = compute_hrs(data_val_filtered,tar_fit_desurv,"DeSurv")
      nmf_df = compute_hrs(data_val_filtered,fit_std_desurvk,"NMF")
      
      nmf_var <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = fit_std_desurvk$W,
        loadings = fit_std_desurvk$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "Std NMF"
      )
      
      desurv_var <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = tar_fit_desurv$W,
        loadings = tar_fit_desurv$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "DeSurv"
      )
      
      nmf_fac = nmf_var$factor[which.max(nmf_var$delta_loglik)]
      desurv_fac = desurv_var$factor[which.max(desurv_var$delta_loglik)]
      
      df = rbind(desurv_df,nmf_df)
      df_select = df #%>% filter(
        # (factor==nmf_fac & method=="NMF") | (factor==desurv_fac & method=="DeSurv"))
      pd <- position_dodge(width = 0.6)
      
      df_select$label <- sprintf("%.2f (%.2f–%.2f)", df_select$HR, df_select$lower, df_select$upper)
      
      ggplot(df_select, aes(x = HR, y = as.factor(factor),color=dataset,group=dataset)) +
        geom_vline(xintercept = 1, linetype = "dashed",
                   linewidth = 0.5, color = "grey60") +
        geom_errorbarh(
          aes(xmin = lower, xmax = upper),
          height = 0.25,
          linewidth = 0.8,
          position=pd
        ) +
        geom_point(size = 2.8,position = pd) +
        # geom_text(
        #   aes(x = max(upper) * 1.05, y=factor,label = label,group=dataset),
        #   inherit.aes = FALSE,
        #   hjust = 0,
        #   size = 3,
        #   position = pd,
        #   show.legend = FALSE
        # ) +
        scale_x_log10(
          breaks = c(0.5, 1, 2, 4),
          labels = c("0.5", "1", "2", "4")
        ) +
        theme_classic(base_size = 12) +
        theme(
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()
        ) +
        coord_cartesian(xlim = c(min(df_select$lower), max(df_select$upper) * 1.4)) +
        labs(
          x = "Hazard ratio (95% CI)",
          y = "Factor label"
        )+facet_wrap(~method)
    }
  ),
  tar_target(
    fig_median_survival_desurv,
    {
      
      desurv_var <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = tar_fit_desurv$W,
        loadings = tar_fit_desurv$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "DeSurv"
      )
      
      desurv_fac = desurv_var$factor[which.max(desurv_var$delta_loglik)]
      
      splot_median(data_val_filtered,tar_fit_desurv,desurv_fac)
    }
  ),
  tar_target(
    fig_median_survival_std_desurvk,
    {
      nmf_var <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = fit_std_desurvk$W,
        loadings = fit_std_desurvk$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "Std NMF"
      )
      
      nmf_fac = nmf_var$factor[which.max(nmf_var$delta_loglik)]
      
      splot_median(data_val_filtered,fit_std_desurvk,nmf_fac)
    }
  )
)
