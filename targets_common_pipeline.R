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
    )
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
    )
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
    )
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
    )
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
    )
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
    )
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
    )
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
              transform_target = val_run_bundle$bo_bundle$data_filtered$transform_target
            )
          }
        ),
        names(datasets_named)
      )
    },
    iteration="list"
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
              transform_target = val_run_bundle$bo_bundle$data_filtered$transform_target
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
  ),
  
  #clustering
  tar_target(
    clusters_desurv_X,
    {
      beta = tar_fit_desurv$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "desurv"
      )
      run_clustering(tops = tar_tops_desurv$top_genes,
                     data = data_val_filtered,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = FALSE)
    },
    pattern = map(data_val_filtered),
    iteration = "list"
  ),
  tar_target(
    nclusters_desurv_X,
    {
      sel = select_nclusters(clusters_desurv_X$clus,k_max=length(clusters_desurv_X$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_desurv_X)
  ),
  tar_target(
    clusters_desurv_X_aligned,
    {
      cluster_list = lapply(1:length(clusters_desurv_X),function(i){
        clusters_desurv_X[[i]]$clus[[nclusters_desurv_X[i]]]$consensusClass
      })
      scores_list = lapply(1:length(clusters_desurv_X),function(i){
        t(clusters_desurv_X[[i]]$Xtemp)
      })

      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                         linkage = "average",
                         similarity_threshold = .5,
                         zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_desurv_WtX,
    {
      beta = tar_fit_desurv$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "desurv"
      )
      run_clustering(tops = tar_tops_desurv$top_genes,
                     data = val_latent_desurv,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = TRUE)
    },
    pattern = map(val_latent_desurv),
    iteration = "list"
  ),
  tar_target(
    nclusters_desurv_WtX,
    {
      sel = select_nclusters(clusters_desurv_WtX$clus,k_max=length(clusters_desurv_WtX$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_desurv_WtX)
  ),
  tar_target(
    clusters_desurv_WtX_aligned,
    {
      cluster_list = lapply(1:length(clusters_desurv_WtX),function(i){
        clusters_desurv_WtX[[i]]$clus[[nclusters_desurv_WtX[i]]]$consensusClass
      })
      scores_list = lapply(1:length(val_latent_desurv),function(i){
        val_latent_desurv[[i]]$Z_scaled
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_desurv_elbowk_X,
    {
      beta = tar_fit_desurv_elbowk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "desurv_elbowk"
      )
      run_clustering(tops = tar_tops_desurv_elbowk$top_genes,
                     data = data_val_filtered_elbowk,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = FALSE)
    },
    pattern = map(data_val_filtered_elbowk),
    iteration = "list"
  ),
  tar_target(
    nclusters_desurv_elbowk_X,
    {
      sel = select_nclusters(clusters_desurv_elbowk_X$clus,k_max=length(clusters_desurv_elbowk_X$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_desurv_elbowk_X)
  ),
  tar_target(
    clusters_desurv_elbowk_X_aligned,
    {
      cluster_list = lapply(1:length(clusters_desurv_elbowk_X),function(i){
        clusters_desurv_elbowk_X[[i]]$clus[[nclusters_desurv_elbowk_X[i]]]$consensusClass
      })
      scores_list = lapply(1:length(clusters_desurv_elbowk_X),function(i){
        t(clusters_desurv_elbowk_X[[i]]$Xtemp)
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_desurv_elbowk_WtX,
    {
      beta = tar_fit_desurv_elbowk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "desurv_elbowk"
      )
      run_clustering(tops = tar_tops_desurv_elbowk$top_genes,
                     data = val_latent_desurv_elbowk,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = TRUE)
    },
    pattern = map(val_latent_desurv_elbowk),
    iteration = "list"
  ),
  tar_target(
    nclusters_desurv_elbowk_WtX,
    {
      sel = select_nclusters(clusters_desurv_elbowk_WtX$clus,
                             k_max=length(clusters_desurv_elbowk_WtX$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_desurv_elbowk_WtX)
  ),
  tar_target(
    clusters_desurv_elbowk_WtX_aligned,
    {
      cluster_list = lapply(1:length(clusters_desurv_elbowk_WtX),function(i){
        clusters_desurv_elbowk_WtX[[i]]$clus[[nclusters_desurv_elbowk_WtX[i]]]$consensusClass
      })
      scores_list = lapply(1:length(val_latent_desurv_elbowk),function(i){
        val_latent_desurv_elbowk[[i]]$Z_scaled
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_std_elbowk_X,
    {
      beta = fit_std_elbowk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "std_elbowk"
      )
      run_clustering(tops = tar_tops_std_elbowk$top_genes,
                     data = data_val_filtered_elbowk,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = FALSE)
    },
    pattern = map(data_val_filtered_elbowk),
    iteration = "list"
  ),
  tar_target(
    nclusters_std_elbowk_X,
    {
      sel = select_nclusters(clusters_std_elbowk_X$clus,k_max=length(clusters_std_elbowk_X$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_std_elbowk_X)
  ),
  tar_target(
    clusters_std_elbowk_X_aligned,
    {
      cluster_list = lapply(1:length(clusters_std_elbowk_X),function(i){
        clusters_std_elbowk_X[[i]]$clus[[nclusters_std_elbowk_X[i]]]$consensusClass
      })
      scores_list = lapply(1:length(clusters_std_elbowk_X),function(i){
        t(clusters_std_elbowk_X[[i]]$Xtemp)
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_std_elbowk_WtX,
    {
      beta = fit_std_elbowk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "std_elbowk"
      )
      run_clustering(tops = tar_tops_std_elbowk$top_genes,
                     data = val_latent_std_elbowk,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = TRUE)
    },
    pattern = map(val_latent_std_elbowk),
    iteration = "list"
  ),
  tar_target(
    nclusters_std_elbowk_WtX,
    {
      sel = select_nclusters(clusters_std_elbowk_WtX$clus,
                             k_max=length(clusters_std_elbowk_WtX$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_std_elbowk_WtX)
  ),
  tar_target(
    clusters_std_elbowk_WtX_aligned,
    {
      cluster_list = lapply(1:length(clusters_std_elbowk_WtX),function(i){
        clusters_std_elbowk_WtX[[i]]$clus[[nclusters_std_elbowk_WtX[i]]]$consensusClass
      })
      scores_list = lapply(1:length(val_latent_std_elbowk),function(i){
        val_latent_std_elbowk[[i]]$Z_scaled
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_std_desurvk_X,
    {
      beta = fit_std_desurvk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "std_desurvk"
      )
      run_clustering(tops = tar_tops_std_desurvk$top_genes,
                     data = data_val_filtered,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = FALSE)
    },
    pattern = map(data_val_filtered),
    iteration = "list"
  ),
  tar_target(
    nclusters_std_desurvk_X,
    {
      sel = select_nclusters(clusters_std_desurvk_X$clus,k_max=length(clusters_std_desurvk_X$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_std_desurvk_X)
  ),
  tar_target(
    clusters_std_desurvk_X_aligned,
    {
      cluster_list = lapply(1:length(clusters_std_desurvk_X),function(i){
        clusters_std_desurvk_X[[i]]$clus[[nclusters_std_desurvk_X[i]]]$consensusClass
      })
      scores_list = lapply(1:length(clusters_std_desurvk_X),function(i){
        t(clusters_std_desurvk_X[[i]]$Xtemp)
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
    }
  ),
  
  tar_target(
    clusters_std_desurvk_WtX,
    {
      beta = fit_std_desurvk$beta
      facs = which(beta != 0)
      base_dir = file.path(
        val_run_bundle$training_results_dir,
        "validation",
        val_config_effective$config_id,
        "std_desurvk"
      )
      run_clustering(tops = tar_tops_std_desurvk$top_genes,
                     data = val_latent_std_desurvk,
                     gene_lists = top_genes,
                     color.lists = colors,
                     facs = facs,
                     base_dir = base_dir,
                     WtX = TRUE)
    },
    pattern = map(val_latent_std_desurvk),
    iteration = "list"
  ),
  tar_target(
    nclusters_std_desurvk_WtX,
    {
      sel = select_nclusters(clusters_std_desurvk_WtX$clus,
                             k_max=length(clusters_std_desurvk_WtX$clus))
      sel$k
    },
    iteration = "vector",
    pattern = map(clusters_std_desurvk_WtX)
  ),
  tar_target(
    clusters_std_desurvk_WtX_aligned,
    {
      cluster_list = lapply(1:length(clusters_std_desurvk_WtX),function(i){
        clusters_std_desurvk_WtX[[i]]$clus[[nclusters_std_desurvk_WtX[i]]]$consensusClass
      })
      scores_list = lapply(1:length(val_latent_std_desurvk),function(i){
        val_latent_std_desurvk[[i]]$Z_scaled
      })
      
      temp=meta_cluster_align(scores_list,cluster_list,similarity = 'cosine',
                              linkage = "average",
                              similarity_threshold = .5,
                              zscore_within_dataset = TRUE)
      temp
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
        labs(
          x = "Number of components (k)",
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
      p = make_gene_overlap_heatmap(tar_fit_desurv,tar_tops_desurv$top_genes,top_genes)
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
                                    top_genes)
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
        method = "Std NMF"
      )
      
      df_desurv <- build_variance_survival_df(
        X = tar_data_filtered$ex,
        scores = tar_fit_desurv$W,
        loadings = tar_fit_desurv$H,
        time = tar_data_filtered$sampInfo$time,
        event = tar_data_filtered$sampInfo$event,
        method = "DeSurv"
      )
      
      df_plot = rbind(df_nmf,df_desurv)
      
      ggplot(df_plot,
             aes(x = variance_explained,
                 y = delta_loglik,
                 label = paste0("F", factor),
                 shape = method)) +
        geom_point(size = 4) +
        geom_text_repel(
          size = 4,
          max.overlaps = Inf,
          box.padding = 0.4,
          point.padding = 0.3,
          segment.size = 0.3
        ) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(
          x = "Fraction of expression variance explained",
          y = "Δ partial log-likelihood (survival)",
          shape = "Method"
        ) +
        theme_classic(base_size = 12)
      
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
    fig_validation_heatmap_desurv,
    {
      browser()
      p = make_expression_heatmap(data_val_filtered,
                              tar_fit_desurv,
                              tar_tops_desurv,
                              clusters_desurv_X_aligned,
                              clusters_desurv_X,
                              nclusters_desurv_X)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_val_heatmap_desurv_%s.pdf", bo_label)
        )
      )
      p
    }

  ),
  
  tar_target(
    fig_validation_heatmap_desurv_elbowk,
    {
      browser()
      p = make_expression_heatmap(data_val_filtered_elbowk,
                              tar_fit_desurv_elbowk,
                              tar_tops_desurv_elbowk,
                              clusters_desurv_elbowk_X_aligned,
                              clusters_desurv_elbowk_X,
                              nclusters_desurv_elbowk_X)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_val_heatmap_desurv_elbowk_%s.pdf", bo_label)
        )
      )
      p
    }
    
  ),
  
  tar_target(
    fig_validation_heatmap_std_elbowk,
    {
      browser()
      p = make_expression_heatmap(data_val_filtered_elbowk,
                                  fit_std_elbowk,
                                  tar_tops_std_elbowk,
                                  clusters_std_elbowk_X_aligned,
                                  clusters_std_elbowk_X,
                                  nclusters_std_elbowk_X)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_val_heatmap_std_elbowk_%s.pdf", bo_label)
        )
      )
      p
    }
    
  ),
  
  tar_target(
    fig_validation_heatmap_std_desurvk,
    {
      browser()
      p = make_expression_heatmap(data_val_filtered,
                                  fit_std_desurvk,
                                  tar_tops_std_desurvk,
                                  clusters_std_desurvk_X_aligned,
                                  clusters_std_desurvk_X,
                                  nclusters_std_desurvk_X)
      save_plot_pdf(
        p,
        file.path(
          FIGURE_CONFIGS$panel_dir,
          sprintf("fig_val_heatmap_std_desurvk_%s.pdf", bo_label)
        )
      )
      p
    }
    
  ),
  
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
  )
)
