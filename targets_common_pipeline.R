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
    tar_params_best,
    standardize_bo_params(desurv_bo_results$overall_best$params)
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
    tar_params_best_alpha0,
    {
      params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
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
    fit_std,
    nmf(
      bo_bundle_selected$data_filtered$ex,
      run_config$std_nmf_k_grid,
      nrun = run_config$std_nmf_nrun,
      method = "lee",
      .options = paste0("p", run_config$std_nmf_nrun)
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
  
  # tar_target(
  #   std_nmf_k_selection_table,
  #   {
  #     save_dir=dirname(std_nmf_k_selection_plots)
  #     path <- file.path(save_dir,
  #                       paste0("std_nmf_k_selection_table.csv")
  #                       )
  #     if(!file.exists(std_nmf_k_selection_plots)){
  #       stop("Standard NMF k selection plots not generated")
  #     }
  #     build_std_nmf_k_selection_table(STD_NMF_K_GRID,path,fit_std)
  #   },
  # ),
  
  # tar_target(
  #   fit_std_selected_k,
  #   {
  #     std_nmf_k_table = read.csv(std_nmf_k_selection_table,
  #                                stringsAsFactors = FALSE)
  #     if(all(std_nmf_k_table$selected==FALSE)){
  #       stop("Please edit std_nmf_k_selection.csv to select a k")
  #     }
  #     k_selected = std_nmf_k_table$rank[which(std_nmf_k_table$selected)]
  #     print(paste0("You have selected ",k_selected," factors for standard NMF"))
  #     fit_std$fit[[as.character(k_selected)]]
  #   }
  # ),
  # 
  # 
  # tar_target(
  #   fit_std_beta,
  #   {
  #     fit_list <- fit_std$fit
  #     if (is.null(fit_list) || !length(fit_list)) {
  #       stop("fit_std does not contain any fits to evaluate")
  #     }
  #     k_labels <- names(fit_list)
  #     if (is.null(k_labels) || any(!nzchar(k_labels))) {
  #       k_labels <- as.character(seq_along(fit_list))
  #       names(fit_list) <- k_labels
  #     }
  #     X <- tar_data_filtered$ex
  #     y <- tar_data_filtered$sampInfo$time
  #     d <- tar_data_filtered$sampInfo$event
  #     strata <- interaction(d, tar_data_filtered$sampInfo$dataset, drop = FALSE)
  #     foldid <- caret::createFolds(strata, NFOLD, list = FALSE)
  #     results <- setNames(vector("list", length(fit_list)), k_labels)
  #     for (k_label in k_labels) {
  #       nmf_fit <- fit_list[[k_label]]
  #       W <- nmf_fit@fit@W
  #       H <- nmf_fit@fit@H
  #       genes <- intersect(rownames(W), rownames(X))
  #       W_use <- W[genes, , drop = FALSE]
  #       X_use <- X[genes, , drop = FALSE]
  #       Z <- t(X_use) %*% W_use
  #       gfit <- list()
  #       mets <- vector("list", length(COXNET_ALPHA_GRID))
  #       for (i in seq_along(COXNET_ALPHA_GRID)) {
  #         alpha <- COXNET_ALPHA_GRID[[i]]
  #         alpha_chr <- as.character(alpha)
  #         cv_fit <- cv.glmnet(
  #           Z,
  #           Surv(y, d),
  #           family = "cox",
  #           type.measure = "C",
  #           alpha = alpha,
  #           lambda = COXNET_LAMBDA_GRID,
  #           foldid = foldid
  #         )
  #         gfit[[alpha_chr]] <- cv_fit
  #         ind <- which(cv_fit$lambda == cv_fit$lambda.1se)
  #         mets[[i]] <- data.frame(
  #           eta = alpha,
  #           lambda = cv_fit$lambda.1se,
  #           cvm = cv_fit$cvm[ind]
  #         )
  #       }
  #       mets_all <- dplyr::bind_rows(mets)
  #       mets_best <- mets_all[which.max(mets_all$cvm), , drop = FALSE]
  #       best_fit <- gfit[[as.character(mets_best$eta)]]
  #       beta <- coef(best_fit, s = best_fit$lambda.1se)
  #       results[[k_label]] <- list(
  #         k = suppressWarnings(as.numeric(k_label)),
  #         W = W_use,
  #         H = H,
  #         beta = beta,
  #         metrics = mets_all,
  #         alpha = mets_best$eta,
  #         lambda = mets_best$lambda
  #       )
  #     }
  #     results
  #   }
  # ),
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
  )
  
  # tar_target(
  #   ora_analysis_std,
  #   {
  #     universe = rownames(tar_data_filtered$ex)
  #     organism <- org.Hs.eg.db
  #     ora(tops_std$top_genes,universe,organism)
  #   }
  # ),
  
  # 
  # 
  # tar_target(
  #   selected_factors_desurv,
  #   {
  #     save_dir = file.path(tar_training_results_dir,"factor_selection","desurv")
  #     dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
  #     path = file.path(save_dir,"factor_selection.csv")
  #     build_factor_selection_table(tar_tops_desurv$top_genes,path)
  #   },
  #   format="file"
  # ),

  # tar_target(
  #   selected_factors_desurv_alpha0,
  #   {
  #     save_dir = file.path(tar_training_results_dir_alpha0,"factor_selection","desurv_alpha0")
  #     dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
  #     path = file.path(save_dir,"factor_selection.csv")
  #     build_factor_selection_table(tar_tops_desurv_alpha0$top_genes,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   selected_factors_std,
  #   {
  #     save_dir = file.path(tar_training_results_dir,"factor_selection","std")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"factor_selection.csv")
  #     build_factor_selection_table(tops_std$top_genes,path)
  #   },
  #   format="file"
  # ),
    # tar_target(
    #   clusters_desurv,
    #   {
    #     sel = read_selected_factor_indices(
    #       selected_factors_desurv,
    #       label = "DeSurv factor selection table"
    #     )
    #     data <- unwrap_validation_dataset(data_val_filtered)
    #     val_dataset <- infer_validation_dataset_name(data)
    #     dir = create_filepath_clustering_output(tar_ngene_value, TOL, MAXIT, PKG_VERSION,
    #                                             GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN,
    #                                             ntop_value, val_dataset, "DeSurv")
    #     run_clustering(
    #       tar_tops_desurv$top_genes, data, top_genes, colors,
    #       facs = sel, plot = FALSE, dir = dir,
    #       maxKcol = 5, maxKrow = 5
    #     )
    #   },
    #   pattern = map(data_val_filtered),
    #   iteration = "list",
    #   resources = tar_resources(
    #     crew = tar_resources_crew(controller = "med_mem")
    #   )
    # ),
    # 
    # tar_target(
    #   clusters_desurv_alpha0,
    #   {
    #     sel = read_selected_factor_indices(
    #       selected_factors_desurv_alpha0,
    #       label = "DeSurv alpha=0 factor selection table"
    #     )
    #     data <- unwrap_validation_dataset(data_val_filtered)
    #     val_dataset <- infer_validation_dataset_name(data)
    #     dir = create_filepath_clustering_output(tar_ngene_value_alpha0, TOL, MAXIT, PKG_VERSION,
    #                                             GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN,
    #                                             ntop_value_alpha0, val_dataset, "DeSurv_alpha0")
    #     run_clustering(
    #       tar_tops_desurv_alpha0$top_genes, data, top_genes, colors,
    #       facs = sel, plot = FALSE, dir = dir,
    #       maxKcol = 5, maxKrow = 5
    #     )
    #   },
    #   pattern = map(data_val_filtered),
    #   iteration = "list",
    #   resources = tar_resources(
    #     crew = tar_resources_crew(controller = "med_mem")
    #   )
    # )

  #   tar_target(
  #     clusters_std,
  #     {
  #       if(!file.exists(selected_factors_std)){
  #         stop("first generate factor selection table for std nmf")
  #       }
  #       tbl = read.csv(selected_factors_std)
  #       sel = tbl$factor[tbl$selected]
  #       if(all(!sel)){
  #         stop("No factor selected for std nmf clustering, please select at least 1")
  #       }
  #       data <- data_val_filtered
  #       val_dataset = data$dataname
  #       dir = create_filepath_clustering_output(tar_ngene_value, TOL, MAXIT, PKG_VERSION, 
  #                                               GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN, 
  #                                               ntop_value, val_dataset, "stdNMF")
  #       run_clustering(tops_std$top_genes,data,top_genes,colors,
  #                      facs=sel,plot=FALSE,dir=dir,
  #                    maxKcol = 5, maxKrow = 5)
  #   },
  #   pattern=map(data_val_filtered),
  #   iteration="list",
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "med_mem")
  #   )
  # ),
  # 
  # tar_target(
  #   selected_nclusters_desurv,
  #   {
  #     save_dir = file.path(tar_training_results_dir,"ncluster_selection","desurv")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"ncluster_selection.csv")
  #     build_ncluster_selection_table(clusters_desurv,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   selected_nclusters_std,
  #   {
  #     save_dir = file.path(tar_training_results_dir,"ncluster_selection","std")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"ncluster_selection.csv")
  #     build_ncluster_selection_table(clusters_std,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   cluster_alignment_plot_desurv,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_desurv)
  #     nclus_vec <- setNames(nclus_tbl$nclus, nclus_tbl$dataset)
  #     save_dir = file.path(tar_training_results_dir,"cluster_alignment","desurv")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"cluster_alignment.pdf")
  #     pdf(path)
  #     plot_cluster_cor(clusters_desurv,tar_tops_desurv$top_genes,nclus_vec)
  #     dev.off()
  #   }
  # ),
  # 
  # tar_target(
  #   cluster_alignment_plot_std,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_std)
  #     nclus_vec <- setNames(nclus_tbl$nclus, nclus_tbl$dataset)
  #     save_dir = file.path(tar_training_results_dir,"cluster_alignment","std")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"cluster_alignment.pdf")
  #     pdf(path)
  #     plot_cluster_cor(clusters_std,tops_std$top_genes,nclus_vec)
  #     dev.off()
  #   }
  # ),
  # 
  # tar_target(
  #   cluster_alignment_table_desurv,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_desurv)
  #     save_dir = file.path(tar_training_results_dir,"cluster_alignment","desurv")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"cluster_alignment.csv")
  #     build_cluster_alignment_table(nclus_tbl,clusters_desurv,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   cluster_alignment_table_std,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_std)
  #     save_dir = file.path(tar_training_results_dir,"cluster_alignment","std")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"cluster_alignment.csv")
  #     build_cluster_alignment_table(nclus_tbl,clusters_std,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   aligned_clusters_desurv,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_desurv)
  #     samp_clus = read.csv(cluster_alignment_table_desurv)
  #     clus=clusters_desurv
  #     datasets <- vapply(clus, function(x) x$data$dataname, character(1))
  #     for(i in seq_along(clus)){
  #       dataname=datasets[[i]]
  #       nclus_i <- nclus_tbl$nclus[nclus_tbl$dataset==dataname]
  #       if(length(nclus_i) != 1 || is.na(nclus_i)){
  #         stop("Selected ncluster table missing entry for dataset: ", dataname)
  #       }
  #       nms = names(clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass)
  #       y = samp_clus[clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass,dataname]
  #       names(y) = nms
  #       clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass=y
  # 
  #       clus[[i]]$data$sampInfo$samp_cluster = clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass
  #     }
  #     clus
  #   }
  # ),
  # 
  # tar_target(
  #   aligned_clusters_std,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_std)
  #     samp_clus = read.csv(cluster_alignment_table_std)
  #      clus=clusters_std
  #     datasets <- vapply(clus, function(x) x$data$dataname, character(1))
  #     for(i in seq_along(clus)){
  #       dataname=datasets[[i]]
  #       nclus_i <- nclus_tbl$nclus[nclus_tbl$dataset==dataname]
  #       if(length(nclus_i) != 1 || is.na(nclus_i)){
  #         stop("Selected ncluster table missing entry for dataset: ", dataname)
  #       }
  #       nms = names(clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass)
  #       y = samp_clus[clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass,dataname]
  #       names(y) = nms
  #       clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass=y
  #       
  #       clus[[i]]$data$sampInfo$samp_cluster = clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass
  #     }
  #     clus
  #   }
  # )
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
    }
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
    }
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
  )
  # tar_target(
  #   val_clusters_desurv,
  #   {
  #     base_dir <- file.path(
  #       val_run_bundle$training_results_dir,
  #       "validation",
  #       val_config$config_id,
  #       "desurv",
  #       "clusters"
  #     )
  #     dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  #     lapply(
  #       val_latent_desurv,
  #       function(entry) {
  #         dataset_dir <- file.path(base_dir, entry$dataset)
  #         dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  #         result <- cluster_validation_latent(
  #           latent_entry = entry,
  #           maxK = val_config$val_cluster_maxk,
  #           reps = val_config$val_cluster_reps,
  #           pItem = val_config$val_cluster_pitem,
  #           pFeature = val_config$val_cluster_pfeature,
  #           seed = val_config$val_cluster_seed,
  #           dir = dataset_dir
  #         )
  #         if (nrow(result$assignments)) {
  #           utils::write.csv(
  #             result$assignments,
  #             file = file.path(dataset_dir, "cluster_assignments.csv"),
  #             row.names = FALSE
  #           )
  #         }
  #         result
  #       }
  #     )
  #   },
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "med_mem")
  #   )
  # ),
  # tar_target(
  #   val_clusters_desurv_alpha0,
  #   {
  #     base_dir <- file.path(
  #       val_run_bundle$training_results_dir_alpha0,
  #       "validation",
  #       val_config$config_id,
  #       "desurv_alpha0",
  #       "clusters"
  #     )
  #     dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  #     lapply(
  #       val_latent_desurv_alpha0,
  #       function(entry) {
  #         dataset_dir <- file.path(base_dir, entry$dataset)
  #         dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  #         result <- cluster_validation_latent(
  #           latent_entry = entry,
  #           maxK = val_config$val_cluster_maxk,
  #           reps = val_config$val_cluster_reps,
  #           pItem = val_config$val_cluster_pitem,
  #           pFeature = val_config$val_cluster_pfeature,
  #           seed = val_config$val_cluster_seed,
  #           dir = dataset_dir
  #         )
  #         if (nrow(result$assignments)) {
  #           utils::write.csv(
  #             result$assignments,
  #             file = file.path(dataset_dir, "cluster_assignments.csv"),
  #             row.names = FALSE
  #           )
  #         }
  #         result
  #       }
  #     )
  #   },
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "med_mem")
  #   )
  # )
)

FIGURE_TARGETS <- list(
  tar_target(
    fig_bo_panels,
    build_fig_bo_panels(
      bo_history_path = desurv_bo_history,
      bo_history_alpha0_path = desurv_bo_history_alpha0,
      bo_results_supervised = desurv_bo_results,
      bo_results_alpha0 = desurv_bo_results_alpha0,
      fit_std = fit_std
    ),
    packages = c("ggplot2", "dplyr", "cowplot", "tibble", "DiceKriging", "NMF")
  ),
  tar_target(fig_bo_panel_a, fig_bo_panels$A),
  tar_target(fig_bo_panel_b, fig_bo_panels$B),
  tar_target(fig_bo_panel_c, fig_bo_panels$C),
  tar_target(fig_bo_panel_d, fig_bo_panels$D),
  tar_target(fig_bo_panel_e, fig_bo_panels$E),
  tar_target(fig_bo_panel_f, fig_bo_panels$F),
  tar_target(fig_bo_panel_g, fig_bo_panels$G),
  tar_target(
    fig_bo_panel_a_file,
    save_plot_pdf(
      fig_bo_panel_a,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_a.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_b_file,
    save_plot_pdf(
      fig_bo_panel_b,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_b.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_c_file,
    save_plot_pdf(
      fig_bo_panel_c,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_c.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_d_file,
    save_plot_pdf(
      fig_bo_panel_d,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_d.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_e_file,
    save_plot_pdf(
      fig_bo_panel_e,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_e.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_f_file,
    save_plot_pdf(
      fig_bo_panel_f,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_f.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_panel_g_file,
    save_plot_pdf(
      fig_bo_panel_g,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bo_%s_%s_panel_g.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bo_plot,
    combine_fig_bo_panels(fig_bo_panels),
    packages = c("cowplot")
  ),
  tar_target(
    fig_bo,
    save_plot_pdf(
      fig_bo_plot,
      file.path(
        FIGURE_CONFIGS$figures_dir,
        sprintf("fig_bo_%s_%s.pdf", run_label, bo_label)
      ),
      width = 6,
      height = 5.5
    ),
    format = "file"
  ),
  tar_target(
    fig_bio_bundle,
    build_fig_bio_panels(
      ora_analysis = ora_analysis_desurv,
      fit_desurv = tar_fit_desurv,
      tops_desurv = tar_tops_desurv,
      top_genes_ref = top_genes
    ),
    packages = c("ggplot2", "dplyr", "stringr", "viridis", "cowplot", "purrr", "enrichplot", "pheatmap")
  ),
  tar_target(fig_bio_panels, fig_bio_bundle$panels),
  tar_target(fig_bio_panel_a, fig_bio_panels$A),
  tar_target(fig_bio_panel_b, fig_bio_panels$B),
  tar_target(fig_bio_panel_c, fig_bio_panels$C),
  tar_target(fig_bio_panel_d, fig_bio_panels$D),
  tar_target(
    fig_bio_panel_a_file,
    save_plot_pdf(
      fig_bio_panel_a,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bio_%s_%s_panel_a.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bio_panel_b_file,
    save_plot_pdf(
      fig_bio_panel_b,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bio_%s_%s_panel_b.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bio_panel_c_file,
    save_plot_pdf(
      fig_bio_panel_c,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bio_%s_%s_panel_c.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bio_panel_d_file,
    save_plot_pdf(
      fig_bio_panel_d,
      file.path(
        FIGURE_CONFIGS$panel_dir,
        sprintf("fig_bio_%s_%s_panel_d.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_bio_plot,
    combine_fig_bio_panels(fig_bio_bundle),
    packages = c("cowplot")
  ),
  tar_target(
    fig_bio,
    save_plot_pdf(
      fig_bio_plot,
      file.path(
        FIGURE_CONFIGS$figures_dir,
        sprintf("fig_bio_%s_%s.pdf", run_label, bo_label)
      ),
      width = 7,
      height = 4.5
    ),
    format = "file"
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
        sprintf("fig_sc_%s_%s_panel_a.pdf", run_label, bo_label)
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
        sprintf("fig_sc_%s_%s_panel_b.pdf", run_label, bo_label)
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
        sprintf("fig_sc_%s_%s_panel_c.pdf", run_label, bo_label)
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
        sprintf("fig_sc_%s_%s_panel_d.pdf", run_label, bo_label)
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
        sprintf("fig_sc_%s_%s_panel_e.pdf", run_label, bo_label)
      )
    ),
    format = "file"
  ),
  tar_target(
    fig_sc_plot,
    combine_fig_sc_panels(fig_sc_panels),
    packages = c("cowplot")
  ),
  tar_target(
    fig_sc,
    save_plot_pdf(
      fig_sc_plot,
      file.path(
        FIGURE_CONFIGS$figures_dir,
        sprintf("fig_sc_%s_%s.pdf", run_label, bo_label)
      ),
      width = 7.5,
      height = 8
    ),
    format = "file"
  )
)

FIGURE_VAL_TARGETS <- list(
  # tar_target(
  #   fig_clus,
  #   save_fig_clus(
  #     data_val_filtered = data_val_filtered,
  #     fit_desurv = tar_fit_desurv,
  #     tops_desurv = tar_tops_desurv,
  #     path = file.path(
  #       FIGURE_CONFIGS$figures_dir,
  #       sprintf("fig_clus_%s_%s_%s.pdf", val_label, run_label, bo_label)
  #     )
  #   ),
  #   format = "file",
  #   packages = c(
  #     "ggplot2",
  #     "dplyr",
  #     "cowplot",
  #     "pheatmap",
  #     "survival",
  #     "survminer",
  #     "gt",
  #     "magick",
  #     "ConsensusClusterPlus"
  #   )
  # )
)
