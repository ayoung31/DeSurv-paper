COMMON_DESURV_TARGETS <- list(
  tar_target(
    desurv_bo_results,
    {
      bounds <- DESURV_BO_BOUNDS
      bounds <- maybe_add_numeric_bound(bounds, NGENE_CONFIG, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, NTOP_CONFIG, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAW_CONFIG, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAH_CONFIG, "lambdaH_grid", log_scale = TRUE)

      bo_fixed <- list(n_starts = NINIT)
      if (!TUNE_NGENE) {
        bo_fixed$ngene <- NGENE_DEFAULT
      }
      if (!TUNE_NTOP) {
        bo_fixed$ntop <- NTOP_DEFAULT
      }
      if (!TUNE_LAMBDAW) {
        bo_fixed$lambdaW_grid <- LAMBDAW_DEFAULT
      }
      if (!TUNE_LAMBDAH) {
        bo_fixed$lambdaH_grid <- LAMBDAH_DEFAULT
      }

      DeSurv::desurv_cv_bayesopt_refine(
        X = data$ex,
        y = data$sampInfo$time,
        d = data$sampInfo$event,
        dataset = data$sampInfo$dataset,
        samp_keeps = data$samp_keeps,
        preprocess = TRUE,
        method_trans_train = METHOD_TRANS_TRAIN,
        engine = "warmstart",
        nfolds = NFOLD,
        tol = TOL,
        maxit = MAXIT,
        coarse_bounds = bounds,
        bo_fixed = bo_fixed,
        max_refinements = BO_MAX_REFINEMENTS,
        tol_gain = BO_TOL_GAIN,
        plateau = BO_PLATEAU,
        top_k = BO_TOP_K,
        shrink_base = BO_SHRINK_BASE,
        importance_gain = BO_IMPORTANCE_GAIN,
        coarse_control = BO_COARSE_CONTROL,
        refine_control = BO_REFINE_CONTROL,
        verbose = TRUE,
        parallel_grid = TRUE,
        ncores_grid = NINIT
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  
  tar_target(
    params_best,
    standardize_bo_params(desurv_bo_results$overall_best$params)
  ),
  
  tar_target(
    ngene_value,
    {
      value <- params_best$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(NGENE_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    ntop_value,
    {
      value <- params_best$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(NTOP_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    lambdaW_value,
    {
      value <- params_best$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAW_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    lambdaH_value,
    {
      value <- params_best$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAH_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),
  
  tar_target(
    training_results_dir,
    results_root_dir(
      ngene = ngene_value,
      tol = TOL,
      maxit = MAXIT,
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = TRAIN_PREFIX,
      method_trans_train = METHOD_TRANS_TRAIN
    )
  ),
  
  tar_target(
    desurv_bo_history,
    {
      path <- file.path(training_results_dir, "desurv_bo_history.csv")
      utils::write.csv(desurv_bo_results$history, path, row.names = FALSE)
      path
    },
    format = "file"
  ),
  
  tar_target(
    data_filtered,
    {
      prep <- DeSurv::preprocess_data(
        X = data$ex,
        y = data$sampInfo$time,
        d = data$sampInfo$event,
        dataset = data$sampInfo$dataset,
        samp_keeps = data$samp_keeps,
        ngene = ngene_value,
        method_trans_train = METHOD_TRANS_TRAIN,
        verbose = FALSE
      )
      prep$dataname <- data$dataname
      prep
    }
  ),

  tar_target(
    desurv_seed_fits,
    {
      seeds <- seq_len(NINIT_FULL)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fit_i <- try(
          desurv_fit(
            X = data_filtered$ex,
            y = data_filtered$sampInfo$time,
            d = data_filtered$sampInfo$event,
            k = params_best$k,
            alpha = params_best$alpha,
            lambda = params_best$lambda,
            nu = params_best$nu,
            lambdaW = lambdaW_value,
            lambdaH = lambdaH_value,
            seed = seeds[i],
            tol = TOL / 100,
            tol_init = TOL,
            maxit = MAXIT,
            imaxit = MAXIT,
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
        X = data_filtered$ex,
        ntop = ntop_value,
        k = params_best$k,
        min_frequency = 0.3*NINIT_FULL
      )
    }
  ),
  
  tar_target(
    fit_desurv,
    {
      init_vals <- desurv_consensus_init
      desurv_fit(
        X = data_filtered$ex,
        y = data_filtered$sampInfo$time,
        d = data_filtered$sampInfo$event,
        k = params_best$k,
        alpha = params_best$alpha,
        lambda = params_best$lambda,
        nu = params_best$nu,
        lambdaW = lambdaW_value,
        lambdaH = lambdaH_value,
        W0 = init_vals$W0,
        H0 = init_vals$H0,
        beta0 = init_vals$beta0,
        seed = NULL,
        tol = TOL / 100,
        maxit = MAXIT,
        verbose = FALSE
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),

  tar_target(
    desurv_bo_results_alpha0,
    {
      bounds <- DESURV_BO_BOUNDS
      bounds$alpha_grid <- NULL
      bounds <- maybe_add_numeric_bound(bounds, NGENE_CONFIG, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, NTOP_CONFIG, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAW_CONFIG, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAH_CONFIG, "lambdaH_grid", log_scale = TRUE)

      bo_fixed <- list(
        n_starts = NINIT,
        alpha_grid = 0
      )
      if (!TUNE_NGENE) {
        bo_fixed$ngene <- NGENE_DEFAULT
      }
      if (!TUNE_NTOP) {
        bo_fixed$ntop <- NTOP_DEFAULT
      }
      if (!TUNE_LAMBDAW) {
        bo_fixed$lambdaW_grid <- LAMBDAW_DEFAULT
      }
      if (!TUNE_LAMBDAH) {
        bo_fixed$lambdaH_grid <- LAMBDAH_DEFAULT
      }

      DeSurv::desurv_cv_bayesopt_refine(
        X = data$ex,
        y = data$sampInfo$time,
        d = data$sampInfo$event,
        dataset = data$sampInfo$dataset,
        samp_keeps = data$samp_keeps,
        preprocess = TRUE,
        method_trans_train = METHOD_TRANS_TRAIN,
        engine = "warmstart",
        nfolds = NFOLD,
        tol = TOL,
        maxit = MAXIT,
        coarse_bounds = bounds,
        bo_fixed = bo_fixed,
        max_refinements = BO_MAX_REFINEMENTS,
        tol_gain = BO_TOL_GAIN,
        plateau = BO_PLATEAU,
        top_k = BO_TOP_K,
        shrink_base = BO_SHRINK_BASE,
        importance_gain = BO_IMPORTANCE_GAIN,
        coarse_control = BO_COARSE_CONTROL,
        refine_control = BO_REFINE_CONTROL,
        verbose = TRUE,
        parallel_grid = TRUE,
        ncores_grid = NINIT
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),

  tar_target(
    params_best_alpha0,
    {
      params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
      params$alpha <- 0
      params
    }
  ),

  tar_target(
    ngene_value_alpha0,
    {
      value <- params_best_alpha0$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(NGENE_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    ntop_value_alpha0,
    {
      value <- params_best_alpha0$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(NTOP_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    lambdaW_value_alpha0,
    {
      value <- params_best_alpha0$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAW_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    lambdaH_value_alpha0,
    {
      value <- params_best_alpha0$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAH_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    training_results_dir_alpha0,
    results_root_dir(
      ngene = ngene_value_alpha0,
      tol = TOL,
      maxit = MAXIT,
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = paste0(TRAIN_PREFIX, "_alpha0"),
      method_trans_train = METHOD_TRANS_TRAIN
    )
  ),

  tar_target(
    desurv_bo_history_alpha0,
    {
      path <- file.path(training_results_dir_alpha0, "desurv_bo_history_alpha0.csv")
      utils::write.csv(desurv_bo_results_alpha0$history, path, row.names = FALSE)
      path
    },
    format = "file"
  ),

  tar_target(
    data_filtered_alpha0,
    {
      prep <- DeSurv::preprocess_data(
        X = data$ex,
        y = data$sampInfo$time,
        d = data$sampInfo$event,
        dataset = data$sampInfo$dataset,
        samp_keeps = data$samp_keeps,
        ngene = ngene_value_alpha0,
        method_trans_train = METHOD_TRANS_TRAIN,
        verbose = FALSE
      )
      prep$dataname <- data$dataname
      prep
    }
  ),

  tar_target(
    desurv_seed_fits_alpha0,
    {
      seeds <- seq_len(NINIT_FULL)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fit_i <- try(
          desurv_fit(
            X = data_filtered_alpha0$ex,
            y = data_filtered_alpha0$sampInfo$time,
            d = data_filtered_alpha0$sampInfo$event,
            k = params_best_alpha0$k,
            alpha = params_best_alpha0$alpha,
            lambda = params_best_alpha0$lambda,
            nu = params_best_alpha0$nu,
            lambdaW = lambdaW_value_alpha0,
            lambdaH = lambdaH_value_alpha0,
            seed = seeds[i],
            tol = TOL / 100,
            tol_init = TOL,
            maxit = MAXIT,
            imaxit = MAXIT,
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
        X = data_filtered_alpha0$ex,
        ntop = ntop_value_alpha0,
        k = params_best_alpha0$k,
        min_frequency = 0.3 * NINIT_FULL
      )
    }
  ),

  tar_target(
    fit_desurv_alpha0,
    {
      init_vals <- desurv_consensus_init_alpha0
      desurv_fit(
        X = data_filtered_alpha0$ex,
        y = data_filtered_alpha0$sampInfo$time,
        d = data_filtered_alpha0$sampInfo$event,
        k = params_best_alpha0$k,
        alpha = params_best_alpha0$alpha,
        lambda = params_best_alpha0$lambda,
        nu = params_best_alpha0$nu,
        lambdaW = lambdaW_value_alpha0,
        lambdaH = lambdaH_value_alpha0,
        W0 = init_vals$W0,
        H0 = init_vals$H0,
        beta0 = init_vals$beta0,
        seed = NULL,
        tol = TOL / 100,
        maxit = MAXIT,
        verbose = FALSE
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),

  tar_target(
    fit_std,
    nmf(data_filtered$ex,STD_NMF_K_GRID,nrun=NINIT,method="lee",.options=paste0("p",NINIT)),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  
  tar_target(
    std_nmf_k_selection_plots,
    {
      save_dir=file.path(training_results_dir,"std_nmf_k_selection")
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
  #     X <- data_filtered$ex
  #     y <- data_filtered$sampInfo$time
  #     d <- data_filtered$sampInfo$event
  #     strata <- interaction(d, data_filtered$sampInfo$dataset, drop = FALSE)
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
  tar_target(
    data_val_filtered,
    {
      datasets_named <- data_val
      if (is.null(names(datasets_named)) || any(!nzchar(names(datasets_named)))) {
        names(datasets_named) <- vapply(
          datasets_named,
          function(dataset) {
            dataname <- dataset$dataname
            if (!is.null(dataname) && nzchar(dataname)) {
              dataname
            } else {
              dataset_ids <- unique(dataset$sampInfo$dataset)
              dataset_ids <- dataset_ids[!is.na(dataset_ids)]
              if (length(dataset_ids)) dataset_ids[[1]] else "validation"
            }
          },
          character(1)
        )
      }
      setNames(
        lapply(
          seq_along(datasets_named),
          function(idx) {
            dataset <- datasets_named[[idx]]
            dataname <- names(datasets_named)[idx]
            gene_arg <- if (USE_TRAIN_GENES_FOR_VAL) rownames(data_filtered$ex) else NULL
            prep <- DeSurv::preprocess_data(
              X = dataset$ex,
              y = dataset$sampInfo$time,
              d = dataset$sampInfo$event,
              dataset = dataset$sampInfo$dataset,
              samp_keeps = dataset$samp_keeps,
              genes = gene_arg,
              ngene = ngene_value,
              method_trans_train = METHOD_TRANS_TRAIN,
              verbose = FALSE
            )
            prep$dataname <- dataname
            prep
          }
        ),
        names(datasets_named)
      )
    }
  ),
  tar_target(
    val_predictions_desurv,
    desurv_predict_validation(
      fit = fit_desurv,
      data_list = data_val_filtered,
      top_genes = tops_desurv$top_genes
    )
  ),
  tar_target(
    val_predictions_desurv_alpha0,
    desurv_predict_validation(
      fit = fit_desurv_alpha0,
      data_list = data_val_filtered,
      top_genes = tops_desurv_alpha0$top_genes
    )
  ),
  tar_target(
    val_latent_desurv,
    {
      latent <- desurv_collect_validation_latent(
        fit = fit_desurv,
        data_list = data_val_filtered,
        top_genes = tops_desurv$top_genes
      )
      write_validation_latent_outputs(
        latent_list = latent,
        base_dir = file.path(training_results_dir, "validation", "desurv")
      )
      latent
    }
  ),
  tar_target(
    val_latent_desurv_alpha0,
    {
      latent <- desurv_collect_validation_latent(
        fit = fit_desurv_alpha0,
        data_list = data_val_filtered,
        top_genes = tops_desurv_alpha0$top_genes
      )
      write_validation_latent_outputs(
        latent_list = latent,
        base_dir = file.path(training_results_dir_alpha0, "validation", "desurv_alpha0")
      )
      latent
    }
  ),
  tar_target(
    val_cindex_desurv,
    {
      summary_tbl <- summarize_validation_cindex(val_latent_desurv)
      if (nrow(summary_tbl)) {
        dir.create(file.path(training_results_dir, "validation"), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(
          summary_tbl,
          file = file.path(training_results_dir, "validation", "val_cindex_desurv.csv"),
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
        dir.create(file.path(training_results_dir_alpha0, "validation"), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(
          summary_tbl,
          file = file.path(training_results_dir_alpha0, "validation", "val_cindex_desurv_alpha0.csv"),
          row.names = FALSE
        )
      }
      summary_tbl
    }
  ),
  tar_target(
    val_clusters_desurv,
    {
      base_dir <- file.path(training_results_dir, "validation", "desurv", "clusters")
      dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
      lapply(
        val_latent_desurv,
        function(entry) {
          dataset_dir <- file.path(base_dir, entry$dataset)
          dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
          result <- cluster_validation_latent(
            latent_entry = entry,
            maxK = VAL_CLUSTER_MAXK,
            reps = VAL_CLUSTER_REPS,
            pItem = VAL_CLUSTER_PITEM,
            pFeature = VAL_CLUSTER_PFEATURE,
            seed = VAL_CLUSTER_SEED,
            dir = dataset_dir
          )
          if (nrow(result$assignments)) {
            utils::write.csv(
              result$assignments,
              file = file.path(dataset_dir, "cluster_assignments.csv"),
              row.names = FALSE
            )
          }
          result
        }
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),
  tar_target(
    val_clusters_desurv_alpha0,
    {
      base_dir <- file.path(training_results_dir_alpha0, "validation", "desurv_alpha0", "clusters")
      dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
      lapply(
        val_latent_desurv_alpha0,
        function(entry) {
          dataset_dir <- file.path(base_dir, entry$dataset)
          dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
          result <- cluster_validation_latent(
            latent_entry = entry,
            maxK = VAL_CLUSTER_MAXK,
            reps = VAL_CLUSTER_REPS,
            pItem = VAL_CLUSTER_PITEM,
            pFeature = VAL_CLUSTER_PFEATURE,
            seed = VAL_CLUSTER_SEED,
            dir = dataset_dir
          )
          if (nrow(result$assignments)) {
            utils::write.csv(
              result$assignments,
              file = file.path(dataset_dir, "cluster_assignments.csv"),
              row.names = FALSE
            )
          }
          result
        }
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),
  tar_target(
    tops_desurv,
    get_top_genes(W = fit_desurv$W, ntop = ntop_value)
  ),
  tar_target(
    gene_overlap_desurv,
    {
      create_table(tops = tops_desurv$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),

  tar_target(
    tops_desurv_alpha0,
    get_top_genes(W = fit_desurv_alpha0$W, ntop = ntop_value_alpha0)
  ),
  tar_target(
    gene_overlap_desurv_alpha0,
    {
      create_table(tops = tops_desurv_alpha0$top_genes, gene_lists = top_genes,
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
      universe = rownames(data_filtered$ex)
      organism <- org.Hs.eg.db
      ora(tops_desurv$top_genes,universe,organism)
    }
  ),

  tar_target(
    ora_analysis_desurv_alpha0,
    {
      universe = rownames(data_filtered_alpha0$ex)
      organism <- org.Hs.eg.db
      ora(tops_desurv_alpha0$top_genes,universe,organism)
    }
  )
  
  # tar_target(
  #   ora_analysis_std,
  #   {
  #     universe = rownames(data_filtered$ex)
  #     organism <- org.Hs.eg.db
  #     ora(tops_std$top_genes,universe,organism)
  #   }
  # ),
  
  # 
  # 
  # tar_target(
  #   selected_factors_desurv,
  #   {
  #     save_dir = file.path(training_results_dir,"factor_selection","desurv")
  #     dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
  #     path = file.path(save_dir,"factor_selection.csv")
  #     build_factor_selection_table(tops_desurv$top_genes,path)
  #   },
  #   format="file"
  # ),

  # tar_target(
  #   selected_factors_desurv_alpha0,
  #   {
  #     save_dir = file.path(training_results_dir_alpha0,"factor_selection","desurv_alpha0")
  #     dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
  #     path = file.path(save_dir,"factor_selection.csv")
  #     build_factor_selection_table(tops_desurv_alpha0$top_genes,path)
  #   },
  #   format="file"
  # ),
  # 
  # tar_target(
  #   selected_factors_std,
  #   {
  #     save_dir = file.path(training_results_dir,"factor_selection","std")
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
    #     dir = create_filepath_clustering_output(ngene_value, TOL, MAXIT, PKG_VERSION,
    #                                             GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN,
    #                                             ntop_value, val_dataset, "DeSurv")
    #     run_clustering(
    #       tops_desurv$top_genes, data, top_genes, colors,
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
    #     dir = create_filepath_clustering_output(ngene_value_alpha0, TOL, MAXIT, PKG_VERSION,
    #                                             GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN,
    #                                             ntop_value_alpha0, val_dataset, "DeSurv_alpha0")
    #     run_clustering(
    #       tops_desurv_alpha0$top_genes, data, top_genes, colors,
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
  #       dir = create_filepath_clustering_output(ngene_value, TOL, MAXIT, PKG_VERSION, 
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
  #     save_dir = file.path(training_results_dir,"ncluster_selection","desurv")
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
  #     save_dir = file.path(training_results_dir,"ncluster_selection","std")
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
  #     save_dir = file.path(training_results_dir,"cluster_alignment","desurv")
  #     dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
  #     path = file.path(save_dir,"cluster_alignment.pdf")
  #     pdf(path)
  #     plot_cluster_cor(clusters_desurv,tops_desurv$top_genes,nclus_vec)
  #     dev.off()
  #   }
  # ),
  # 
  # tar_target(
  #   cluster_alignment_plot_std,
  #   {
  #     nclus_tbl = read.csv(selected_nclusters_std)
  #     nclus_vec <- setNames(nclus_tbl$nclus, nclus_tbl$dataset)
  #     save_dir = file.path(training_results_dir,"cluster_alignment","std")
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
  #     save_dir = file.path(training_results_dir,"cluster_alignment","desurv")
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
  #     save_dir = file.path(training_results_dir,"cluster_alignment","std")
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
