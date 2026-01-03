source("targets_configs.R")
BO_CONFIGS_RAW <- targets_bo_configs()
RUN_CONFIGS_RAW <- targets_run_configs()
VAL_CONFIGS_RAW <- targets_val_configs()
DEFAULT_NINIT <- if (length(BO_CONFIGS_RAW)) {
  max(vapply(BO_CONFIGS_RAW, function(cfg) if (is.null(cfg$ninit)) 50 else cfg$ninit, numeric(1)))
} else {
  50
}
DEFAULT_NINIT_FULL <- if (length(RUN_CONFIGS_RAW)) {
  max(vapply(RUN_CONFIGS_RAW, function(cfg) if (is.null(cfg$ninit_full)) 100 else cfg$ninit_full, numeric(1)))
} else {
  100
}

source("targets_setup.R")
source("targets_common_pipeline.R")

RESOLVED_BO_CONFIGS <- resolve_desurv_bo_configs(BO_CONFIGS_RAW)
RESOLVED_RUN_CONFIGS <- resolve_desurv_run_configs(RUN_CONFIGS_RAW)
RESOLVED_VAL_CONFIGS <- resolve_desurv_val_configs(VAL_CONFIGS_RAW)
validate_desurv_configs(RESOLVED_BO_CONFIGS, RESOLVED_RUN_CONFIGS, RESOLVED_VAL_CONFIGS)


# ---- Targets ----
targets_list <- list(
  tar_target(
    bo_config,
    RESOLVED_BO_CONFIGS,
    iteration = "list"
  ),
  tar_target(
    run_config,
    RESOLVED_RUN_CONFIGS,
    iteration = "list"
  ),
  tar_target(
    val_config,
    RESOLVED_VAL_CONFIGS,
    iteration = "list"
  ),

  ################# Training ###################

  tar_target(
    raw_data,
    {
      if (bo_config$data_mode == "external") {
        train_datasets <- bo_config$train_datasets
        c(
          file.path("data/original", paste0(train_datasets, ".rds")),
          file.path("data/original", paste0(train_datasets, ".survival_data.rds")),
          file.path("data/original", paste0(train_datasets, "_subtype.csv"))
        )
      } else {
        bo_config$split_raw_files
      }
    },
    format = "file"
  ),
  tar_target(
    data_input,
    {
      raw_data
      if (bo_config$data_mode == "split") {
        loader <- get(bo_config$data_loader, mode = "function")
        loader(raw_data)
      } else {
        NULL
      }
    }
  ),
  tar_target(
    data_split,
    {
      if (bo_config$data_mode == "split") {
        split_train_validation(
          data = data_input,
          train_frac = bo_config$split_train_frac,
          seed = bo_config$split_seed,
          strata_vars = bo_config$split_strata_vars
        )
      } else {
        NULL
      }
    }
  ),
  tar_target(
    data,
    { 
      raw_data
      if (bo_config$data_mode == "split") {
        data_split$train
      } else {
        loader <- get(bo_config$data_loader, mode = "function")
        loader(bo_config$train_datasets)
      }
    }
  ),
  tar_target(
    val_datasets,
    val_config$val_datasets
  ),
  tar_target(
    val_dataset_name,
    val_datasets,
    iteration = "vector"
  ),
  
  tar_target(
    raw_data_val,
    if (val_config$mode == "external") {
      c(
        file.path("data/original", paste0(val_dataset_name, ".rds")),
        file.path("data/original", paste0(val_dataset_name, ".survival_data.rds")),
        file.path("data/original", paste0(val_dataset_name, "_subtype.csv"))
      )
    } else {
      character(0)
    },
    format = "file",
    pattern = map(val_dataset_name)
  ),
  tar_target(
    data_val,
    {
      if (val_config$mode == "train_split") {
        dataset <- val_run_bundle$bo_bundle$data_split$test
        dataname <- infer_dataset_name(dataset)
        setNames(list(dataset), dataname)
      } else {
        setNames(
          lapply(val_datasets, load_data),
          val_datasets
        )
      }
    }
  ),
  
  tar_target(
    data_val_comb,
    {
      if (val_config$mode == "external") {
        raw_data_val
        load_data(val_datasets)
      } else {
        NULL
      }
    },
  ),
  
  tar_target(
    data_val_comb_filtered,
    {
      if (val_config$mode == "external") {
        genes_train = rownames(val_run_bundle$bo_bundle$data_filtered$ex)
        preprocess_validation_data(
          dataset = data_val_comb,
          genes = genes_train,
          method_trans_train = val_run_bundle$bo_bundle$config$method_trans_train,
          dataname = data_val_comb$dataname
        )
      } else {
        NULL
      }
    },
  )
)

c(
  targets_list,
  COMMON_DESURV_BO_TARGETS,
  COMMON_DESURV_RUN_TARGETS#,
  # tarchetypes::tar_render(
  #   paper,
  #   "paper/paper.Rmd",
  #   quiet = FALSE,
  #   deps = c(
  #     "desurv_bo_history",
  #     "desurv_bo_history_alpha0",
  #     "params_best",
  #     "ngene_value",
  #     "ntop_value",
  #     "fit_desurv",
  #     "fit_desurv_alpha0",
  #     "fit_std",
  #     "fit_std_beta",
  #     "tops_desurv",
  #     "tops_desurv_alpha0",
  #     "tops_std",
  #     "gene_overlap_desurv",
  #     "gene_overlap_desurv_alpha0",
  #     "gene_overlap_std",
  #     "ora_analysis_desurv",
  #     "ora_analysis_desurv_alpha0",
  #     "ora_analysis_std",
  #     "selected_factors_desurv",
  #     "selected_factors_desurv_alpha0",
  #     "selected_factors_std",
  #     "clusters_desurv",
  #     "clusters_std",
  #     "aligned_clusters_desurv",
  #     "aligned_clusters_std",
  #     "cluster_alignment_plot_desurv",
  #     "cluster_alignment_plot_std",
  #     "cluster_alignment_table_desurv",
  #     "cluster_alignment_table_std",
  #     "data_filtered",
  #     "data_val",
  #     "data_val_filtered",
  #     "training_results_dir"
  #   )
  # )
)

  # 
  # tar_render(
  #   paper,
  #   "paper/paper.Rmd"
  # )
  
  
  # 
  #   tar_target(
  #     best_init_per_param_combo,
  #     {
  #       df = readRDS(inits)
  #       select_best_init(df = df, method_select = METHOD_SELECT_INIT)
  #     },
  #     pattern   = map(inits),
  #     iteration = "list"
  #   ),
  #   
  #   tar_target(
  #     model_runs,
  #     {
  #       print("running coldstarts...")
  # 
  #       path = create_filepath_coldstart_runs(params = best_init_per_param_combo)
  # 
  #       run_coldstarts(
  #         X = data_filtered$ex, y = data_filtered$sampInfo$time, delta = data_filtered$sampInfo$event,
  #         params = best_init_per_param_combo,
  #         verbose = FALSE,
  #         path = path
  #       )
  # 
  #     },
  #     pattern   = map(best_init_per_param_combo),
  #     iteration = "list",
  #     format = "file",
  #     resources = tar_resources(
  #       crew = tar_resources_crew(controller = "model_runs")
  #     )
  #   ),
  #   # #
  #   # Compute model metrics
  #   tar_target(
  #     metrics,
  #     compute_metrics(model_runs,
  #                     data_filtered$ex,
  #                     data_filtered$sampInfo$time,
  #                     data_filtered$sampInfo$event),
  #     pattern   = map(model_runs),
  #     iteration = "list"
  #   ),
  #   # 
  #   tar_target(
  #     metrics_table,
  #     dplyr::bind_rows(metrics)
  #   )
  # 
  
  
  # run and save initializations for each parameter combo
  
  # 
  # # read in all initializations
  # tar_target(
  #   alpha0_scores,
  #   {
  #     df = dplyr::bind_rows(lapply(alpha0_score_file, function(path) {
  #       readRDS(path)
  #     }))
  #     df
  #   }
  # ),
  # 
  # tarchetypes::tar_group_by(
  #   alpha0_groups,                 # = grouped target name
  #   alpha0_scores,                 # = your tibble from reading *.rds
  #   groups = c(k, lambda, eta, lambdaW, lambdaH)
  # ),
  # 
  
  
  
  #
  # # initialize full model run for selected params
  
  
  ##### Summaries #####
