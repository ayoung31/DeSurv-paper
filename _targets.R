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
targets_list <- tar_map(
  values = data.frame(
    bo_config_value = I(RESOLVED_BO_CONFIGS),
    bo_label = names(RESOLVED_BO_CONFIGS),
    stringsAsFactors = FALSE
  ),
  names = "bo_label",
  tar_target(
    bo_config,
    bo_config_value
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
    tar_data_split,
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
    tar_data,
    { 
      raw_data
      if (bo_config$data_mode == "split") {
        tar_data_split$train
      } else {
        loader <- get(bo_config$data_loader, mode = "function")
        loader(bo_config$train_datasets)
      }
    }
  ),
  COMMON_DESURV_BO_TARGETS,
  tar_map(
    values = data.frame(
      run_config_value = I(RESOLVED_RUN_CONFIGS),
      run_label = names(RESOLVED_RUN_CONFIGS),
      stringsAsFactors = FALSE
    ),
    names = "run_label",
    tar_target(
      run_config,
      run_config_value
    ),
    COMMON_DESURV_RUN_TARGETS,
    tar_map(
          values = data.frame(
            val_config_value = I(unname(RESOLVED_VAL_CONFIGS)),
            val_label = names(RESOLVED_VAL_CONFIGS),
            stringsAsFactors = FALSE
          ),
      names = "val_label",
      tar_target(
        val_config,
        val_config_value
      ),
      tar_target(
        val_config_effective,
        {
          cfg <- val_config
          if (identical(bo_config$data_mode, "split")) {
            cfg$mode <- "train_split"
            cfg$val_datasets <- character(0)
          } else if (identical(bo_config$data_mode, "external")) {
            cfg$mode <- "external"
            val_datasets <- cfg$val_datasets
            if (is.null(val_datasets)) {
              val_datasets <- character(0)
            }
            if (!length(val_datasets)) {
              stop(sprintf(
                "val_config '%s' must define val_datasets for external validation.",
                cfg$label
              ))
            }
            overlap <- intersect(val_datasets, bo_config$train_datasets)
            if (length(overlap)) {
              val_datasets <- setdiff(val_datasets, overlap)
            }
            if (!length(val_datasets)) {
              stop(sprintf(
                "val_config '%s' has no external datasets after removing training data.",
                cfg$label
              ))
            }
            cfg$val_datasets <- val_datasets
          } else {
            stop(sprintf(
              "bo_config '%s' has unsupported data_mode '%s'.",
              bo_config$label,
              bo_config$data_mode
            ))
          }

          cfg$use_train_genes_for_val <- NULL
          cfg$val_cluster_maxk <- NULL
          cfg$val_cluster_reps <- NULL
          cfg$val_cluster_pitem <- NULL
          cfg$val_cluster_pfeature <- NULL
          cfg$val_cluster_seed <- NULL

          cfg$config_id <- desurv_val_config_hash(cfg)
          cfg$short_id <- substr(cfg$config_id, 1, 8)
          cfg$path_tag <- build_config_tag(cfg$label, cfg$config_id)
          cfg
        }
      ),
      tar_target(
        val_datasets_raw,
        val_config_effective$val_datasets
      ),
  tar_target(
    tar_val_datasets,
    if (length(val_datasets_raw)) val_datasets_raw else NA_character_
  ),
  tar_target(
    val_dataset_name,
    tar_val_datasets,
    iteration = "vector"
  ),
      tar_target(
        raw_data_val,
        if (val_config_effective$mode == "external") {
          if (is.na(val_dataset_name) || !nzchar(val_dataset_name)) {
            character(0)
          } else {
            c(
              file.path("data/original", paste0(val_dataset_name, ".rds")),
              file.path("data/original", paste0(val_dataset_name, ".survival_data.rds")),
              file.path("data/original", paste0(val_dataset_name, "_subtype.csv"))
            )
          }
        } else {
          character(0)
        },
        format = "file",
        pattern = map(val_dataset_name)
      ),
      tar_target(
        data_val,
        {
          if (val_config_effective$mode == "train_split") {
            dataset <- val_run_bundle$bo_bundle$data_split$test
            dataname <- infer_dataset_name(dataset)
            setNames(list(dataset), dataname)
          } else {
            if (length(val_datasets_raw)) {
              # Reference raw_data_val to establish dependency so file changes invalidate this target
              raw_data_val
              setNames(
                lapply(val_datasets_raw, load_data),
                val_datasets_raw
              )
            } else {
              list()
            }
          }
        }
      ),
      tar_target(
        data_val_comb,
        {
          if (val_config_effective$mode == "external") {
            if (length(val_datasets_raw)) {
              raw_data_val
              load_data(val_datasets_raw)
            } else {
              NULL
            }
          } else {
            NULL
          }
        },
      ),
      tar_target(
        data_val_comb_filtered,
        {
          if (val_config_effective$mode == "external") {
            genes_train = rownames(val_run_bundle$bo_bundle$data_filtered$ex)
            preprocess_validation_data(
              dataset = data_val_comb,
              genes = genes_train,
              method_trans_train = val_run_bundle$bo_bundle$config$method_trans_train,
              dataname = data_val_comb$dataname,
              transform_target = val_run_bundle$bo_bundle$data_filtered$transform_target
            )
          } else {
            NULL
          }
        },
      ),
      COMMON_DESURV_VAL_TARGETS
    )
  )
)

c(
  # tar_target(
  #   version_info,
  #   {
  #     list(
  #       desurv = list(
  #         branch = gert::git_branch(repo = "../DeSurv"),
  #         commit = gert::git_log(repo = "../DeSurv", max = 1)$commit,
  #         date = Sys.time()
  #       ),
  #       paper = list(
  #         branch = gert::git_branch(),
  #         commit = gert::git_log(max = 1)$commit
  #       )
  #     )
  #   }
  # ),
  targets_list
  # COMMON_DESURV_BO_TARGETS,
  # COMMON_DESURV_RUN_TARGETS,
  # COMMON_DESURV_VAL_TARGETS#,
  # tarchetypes::tar_render(
  #   paper,
  #   "paper/paper.Rmd",
  #   quiet = FALSE,
  #   deps = c(
  #     "desurv_bo_history",
  #     "desurv_bo_history_alpha0",
  #     "tar_params_best",
  #     "tar_ngene_value",
  #     "ntop_value",
  #     "tar_fit_desurv",
  #     "tar_fit_desurv_alpha0",
  #     "fit_std",
  #     "fit_std_beta",
  #     "tar_tops_desurv",
  #     "tar_tops_desurv_alpha0",
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
  #     "tar_data_filtered",
  #     "data_val",
  #     "data_val_filtered",
  #     "tar_training_results_dir"
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
  #         X = tar_data_filtered$ex, y = tar_data_filtered$sampInfo$time, delta = tar_data_filtered$sampInfo$event,
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
  #                     tar_data_filtered$ex,
  #                     tar_data_filtered$sampInfo$time,
  #                     tar_data_filtered$sampInfo$event),
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
