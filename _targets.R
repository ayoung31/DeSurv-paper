NINIT <- 50
NINIT_FULL <- 100
BO_N_INIT <- 20
BO_N_ITER <- 100
BO_CANDIDATE_POOL <- 2000
BO_MAX_REFINEMENTS <- 3
BO_TOL_GAIN <- 0.002
BO_PLATEAU <- 1
BO_TOP_K <- 10
BO_SHRINK_BASE <- 0.5
BO_IMPORTANCE_GAIN <- 0.3
BO_COARSE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 123,
  cv_verbose = FALSE
)
BO_REFINE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 456,
  cv_verbose = FALSE
)

VAL_DATASETS       = c("Dijk","Moffitt_GEO_array",
                       "PACA_AU_array","PACA_AU_seq","Puleo_array")
TRAIN_DATASETS     = c("TCGA_PAAD","CPTAC")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
USE_TRAIN_GENES_FOR_VAL <- FALSE

##### hyperparameters #####

# always tuned hyperparams
DESURV_BO_BOUNDS   = list(
  k_grid = list(lower = 2L, upper = 10L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

# sometimes tuned hyperparams
# for the following set length 1 to fix, or range to tune via BO
NGENE_CONFIG       = c(500,5000) 
NTOP_CONFIG        = c(50,250) 
LAMBDAW_CONFIG     = c(0)
LAMBDAH_CONFIG     = c(1e-7,1e2)


source("targets_setup.R")
source("targets_common_pipeline.R")


# ---- Targets ----
targets_list <- list(

  ################# Training ###################

  tar_target(
    raw_data,
    c(
      file.path("data/original", paste0(TRAIN_DATASETS, ".rds")),
      file.path("data/original", paste0(TRAIN_DATASETS, ".survival_data.rds")),
      file.path("data/original", paste0(TRAIN_DATASETS, "_subtype.csv"))
    ),
    format = "file"
  ),
  tar_target(
    data,
    { 
      raw_data
      load_data(TRAIN_DATASETS) 
    }
  ),
  tar_target(
    val_datasets,
    VAL_DATASETS
  ),
  tar_target(
    val_dataset_name,
    val_datasets,
    iteration = "vector"
  ),
  
  tar_target(
    raw_data_val,
    c(
      file.path("data/original", paste0(val_dataset_name, ".rds")),
      file.path("data/original", paste0(val_dataset_name, ".survival_data.rds")),
      file.path("data/original", paste0(val_dataset_name, "_subtype.csv"))
    ),
    format = "file",
    pattern = map(val_dataset_name)
  ),
  tar_target(
    data_val,
    setNames(
      lapply(val_datasets, load_data),
      val_datasets
    )
  ),
  
  tar_target(
    data_val_comb,
    {
      raw_data_val
      load_data(val_datasets)
    },
  ),
  
  tar_target(
    data_val_comb_filtered,
    {
      genes_train = rownames(data_filtered$ex)
      prep <- DeSurv::preprocess_data(
        X = data_val_comb$ex,
        y = data_val_comb$sampInfo$time,
        d = data_val_comb$sampInfo$event,
        dataset = data_val_comb$sampInfo$dataset,
        samp_keeps = data_val_comb$samp_keeps,
        genes = genes_train,
        method_trans_train = METHOD_TRANS_TRAIN,
        verbose = FALSE
      )
      prep$dataname <- data_val_comb$dataname
      prep
    },
  )
)

c(
  targets_list,
  COMMON_DESURV_TARGETS,
  tarchetypes::tar_render(
    paper,
    "paper/paper.Rmd",
    quiet = FALSE,
    deps = c(
      "desurv_bo_history",
      "desurv_bo_history_alpha0",
      "params_best",
      "ngene_value",
      "ntop_value",
      "fit_desurv",
      "fit_desurv_alpha0",
      "fit_std",
      "fit_std_beta",
      "tops_desurv",
      "tops_desurv_alpha0",
      "tops_std",
      "gene_overlap_desurv",
      "gene_overlap_desurv_alpha0",
      "gene_overlap_std",
      "ora_analysis_desurv",
      "ora_analysis_desurv_alpha0",
      "ora_analysis_std",
      "selected_factors_desurv",
      "selected_factors_desurv_alpha0",
      "selected_factors_std",
      "clusters_desurv",
      "clusters_std",
      "aligned_clusters_desurv",
      "aligned_clusters_std",
      "cluster_alignment_plot_desurv",
      "cluster_alignment_plot_std",
      "cluster_alignment_table_desurv",
      "cluster_alignment_table_std",
      "data_filtered",
      "data_val",
      "data_val_filtered",
      "training_results_dir"
    )
  )
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
