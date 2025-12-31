
# ---- Training parameters ----
TRAIN_DATASETS     = c("Bladder")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
TRAIN_PERC         = .7
METHOD_TRANS_TRAIN = "none"
USE_TRAIN_GENES_FOR_VAL <- TRUE

# always tuned hyper params
DESURV_BO_BOUNDS   = list(
  k_grid = list(lower = 2L, upper = 10L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

# sometimes tuned hyperparams
# for the following set length 1 to fix, or range to tune via BO
NGENE_CONFIG       = c(500,5000) 
NTOP_CONFIG        = c(50,250) 
LAMBDAW_CONFIG     = c(0)
LAMBDAH_CONFIG     = 0#c(1e-5,1e5)
NINIT <- 50
NINIT_FULL <- 100
BO_N_INIT <- 15
BO_N_ITER <- 60
BO_CANDIDATE_POOL <- 4000
BO_MAX_REFINEMENTS <- 2
BO_TOL_GAIN <- 0.002
BO_PLATEAU <- 1
BO_TOP_K <- 10
BO_SHRINK_BASE <- 0.3
BO_IMPORTANCE_GAIN <- 0.1
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

source("targets_setup.R")
source("targets_common_pipeline.R")

# ---- Targets ----
targets_list <- list(

  ################# Training ###################

  tar_target(
    raw_data,
      "data/original/IMVigor210_sampInfo_ex_clinical.rds",
    format = "file"
  ),
  tar_target(
    data_input,
    { 
      load_data_bladder_vig(raw_data)
    }
  ),
  
  tar_target(
    data_split,
    {
      data = data_input
      samp_keeps = data$samp_keeps
      data$sampInfo = data$sampInfo[samp_keeps,]
      data$ex = data$ex[,samp_keeps]
      strata = interaction(data$sampInfo$event,data$sampInfo$dataset,drop = FALSE)
      idx_train = caret::createDataPartition(strata,p=TRAIN_PERC)[[1]]
      idx_test = setdiff(1:length(strata),idx_train)
      
      data_train = data
      data_train$sampInfo = data_train$sampInfo[idx_train,]
      data_train$ex = data_train$ex[,idx_train]
      data_train$samp_keeps = 1:nrow(data_train$sampInfo)
      
      data_test = data
      data_test$sampInfo = data_test$sampInfo[idx_test,]
      data_test$ex = data_test$ex[,idx_test]
      data_test$samp_keeps = 1:nrow(data_test$sampInfo)
      
      list(train=data_train,test=data_test)
    }
  ),
  
  tar_target(
    data,
    data_split$train
  ),
  
  tar_target(
    data_val,
    {
      dataset <- data_split$test
      dataname <- dataset$dataname
      if (is.null(dataname) || !nzchar(dataname)) {
        dataset_ids <- unique(dataset$sampInfo$dataset)
        dataset_ids <- dataset_ids[!is.na(dataset_ids)]
        dataname <- if (length(dataset_ids)) dataset_ids[[1]] else "validation"
      }
      setNames(list(dataset), dataname)
    }
  )
)

c(
  targets_list,
  COMMON_DESURV_TARGETS
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
  
  
