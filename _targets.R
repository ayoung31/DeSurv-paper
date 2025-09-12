# _targets.R
library(targets)
library(tarchetypes)  
library(crew.cluster)
library(crew)
suppressWarnings(suppressMessages(library(dplyr)))

# ------ Slurm controllers ------
default_controller = crew_controller_sequential()

model_runs_controller = crew_controller_slurm(
  name = "model_runs",
  workers = 200,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 300,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

easy_comp_controller = crew_controller_slurm(
  name = "inits",
  workers = 200,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 300,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

val_comp_controller = crew_controller_slurm(
  name = "cv_validation",
  workers = 300,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 300,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

# ---- Global options ----
tar_option_set(
  packages = c("coxNMF","tidyverse","survival","cvwrapr","rmarkdown","dplyr"),
  format = "rds",
  controller = crew_controller_group(default_controller, 
                                     model_runs_controller,
                                     easy_comp_controller,
                                     val_comp_controller),
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# ---- Reusable constants ----

VAL_DATASETS       = c("CPTAC","Dijk","Linehan","Moffitt_GEO_array",
                       "PACA_AU_array","PACA_AU_seq","Puleo_array")
METHOD_SELECT_INIT = "surv" ### method for selecting best initialization surv=PL, nmf=recon
ALPHA_VALS         = seq(0, .95, by = .05)

# ---- Training parameters ----
TRAIN_DATASETS     = c("TCGA_PAAD")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
METHOD_TRANS_TRAIN = "quant"
NGENE              = 1000
IMAXIT             = 6000
TOL                = 1e-6
MAXIT              = 6000
NINIT              = 10
K_VALS             = 2:12#2:16   #= c(2,3,4,5)
LAMBDA_VALS        = 10^seq(-3,3)#10^seq(-4,4)
ETA_VALS           = c(.1,.5,.9)#seq(.1,.9,by=.1)
LAMBDAW_VALS       = 10^seq(-3,3)#10^seq(-4,4)
LAMBDAH_VALS       = 10^seq(-3,3)#10^seq(-4,4)
NFOLD              = 5
NTOP               = 25



# ---- Targets ----
list(
  ##### Training #####
  # Read data
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

  # Preprocess data
  tar_target(data_filtered, preprocess_data(data = data,
                                            ngene = NGENE,
                                            method_trans_train = METHOD_TRANS_TRAIN)),
  
  # parameters
  tar_target(alpha,
             {
               vals = ALPHA_VALS
               if(!(0 %in% vals)){
                 vals = c(0,vals)
               }
               o = order(vals)
               vals[o]
             }),
  
  # tar_target(param_grid,
  #            create_param_grid(TRAIN_PREFIX = TRAIN_PREFIX, 
  #                              METHOD_TRANS_TRAIN = METHOD_TRANS_TRAIN, 
  #                              NGENE = NGENE, 
  #                              MAXIT = MAXIT, 
  #                              TOL = TOL, 
  #                              IMAXIT = IMAXIT, 
  #                              K_VALS = K_VALS, 
  #                              LAMBDA_VALS = LAMBDA_VALS, 
  #                              ETA_VALS = ETA_VALS,
  #                              LAMBDAW_VALS = LAMBDAW_VALS, 
  #                              LAMBDAH_VALS = LAMBDAH_VALS)
  #            ),
  

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
  
  ### diagnostics for model runs ###
  
  ### training metrics plots ####
  # tar_target(
  #   
  # )
  
  # Summarise model metrics
  # tarchetypes::tar_render(
  #   report_metrics,
  #   path        = "reports/metrics.Rmd",   # create this file below
  #   output_file = "metrics.html",
  #   output_dir  = "reports/_site"
  # )
  
  ##### Cross Validation #####
  
  # split data into folds
  tar_target(data_folds,
             set_folds(data = data_filtered, nfold = NFOLD)
             ),
  
  # create param grid
  tar_target(param_grid_CV,
             create_param_grid_CV(TRAIN_PREFIX = TRAIN_PREFIX, 
                               METHOD_TRANS_TRAIN = METHOD_TRANS_TRAIN, 
                               NGENE = NGENE, 
                               MAXIT = MAXIT, 
                               TOL = TOL, 
                               IMAXIT = IMAXIT, 
                               K_VALS = K_VALS, 
                               LAMBDA_VALS = LAMBDA_VALS, 
                               ETA_VALS = ETA_VALS,
                               LAMBDAW_VALS = LAMBDAW_VALS, 
                               LAMBDAH_VALS = LAMBDAH_VALS, 
                               NFOLD = NFOLD)
  ), 
  
  # initializations
  tar_target(
    alpha0_inits_CV,
    {
      print("running inits...")
      
      path = create_filepath_init_alpha0_CV(param_grid=param_grid_CV)
      
      data_train = data_folds$data_train[[param_grid_CV$fold]]
      
      init_alpha0_CV(
        X = data_train$ex, 
        y = data_train$sampInfo$time, 
        delta = data_train$sampInfo$event,
        param_grid = param_grid_CV,
        path = path,
        NINIT = NINIT
      )
    },
    pattern   = map(param_grid_CV),
    format    = "file",
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "inits")
    ),
    cue = tar_cue(mode = "never")
  ),
  # 
  #select best initializations
  tar_target(
    best_init_per_param_combo_CV,
    {
      df = readRDS(alpha0_inits_CV)
      select_best_init(df = df, method_select = METHOD_SELECT_INIT)
    },
    pattern   = map(alpha0_inits_CV),
    iteration = "list",
    cue = tar_cue(mode = "never")
  ),
  
  tar_target(
    best_inits_cv_feasible_params,
    {
      inits <- dplyr::bind_rows(best_init_per_param_combo_CV, .id = "src")
      
      # 2) Validate names early (paranoid but useful)
      inits <- tibble::as_tibble(inits, .name_repair = "check_unique")
      
      # 3) Filter feasible and add row ids
      feasible <- dplyr::filter(inits, !flag_nan)
      feasible <- dplyr::mutate(feasible, id = dplyr::row_number())
    }
  ),
  
  # 
  #run warm starts
  tar_target(
    warmstarts_files_CV,
    {
      print("running warmstarts...")

      fold = best_inits_cv_feasible_params$fold

      path = create_filepath_warmstart_runs_CV(params = best_inits_cv_feasible_params)

      data_train = data_folds$data_train[[fold]]

      run_warmstarts_cv(
        X = data_train$ex, y = data_train$sampInfo$time, delta = data_train$sampInfo$event,
        params = best_inits_cv_feasible_params,
        alpha_vec = alpha,
        verbose = FALSE,
        path = path
      )

    },
    pattern   = map(best_inits_cv_feasible_params),
    iteration = "list",
    format = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "model_runs")
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    CV_metrics_full_file,
    {
      compute_metrics_CV(path = warmstarts_files_CV,
                         data_folds = data_folds,
                         ntop = NTOP)
    },
    pattern = map(warmstarts_files_CV),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv_validation")
    )
  ),


  tar_target(
    CV_metrics,
    {
      mets = dplyr::bind_rows(CV_metrics_full_file)
      mets %>%
        group_by(alpha,lambda,eta,lambdaW,lambdaH) %>%
        summarise(bic_mean = mean(bic,na.rm=TRUE),
                  bic_sd = sd(bic,na.rm=TRUE)) %>%
        ungroup()

    }

  ),

  # find the param combo with min BIC
  tar_target(
    selected_params,
    CV_metrics %>%
    slice_min(order_by = bic_mean, n = 1, with_ties = FALSE) %>%
      left_join(param_grid_CV) %>%
      filter(fold == 1)
  )
  #
  # # # initialize full model run for selected params
  # tar_target(
  #   alpha0_inits,
  #   {
  #     print("running inits...")
  # 
  #     path = create_filepath_init_alpha0(param_grid=selected_params)
  # 
  # 
  #     init_alpha0(
  #       X = data_filtered$ex,
  #       y = data_filtered$sampInfo$time,
  #       delta = data_filtered$sampInfo$event,
  #       param_grid = selected_params,
  #       path = path,
  #       NINIT = NINIT
  #     )
  #   },
  #   pattern = map(selected_params),
  #   iteration = "list",
  #   format    = "file",
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "inits")
  #   )
  # ),
  # # 
  # tar_target(
  #   best_init_per_param_combo,
  #   {
  #     df = readRDS(alpha0_inits)
  #     select_best_init(df = df, method_select = METHOD_SELECT_INIT)
  #   },
  #   pattern   = map(alpha0_inits),
  #   iteration = "list"
  # ),
  # # 
  # tar_target(
  #   warmstarts_files,
  #   {
  #     print("running warmstarts...")
  # 
  #     path = create_filepath_warmstart_runs(params = best_init_per_param_combo)
  # 
  #     run_warmstarts(
  #       X = data_filtered$ex, y = data_filtered$sampInfo$time, delta = data_filtered$sampInfo$event,
  #       params = best_init_per_param_combo,
  #       alpha_vec = alpha,
  #       verbose = FALSE,
  #       path = path
  #     )
  # 
  #   },
  #   pattern   = map(best_init_per_param_combo),
  #   iteration = "list",
  #   format = "file",
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "model_runs")
  #   )
  # ),
  # # 
  # # Compute model metrics
  # tar_target(
  #   training_metrics,
  #   compute_metrics(warmstarts_files,
  #                   data_filtered$ex,
  #                   data_filtered$sampInfo$time,
  #                   data_filtered$sampInfo$event),
  #   pattern   = map(warmstarts_files),
  #   iteration = "list"
  # ),
  # 
  # tar_target(
  #   training_metrics_table,
  #   dplyr::bind_rows(training_metrics)
  # )
  
  
  ##### Summaries #####
)
