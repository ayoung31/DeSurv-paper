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

# ---- Global options ----
tar_option_set(
  packages = c("coxNMF","tidyverse","survival","cvwrapr","rmarkdown","dplyr"),
  format = "rds",
  controller = crew_controller_group(default_controller, model_runs_controller,easy_comp_controller),
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# ---- Reusable constants ----

VAL_DATASETS       = c("CPTAC","Dijk","Linehan","Moffitt_GEO_array",
                       "PACA_AU_array","PACA_AU_seq","Puleo_array")
METHOD_SELECT_INIT = "surv" ### method for selecting best initialization surv=PL, nmf=recon
ALPHA_VALS         = seq(0, .95, by = .05)
N_SHARDS = 200

# ---- Training parameters ----
TRAIN_DATASETS     = c("TCGA_PAAD")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
METHOD_TRANS_TRAIN = "quant"
NGENE              = 1000
IMAXIT             = 500
TOL                = 1e-6
MAXIT              = 6000
NINIT              = 30
K_VALS             = 2:12#2:16   #= c(2,3,4,5)
LAMBDA_VALS        = .1#10^seq(-3,3)#10^seq(-4,4)
ETA_VALS           = .1#c(.1,.5,.9)#seq(.1,.9,by=.1)
LAMBDAW_VALS       = 100#10^seq(-3,3)#10^seq(-4,4)
LAMBDAH_VALS       = 1e-4#10^seq(-3,3)#10^seq(-4,4)
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
  
  tar_target(param_grid,
             create_param_grid(K_VALS = K_VALS,
                               LAMBDA_VALS = LAMBDA_VALS,
                               ETA_VALS = ETA_VALS,
                               LAMBDAW_VALS = LAMBDAW_VALS,
                               LAMBDAH_VALS = LAMBDAH_VALS)
             ),

  tar_target(grid_sharded, assign_shards(param_grid, N_SHARDS)),
  tar_target(shard_ids, sort(unique(grid_sharded$shard))),
  
  # run and save initializations for each parameter combo
  tar_target(
    metrics_files,
    {
      print("running inits...")
      
      g <- dplyr::filter(grid_sharded, shard == shard_ids)
      out <- vector("list", nrow(g))
      for (i in seq_len(nrow(g))) {
        p <- as.list(g[i, c("k","lambda","eta","lambdaW","lambdaH")])
        paths <- paths_for_combo(VERSION, NGENE, p$k, p$lambda, p$eta, p$lambdaW, p$lambdaH)
        out[[i]] <- process_combo(X = data_filtered$ex,
                                  y = data_filtered$sampInfo$time,
                                  delta = data_filtered$sampInfo$event,
                                  p = p, paths = paths)
      }
      
      path = create_filepath_init_alpha0(param_grid=param_grid)
      
      init_alpha0(
        X = data_filtered$ex,
        y = data_filtered$sampInfo$time,
        delta = data_filtered$sampInfo$event,
        param_grid = param_grid,
        path = path,
        NINIT = NINIT
      )
    },
    pattern = map(shard_ids),
    iteration = "list",
    format    = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "inits")
    )
  ),
  #
  # read in all initializations
  tar_target(
    alpha0_scores,
    {
      df = dplyr::bind_rows(lapply(alpha0_score_file, function(path) {
        readRDS(path)
      }))
      df
    }
  ),
  #
  tarchetypes::tar_group_by(
    alpha0_groups,                 # = grouped target name
    alpha0_scores,                 # = your tibble from reading *.rds
    groups = c(k, lambda, eta, lambdaW, lambdaH)
  ),
  
  tar_target(
    best_init_per_param_combo,
    {
      df = readRDS(alpha0_inits)
      select_best_init(df = df, method_select = METHOD_SELECT_INIT)
    },
    pattern   = map(alpha0_inits),
    iteration = "list"
  )


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
  


  #

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
