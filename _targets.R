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
NINIT              = 2
K_VALS             = 2:12#2:16   #= c(2,3,4,5)
LAMBDA_VALS        = c(0,.1,1)#10^seq(-3,3)#10^seq(-4,4)
ETA_VALS           = .01#c(.01,.1,.5,.9)#seq(.1,.9,by=.1)
LAMBDAW_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
LAMBDAH_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
NTOP               = 25
PKG_VERSION        = utils::packageDescription("coxNMF", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()



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


  tar_target(param_grid,
             create_param_grid_coldstarts()
             ),
  
  tar_target(
    inits,
    {
      print("running inits...")
      
      path = create_filepath_inits_coldstarts(param_grid=param_grid)
      
      
      init_coldstarts(
        X = data_filtered$ex,
        y = data_filtered$sampInfo$time,
        delta = data_filtered$sampInfo$event,
        param_grid = param_grid,
        path = path
      )
    },
    pattern = map(param_grid),
    iteration = "list",
    format    = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "inits")
    )
  ),

  tar_target(
    best_init_per_param_combo,
    {
      df = readRDS(inits)
      select_best_init(df = df, method_select = METHOD_SELECT_INIT)
    },
    pattern   = map(inits),
    iteration = "list"
  ),
  
  tar_target(
    model_runs,
    {
      print("running coldstarts...")

      path = create_filepath_coldstart_runs(params = best_init_per_param_combo)

      run_coldstarts(
        X = data_filtered$ex, y = data_filtered$sampInfo$time, delta = data_filtered$sampInfo$event,
        params = best_init_per_param_combo,
        verbose = FALSE,
        path = path
      )

    },
    pattern   = map(best_init_per_param_combo),
    iteration = "list",
    format = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "model_runs")
    )
  )
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
)
