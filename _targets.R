# _targets.R
library(targets)
library(tarchetypes)  
library(crew.cluster)
library(crew)
suppressWarnings(suppressMessages(library(dplyr)))

NINIT=20

# ------ Slurm controllers ------
default_controller = crew_controller_sequential()

low_mem_controller = crew_controller_slurm(
  name = "low_mem",
  workers = 220,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 120,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

cv_comp_controller = crew_controller_slurm(
  name = "cv",
  workers = 80,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 12,
    cpus_per_task = NINIT,
    time_minutes = 600,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)


# ---- Global options ----
tar_option_set(
  packages = c("coxNMF","tidyverse","survival","cvwrapr","rmarkdown","dplyr",
               "parallel","foreach", "doParallel", "doMC"),
  format = "rds",
  controller = crew_controller_group(default_controller, 
                                     low_mem_controller,
                                     cv_comp_controller),
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# ---- Reusable constants ----

VAL_DATASETS       = c("CPTAC","Dijk","Linehan","Moffitt_GEO_array",
                       "PACA_AU_array","PACA_AU_seq","Puleo_array")
METHOD_SELECT_INIT = "surv" ### method for selecting best initialization surv=PL, nmf=recon
ALPHA         = seq(0, 1, by = .05)

# ---- Training parameters ----
TRAIN_DATASETS     = c("TCGA_PAAD","CPTAC")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
METHOD_TRANS_TRAIN = "rank"
NGENE              = 1000
TOL                = 1e-5
MAXIT              = 6000
K_VALS             = 2:12#2:16   #= c(2,3,4,5)
LAMBDA_VALS        = c(0,1e-3,1e-2,.1)#10^seq(-3,3)#10^seq(-4,4)
ETA_VALS           = c(0,.01)#c(.01,.1,.5,.9)#seq(.1,.9,by=.1)
LAMBDAW_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
LAMBDAH_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
NTOP               = 25
NFOLD              = 5
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
  
  tar_target(
    data_folds,
    {
      set_folds(data,NFOLD,seed=123)
    }
  ),

  # Preprocess data
  tar_target(data_filtered, preprocess_data(data = data,
                                            ngene = NGENE,
                                            method_trans_train = METHOD_TRANS_TRAIN)),


  tar_target(param_grid,
             create_param_grid_cv()
             ),
  
  tar_target(
    cv_runs,
    {
      print("running cv...")
      path_fits = create_filepath_warmstart_runs_CV(params=param_grid)
      
      f = param_grid$fold
      X = data_folds$data_train[[f]]$ex
      y = data_folds$data_train[[f]]$sampInfo$time
      d = data_folds$data_train[[f]]$sampInfo$event
      
      run_warmstarts_cv(X=X,y=y,delta=d,
                     params=param_grid,verbose=FALSE,
                     path_fits = path_fits)

    },
    pattern = map(param_grid),
    iteration = "list",
    format    = "file",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),

  tar_target(
    cv_metrics_list,
    {
      bundle = readRDS(cv_runs)
      path_mets = create_filepath_CV_metrics(params = param_grid)
      f = param_grid$fold
      k = param_grid$k
      lambda = param_grid$lambda
      mets_tr=list()
      mets_te=list()
      j=1
      for(i in 1:length(bundle)){
        b = bundle[[i]]
        for(a in names(b$fits)){
          alpha = as.numeric(a)
          fit = b$fits[[a]]

          X = data_folds$data_train[[f]]$ex
          y = data_folds$data_train[[f]]$sampInfo$time
          d = data_folds$data_train[[f]]$sampInfo$event

          mets_tr[[j]] = tryCatch({
            m = compute_metrics(fit,X,y,d,alpha,f,test=FALSE)
            m$k = k
            m$lambda = lambda
            m$seed = i
            m
          },error=function(e) NULL)

          Xtest = data_folds$data_test[[f]]$ex
          ytest = data_folds$data_test[[f]]$sampInfo$time
          dtest = data_folds$data_test[[f]]$sampInfo$event

          mets_te[[j]] = tryCatch({
            m = compute_metrics(fit,Xtest,ytest,dtest,alpha,f,test=TRUE)
            m$k = k
            m$lambda = lambda
            m$seed = i
            m
          },error=function(e) NULL)

          j=j+1
        }
      }
      mets_train = dplyr::bind_rows(mets_tr)
      mets_test = dplyr::bind_rows(mets_te)
      setNames(list(list(mets_train=mets_train, mets_test=mets_test)), 
               paste0(param_grid$k,param_grid$lambda,param_grid$eta,
                      param_grid$lambdaW,param_grid$lambdaH,param_grid$fold))
    },
    pattern=map(cv_runs,param_grid),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "low_mem")
    )
  ),

  tar_target(
    cv_metrics,
    {
      mets_train = dplyr::bind_rows(lapply(cv_metrics_list,`[[`,"mets_train"))
      mets_test = dplyr::bind_rows(lapply(cv_metrics_list,`[[`,"mets_test"))
      list(mets_train=mets_train,mets_test=mets_test)
    }

  )
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
)
