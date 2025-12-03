# _targets.R
library(targets)
library(tarchetypes)  
LOCAL_RENDER <- identical(Sys.getenv("DESURV_LOCAL_RENDER"), "1")
library(crew)
if (!LOCAL_RENDER) {
  library(crew.cluster)
}
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION        = "HEAD"#utils::packageDescription("DeSurv", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()

pipeline_param <- function(name, default) {
  if (exists(name, inherits = TRUE)) {
    get(name, inherits = TRUE)
  } else {
    default
  }
}

NINIT <- pipeline_param("NINIT", 50)
NINIT_FULL <- pipeline_param("NINIT_FULL", 100)
BO_N_INIT <- pipeline_param("BO_N_INIT", 20)
BO_N_ITER <- pipeline_param("BO_N_ITER", 100)
BO_CANDIDATE_POOL <- pipeline_param("BO_CANDIDATE_POOL", 2000)
BO_MAX_REFINEMENTS <- pipeline_param("BO_MAX_REFINEMENTS", 2)
BO_TOL_GAIN <- pipeline_param("BO_TOL_GAIN", 0.002)
BO_PLATEAU <- pipeline_param("BO_PLATEAU", 1)
BO_TOP_K <- pipeline_param("BO_TOP_K", 10)
BO_SHRINK_BASE <- pipeline_param("BO_SHRINK_BASE", 0.5)
BO_IMPORTANCE_GAIN <- pipeline_param("BO_IMPORTANCE_GAIN", 0.3)
BO_COARSE_CONTROL <- pipeline_param(
  "BO_COARSE_CONTROL",
  list(
    n_init = BO_N_INIT,
    n_iter = BO_N_ITER,
    candidate_pool = BO_CANDIDATE_POOL,
    exploration_weight = 0.01,
    seed = 123,
    cv_verbose = FALSE
  )
)
BO_REFINE_CONTROL <- pipeline_param(
  "BO_REFINE_CONTROL",
  list(
    n_init = BO_N_INIT,
    n_iter = BO_N_ITER,
    candidate_pool = BO_CANDIDATE_POOL,
    exploration_weight = 0.01,
    seed = 456,
    cv_verbose = FALSE
  )
)

# ------ Slurm controllers ------
default_controller = crew_controller_sequential()

if (LOCAL_RENDER) {
  active_controller <- default_controller
} else {
  low_mem_controller = crew_controller_slurm(
    name = "low_mem",
    workers = 202,
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
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_per_cpu = 2,
      cpus_per_task = NINIT,
      time_minutes = 720,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = "module load r/4.4.0"
    )
  )
  
  full_run_controller = crew_controller_slurm(
    name = "full",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = 32,
      cpus_per_task = NINIT_FULL,
      time_minutes = 600,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = "module load r/4.4.0"
    )
  )
  
  med_mem_controller = crew_controller_slurm(
    name = "med_mem",
    workers = 120,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_required = 8,
      cpus_per_task = 1,
      time_minutes = 200,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = "module load r/4.4.0"
    )
  )
  
  active_controller <- crew_controller_group(default_controller, 
                                             low_mem_controller,
                                             cv_comp_controller,
                                             full_run_controller,
                                             med_mem_controller)
}

# ---- Global options ----
TARGET_PACKAGES = c(
  "clusterProfiler","org.Hs.eg.db","DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr",
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2"
)

tar_option_set(
  packages = TARGET_PACKAGES,
  format = "rds",
  controller = active_controller,
  error = "continue"
)


# ---- Training parameters ----
METHOD_TRANS_TRAIN = "rank"

NGENE_DEFAULT      = NGENE_CONFIG[[1]]
TUNE_NGENE         = length(unique(NGENE_CONFIG)) > 1

NTOP_DEFAULT       = NTOP_CONFIG[[1]]
TUNE_NTOP          = length(unique(NTOP_CONFIG)) > 1

LAMBDAW_DEFAULT    = LAMBDAW_CONFIG[[1]]
TUNE_LAMBDAW       = length(unique(LAMBDAW_CONFIG)) > 1

LAMBDAH_DEFAULT    = LAMBDAH_CONFIG[[1]]
TUNE_LAMBDAH       = length(unique(LAMBDAH_CONFIG)) > 1



# convergence
TOL                = 1e-5
MAXIT              = 4000

# standard NMF
STD_NMF_K_GRID     = 2:12#2:16   #= c(2,3,4,5)
COXNET_LAMBDA_GRID = c(1e-4,1e-3,1e-2,.1,1,10)#,seq(.2,.9,by=.1)   #10^seq(-3,3)#10^seq(-4,4)
COXNET_ALPHA_GRID  = seq(0,1,by=.1)#c(0,.01,.1,.5,.9)#c(0,.01)#c(.01,.1,.5,.9)#seq(.1,.9,by=.1)

# cross validation
NFOLD              = 5

DESURV_PARALLEL_GRID <- pipeline_param("DESURV_PARALLEL_GRID", TRUE)
DESURV_NCORES_GRID <- pipeline_param("DESURV_NCORES_GRID", NINIT)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)
load("data/derv/cmbSubtypes_formatted.RData")
