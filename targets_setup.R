# _targets.R
library(targets)
library(tarchetypes)  
library(crew)
library(crew.cluster)
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION        = "HEAD"#utils::packageDescription("DeSurv", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()

DEFAULT_NINIT <- if (exists("DEFAULT_NINIT", inherits = TRUE)) DEFAULT_NINIT else 50
DEFAULT_NINIT_FULL <- if (exists("DEFAULT_NINIT_FULL", inherits = TRUE)) DEFAULT_NINIT_FULL else 100

# ------ Slurm controllers ------
default_controller = crew_controller_sequential()

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
    cpus_per_task = DEFAULT_NINIT,
    time_minutes = 1440,
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
    cpus_per_task = DEFAULT_NINIT_FULL,
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

# ---- Global options ----
TARGET_PACKAGES = c(
  "clusterProfiler","org.Hs.eg.db","DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr","digest",
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2","caret"
)

tar_option_set(
  packages = TARGET_PACKAGES,
  format = "rds",
  controller = active_controller,
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)
load("data/derv/cmbSubtypes_formatted.RData")
