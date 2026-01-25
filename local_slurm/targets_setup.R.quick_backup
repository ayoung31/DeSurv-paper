# _targets.R - Local Slurm version for desktop (20 CPUs, 31GB RAM)
# This file is configured for local desktop Slurm clusters without module system
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION <- "HEAD" # utils::packageDescription("DeSurv", fields = "RemoteRef")

# Get git branch with fallback if gert is not installed
GIT_BRANCH <- tryCatch(

  gert::git_branch(),
  error = function(e) {
    # Fallback: try git command directly, or use "unknown"
    branch <- tryCatch(
      trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE, ignore.stderr = TRUE)),
      error = function(e2) "unknown"
    )
    if (length(branch) == 0 || !nzchar(branch)) branch <- "unknown"
    branch
  }
)

# Local desktop: reduced from 50/100 to 19 to fit within 20 CPU limit
DEFAULT_NINIT <- if (exists("DEFAULT_NINIT", inherits = TRUE)) DEFAULT_NINIT else 19
DEFAULT_NINIT_FULL <- if (exists("DEFAULT_NINIT_FULL", inherits = TRUE)) DEFAULT_NINIT_FULL else 19

# ------ Slurm controllers (local desktop version) ------
default_controller = crew_controller_sequential()

low_mem_controller = crew_controller_slurm(
  name = "low_mem",
  workers = 20,

  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 120,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out"
  )
)

cv_comp_controller = crew_controller_slurm(
  name = "cv",
  workers = 20,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    cpus_per_task = DEFAULT_NINIT,
    time_minutes = 1440,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out"
  )
)

full_run_controller = crew_controller_slurm(
  name = "full",
  workers = 20,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 16,
    cpus_per_task = DEFAULT_NINIT_FULL,
    time_minutes = 600,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out"
  )
)

med_mem_controller = crew_controller_slurm(
  name = "med_mem",
  workers = 20,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 8,
    cpus_per_task = 1,
    time_minutes = 200,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out"
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

# Load subtype data if available (optional - only needed for certain downstream analyses)
.cmb_subtypes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (file.exists(.cmb_subtypes_path)) {
  load(.cmb_subtypes_path)
} else {
  message("Note: ", .cmb_subtypes_path, " not found. Some analyses may be unavailable.")
}
rm(.cmb_subtypes_path)
