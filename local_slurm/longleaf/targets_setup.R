# targets_setup.R - UNC Longleaf HPC mode
# Matched to student's master branch for exact reproduction,
# plus nmf controller for targets_common_pipeline.R compatibility.
#
# Controller assignments:
#   cv       - BO targets: mclapply with ncores_grid=50
#   full     - Full-model seed fits: ninit_full=100 CPUs, 32 GB
#   med_mem  - Medium memory targets
#   low_mem  - Light tasks
#   nmf      - NMF k-selection (kept for targets_common_pipeline.R compatibility)
#   default  - Sequential fallback
#
# Resource budget (Longleaf):
#   cv:      202 workers × 50 CPUs (BO with parallel CV)
#   full:    202 workers × 100 CPUs × 32 GB (seed fits)
#   med_mem: 120 workers × 1 CPU × 8 GB
#   low_mem: 202 workers × 1 CPU × 1 GB
#   nmf:      10 workers × 10 CPUs (NMF k-selection)

library(targets)
library(tarchetypes)
LOCAL_RENDER <- identical(Sys.getenv("DESURV_LOCAL_RENDER"), "1")
library(crew)
if (!LOCAL_RENDER) {
  library(crew.cluster)
}
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION <- "HEAD"

GIT_BRANCH <- tryCatch(
  gert::git_branch(),
  error = function(e) {
    branch <- tryCatch(
      trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE, ignore.stderr = TRUE)),
      error = function(e2) "unknown"
    )
    if (length(branch) == 0 || !nzchar(branch)) branch <- "unknown"
    branch
  }
)

DEFAULT_NINIT <- 50
DEFAULT_NINIT_FULL <- 100

# ------ Slurm controllers (Longleaf HPC) ------
# Longleaf uses module system for R; script_lines ensures each worker loads R.
SLURM_SCRIPT_LINES <- c("module load r/4.4.0")

default_controller <- crew_controller_sequential()

if (LOCAL_RENDER) {
  active_controller <- default_controller
} else {
  low_mem_controller <- crew_controller_slurm(
    name = "low_mem",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_per_cpu = 1,
      time_minutes = 120,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  cv_comp_controller <- crew_controller_slurm(
    name = "cv",
    workers = 202,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      memory_gigabytes_per_cpu = 2,
      cpus_per_task = DEFAULT_NINIT,
      time_minutes = 720,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  full_run_controller <- crew_controller_slurm(
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
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  med_mem_controller <- crew_controller_slurm(
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
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  # NMF controller: not in student's master, but needed by targets_common_pipeline.R
  # which assigns fit_std to controller "nmf". Kept for compatibility (Option A).
  nmf_controller <- crew_controller_slurm(
    name = "nmf",
    workers = 10,
    seconds_idle = 120,
    seconds_interval = 0.25,
    options_cluster = crew_options_slurm(
      cpus_per_task = 10,
      memory_gigabytes_per_cpu = 1,
      time_minutes = 60,
      log_error = "logs/crew_log_%A.err",
      log_output = "logs/crew_log_%A.out",
      script_lines = SLURM_SCRIPT_LINES
    )
  )

  active_controller <- crew_controller_group(
    default_controller,
    low_mem_controller,
    cv_comp_controller,
    full_run_controller,
    med_mem_controller,
    nmf_controller
  )
}

# ---- Global options ----
TARGET_PACKAGES <- c(
  "clusterProfiler","org.Hs.eg.db","DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr",
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2"
)

tar_option_set(
  packages = TARGET_PACKAGES,
  format = "rds",
  controller = active_controller,
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# Load subtype data if available
CMB_SUBTYPES_PATH <- "data/derv/cmbSubtypes_formatted.RData"
if (file.exists(CMB_SUBTYPES_PATH)) {
  load(CMB_SUBTYPES_PATH)
}
