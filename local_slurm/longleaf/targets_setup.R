# targets_setup.R - UNC Longleaf HPC mode
# Configured for Longleaf cluster with large worker pools and module system.
#
# Controller assignments (must match targets_common_pipeline.R):
#   cv       - BO targets (desurv_bo_results_*): mclapply with ncores_grid=50
#   nmf      - NMF k-selection (fit_std_*): NMF forks up to 30 procs
#   med_mem  - Seed fits (desurv_seed_fits_*, tar_fit_desurv_*): sequential, memory-heavy
#   default  - Everything else: figures, validation, ORA, clustering
#
# Resource budget (Longleaf: hundreds of cores, TBs of RAM):
#   cv:      202 workers × 50 CPUs = up to 10,100 CPUs (Slurm queues dynamically)
#   nmf:      10 workers × 10 CPUs = up to 100 CPUs
#   med_mem: 120 workers × 1 CPU × 8 GB = up to 960 GB
#   default:  50 workers × 1 CPU × 2 GB = up to 100 GB

library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
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

# Get git commit for paper repo
GIT_COMMIT <- tryCatch(
  gert::git_info()$commit,
  error = function(e) {
    commit <- tryCatch(
      trimws(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE)),
      error = function(e2) "unknown"
    )
    if (length(commit) == 0 || !nzchar(commit)) commit <- "unknown"
    commit
  }
)

# Get DeSurv package git info
DESURV_GIT_BRANCH <- tryCatch({
  pkg_path <- system.file(package = "DeSurv")
  if (nzchar(pkg_path) && dir.exists(file.path(pkg_path, ".git"))) {
    trimws(system(paste("git -C", shQuote(pkg_path), "rev-parse --abbrev-ref HEAD"),
                  intern = TRUE, ignore.stderr = TRUE))
  } else {
    desc <- utils::packageDescription("DeSurv")
    if (!is.null(desc$RemoteSha)) substr(desc$RemoteSha, 1, 7) else "installed"
  }
}, error = function(e) "unknown")

DESURV_GIT_COMMIT <- tryCatch({
  pkg_path <- system.file(package = "DeSurv")
  if (nzchar(pkg_path) && dir.exists(file.path(pkg_path, ".git"))) {
    trimws(system(paste("git -C", shQuote(pkg_path), "rev-parse HEAD"),
                  intern = TRUE, ignore.stderr = TRUE))
  } else {
    desc <- utils::packageDescription("DeSurv")
    if (!is.null(desc$RemoteSha)) desc$RemoteSha else "unknown"
  }
}, error = function(e) "unknown")

# Longleaf HPC: override NINIT values (uncapped by LOCAL_CPU_LIMIT since
# DESURV_CPU_LIMIT=200 is set in the sbatch script)
DEFAULT_NINIT <- 50
DEFAULT_NINIT_FULL <- 100

# ------ Slurm controllers (Longleaf HPC) ------
# Longleaf uses module system for R; script_lines ensures each worker loads R.
# seconds_idle = 120: release Slurm resources after 2 min idle.
SLURM_SCRIPT_LINES <- c("module load r/4.4.0")

default_controller <- crew_controller_slurm(
  name = "default",
  workers = 50,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_per_cpu = 2,
    time_minutes = 120,
    script_lines = SLURM_SCRIPT_LINES,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# BO targets: desurv_cv_bayesopt uses mclapply(ncores_grid=50) internally.
# Each worker needs 50 CPUs. Longest BO run: 4-12 hrs depending on dataset.
cv_comp_controller <- crew_controller_slurm(
  name = "cv",
  workers = 202,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = DEFAULT_NINIT,
    memory_gigabytes_per_cpu = 2,
    time_minutes = 720,
    script_lines = SLURM_SCRIPT_LINES,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# NMF k-selection: NMF::nmf(.options="p30") forks up to 30 processes.
# 10 CPUs allocated; NMF parallelism capped by cgroups.
nmf_controller <- crew_controller_slurm(
  name = "nmf",
  workers = 10,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = 10,
    memory_gigabytes_per_cpu = 1,
    time_minutes = 60,
    script_lines = SLURM_SCRIPT_LINES,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# Seed fits: sequential for-loop (ninit_full=100, parallel_init=FALSE).
# Single-threaded but memory-heavy (up to 350 MB output per target).
med_mem_controller <- crew_controller_slurm(
  name = "med_mem",
  workers = 120,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_required = 8,
    time_minutes = 240,
    script_lines = SLURM_SCRIPT_LINES,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

active_controller <- crew_controller_group(
  default_controller,
  cv_comp_controller,
  nmf_controller,
  med_mem_controller
)

# ---- Global options ----
TARGET_PACKAGES <- c(
  "clusterProfiler","org.Hs.eg.db","DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr","digest",
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2","caret","ggrepel","survminer"
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
.cmb_subtypes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (file.exists(.cmb_subtypes_path)) {
  load(.cmb_subtypes_path)
} else {
  message("Note: ", .cmb_subtypes_path, " not found. Some analyses may be unavailable.")
}
rm(.cmb_subtypes_path)
