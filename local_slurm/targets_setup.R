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

# Get DeSurv package git info (from installed package location)
DESURV_GIT_BRANCH <- tryCatch({
  pkg_path <- system.file(package = "DeSurv")
  if (nzchar(pkg_path) && dir.exists(file.path(pkg_path, ".git"))) {
    trimws(system(paste("git -C", shQuote(pkg_path), "rev-parse --abbrev-ref HEAD"),
                  intern = TRUE, ignore.stderr = TRUE))
  } else {
    # Package installed from git but no .git dir - use description
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

# Local desktop: reduced from 50/100 to 19 to fit within 20 CPU limit
DEFAULT_NINIT <- if (exists("DEFAULT_NINIT", inherits = TRUE)) DEFAULT_NINIT else 19
DEFAULT_NINIT_FULL <- if (exists("DEFAULT_NINIT_FULL", inherits = TRUE)) DEFAULT_NINIT_FULL else 19

# ------ Slurm controllers (resource-profiled for 20 CPU / 31 GB desktop) ------
# Slurm + cgroups enforces CPU/memory limits per worker.
# Workers queue until resources are available, preventing OOM crashes.
# seconds_idle = 300 for all: workers release Slurm resources after 5 min idle,
# freeing CPUs/RAM between pipeline phases (BO -> seed_fits -> downstream).
#
# Resource budget (from tar_meta() profiling of previous full run):
#   cv (BO):      6 CPUs × 6 GB = fits 3 concurrent (18 CPUs, 18 GB)
#   nmf:         10 CPUs × 10 GB = fits 2 concurrent (20 CPUs, 20 GB)
#   med_mem:      1 CPU  × 4 GB = fits 4 concurrent (4 CPUs, 16 GB)
#   default:      1 CPU  × 2 GB = fits 10 concurrent (10 CPUs, 20 GB)
# Slurm schedules dynamically — these are per-worker requests, not fixed reservations.

default_controller <- crew_controller_slurm(
  name = "default",
  workers = 10,
  seconds_idle = 300,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_per_cpu = 2,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# BO targets: desurv_cv_bayesopt uses mclapply(ncores_grid=5) internally.
# Each worker needs 6 CPUs (1 parent + 5 forks). Longest BO run: 4.2 hrs.
cv_comp_controller <- crew_controller_slurm(
  name = "cv",
  workers = 3,
  seconds_idle = 300,
  options_cluster = crew_options_slurm(
    cpus_per_task = 6,
    memory_gigabytes_per_cpu = 1,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# NMF k-selection: NMF::nmf(.options="p30") forks up to 30 processes.
# With 10 CPUs allocated, NMF parallelism is capped by cgroups.
# Only 2 of 6 fit_std targets are slow (20 min); rest are <6 sec.
nmf_controller <- crew_controller_slurm(
  name = "nmf",
  workers = 2,
  seconds_idle = 300,
  options_cluster = crew_options_slurm(
    cpus_per_task = 10,
    memory_gigabytes_per_cpu = 1,
    log_error = "logs/crew_%A.err",
    log_output = "logs/crew_%A.out"
  )
)

# Seed fits: sequential for-loop (ninit_full=100, parallel_init=FALSE).
# Single-threaded but memory-heavy (up to 350 MB output). Longest: 29 min.
med_mem_controller <- crew_controller_slurm(
  name = "med_mem",
  workers = 4,
  seconds_idle = 300,
  options_cluster = crew_options_slurm(
    cpus_per_task = 1,
    memory_gigabytes_required = 4,
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
TARGET_PACKAGES = c(
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

# Load subtype data if available (optional - only needed for certain downstream analyses)
.cmb_subtypes_path <- "data/derv/cmbSubtypes_formatted.RData"
if (file.exists(.cmb_subtypes_path)) {
  load(.cmb_subtypes_path)
} else {
  message("Note: ", .cmb_subtypes_path, " not found. Some analyses may be unavailable.")
}
rm(.cmb_subtypes_path)
