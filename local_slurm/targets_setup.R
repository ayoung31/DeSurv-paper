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

# ------ Local multicore controllers (bypassing Slurm due to accounting issues) ------
# Using crew_controller_local instead of crew_controller_slurm
# This runs tasks locally with multiple processes
#
# Memory budget (20 CPUs, 30GB RAM):
#   OS/system:  ~3 GB reserved
#   cv:         2 workers × ~1.5 GB (parent + 5 mclapply forks via ncores_grid) = 3 GB, 12 CPUs
#   default:    4 workers × ~0.5 GB (figures, clustering, validation) = 2 GB, 4 CPUs
#   med_mem:    2 workers × ~1.0 GB (seed fits, full model runs) = 2 GB, 2 CPUs
#   full:       2 workers × ~0.5 GB = 1 GB, 2 CPUs
#   low_mem:    2 workers × ~0.3 GB = 0.6 GB
#   Peak estimate: ~12 GB (2× safety = 24 GB, under 30 GB limit)
#
# CPU constraint: each cv worker forks ncores_grid=5 sub-processes via
# parallel::mclapply inside DeSurv. 2 cv workers = 12 CPUs for BO alone.
# 3 cv workers would use 18/20 CPUs, starving the system.

default_controller = crew_controller_local(
  name = "default",
  workers = 2,  # Reduced from 4 to prevent OOM crashes (val_cindex_desurv_bladder)
  seconds_idle = Inf  # Never timeout: workers stay alive for entire pipeline run.
                      # Idle workers use no CPU and minimal RSS. tar_make() kills all on exit.
)

# Local multicore controller for low memory tasks
low_mem_controller = crew_controller_local(
  name = "low_mem",
  workers = 2,
  seconds_idle = Inf  # Never timeout (see default_controller comment)
)

# Local multicore controller for CV tasks (main BO computation)
# Each worker forks ncores_grid=5 sub-processes (parallel::mclapply),
# so 2 workers = 12 CPUs. Do not increase without reducing ncores_grid.
cv_comp_controller = crew_controller_local(
  name = "cv",
  workers = 1,  # Reduced from 2: each worker forks 5 mclapply children (6 procs total).
                # 2 workers = 12 procs, exhausts inotify max_user_instances (128) with system.
  seconds_idle = Inf  # Never timeout (see default_controller comment)
)

# Local multicore controller for full model runs
full_run_controller = crew_controller_local(
  name = "full",
  workers = 2,  # Seed fits are memory-heavy (~1 GB each)
  seconds_idle = Inf  # Never timeout (see default_controller comment)
)

# Local multicore controller for medium memory tasks (seed fits, consensus)
med_mem_controller = crew_controller_local(
  name = "med_mem",
  workers = 2,  # Each seed fit loads full expression matrix (~1 GB)
  seconds_idle = Inf  # Never timeout (see default_controller comment)
)

active_controller <- crew_controller_group(default_controller,
                                           low_mem_controller,
                                           cv_comp_controller,
                                           full_run_controller,
                                           med_mem_controller)

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
