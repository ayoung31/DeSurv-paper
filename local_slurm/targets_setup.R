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

default_controller = crew_controller_sequential()

# Local multicore controller for low memory tasks
low_mem_controller = crew_controller_local(
  name = "low_mem",
  workers = 4,
  seconds_idle = 120
)

# Local multicore controller for CV tasks (main BO computation)
# Limited workers to avoid memory exhaustion (BO tasks are memory-intensive)
cv_comp_controller = crew_controller_local(
  name = "cv",
  workers = 2,  # Run 2 BO tasks concurrently (each is memory intensive)
  seconds_idle = 300  # Longer idle time for long-running tasks
)

# Local multicore controller for full model runs
full_run_controller = crew_controller_local(
  name = "full",
  workers = 4,  # Multiple seed fits can run in parallel
  seconds_idle = 300
)

# Local multicore controller for medium memory tasks
med_mem_controller = crew_controller_local(
  name = "med_mem",
  workers = 4,
  seconds_idle = 120
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
