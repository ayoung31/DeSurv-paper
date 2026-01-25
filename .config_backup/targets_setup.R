# _targets.R
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
suppressWarnings(suppressMessages(library(dplyr)))

PKG_VERSION <- tryCatch(
  as.character(utils::packageVersion("DeSurv")),
  error = function(e) {
    "unknown"
  }
)

DESURV_REPO_PATH <- normalizePath("../DeSurv", winslash = "/", mustWork = FALSE)

DESURV_GIT_BRANCH <- tryCatch(
  {
    if (!dir.exists(DESURV_REPO_PATH)) stop("DeSurv repo not found")
    gert::git_branch(repo = DESURV_REPO_PATH)
  },
  error = function(e) {
    branch <- tryCatch(
      trimws(system(
        sprintf("git -C %s rev-parse --abbrev-ref HEAD", shQuote(DESURV_REPO_PATH)),
        intern = TRUE,
        ignore.stderr = TRUE
      )),
      error = function(e2) "unknown"
    )
    if (length(branch) == 0 || !nzchar(branch)) branch <- "unknown"
    branch
  }
)

DESURV_GIT_COMMIT <- tryCatch(
  {
    if (!dir.exists(DESURV_REPO_PATH)) stop("DeSurv repo not found")
    log <- gert::git_log(repo = DESURV_REPO_PATH, max = 1)
    if (nrow(log)) log$commit[[1]] else "unknown"
  },
  error = function(e) {
    commit <- tryCatch(
      trimws(system(
        sprintf("git -C %s rev-parse HEAD", shQuote(DESURV_REPO_PATH)),
        intern = TRUE,
        ignore.stderr = TRUE
      )),
      error = function(e2) "unknown"
    )
    if (length(commit) == 0 || !nzchar(commit)) commit <- "unknown"
    commit
  }
)

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

GIT_COMMIT <- tryCatch(
  {
    log <- gert::git_log(max = 1)
    if (nrow(log)) log$commit[[1]] else "unknown"
  },
  error = function(e) {
    commit <- tryCatch(
      trimws(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE)),
      error = function(e2) "unknown"
    )
    if (length(commit) == 0 || !nzchar(commit)) commit <- "unknown"
    commit
  }
)

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
  "parallel","foreach","doParallel","doMC","pec","glmnet","webshot2","caret","ggrepel","survminer","cowplot"
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
