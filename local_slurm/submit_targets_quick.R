# submit_targets_quick.R - Run pipeline with easy/quick configs for local testing
# This script temporarily uses the local_slurm configs, then restores originals

library(targets)

# Store paths for backup/restore
main_configs <- "targets_configs.R"
main_setup <- "targets_setup.R"
local_configs <- "local_slurm/targets_configs.R"
local_setup <- "local_slurm/targets_setup.R"

# Backup original files if they exist
backup_configs <- NULL
backup_setup <- NULL

if (file.exists(main_configs)) {

  backup_configs <- readLines(main_configs)
}
if (file.exists(main_setup)) {

  backup_setup <- readLines(main_setup)
}

# Copy local configs to main directory
message("Setting up quick mode configs...")
file.copy(local_configs, main_configs, overwrite = TRUE)
file.copy(local_setup, main_setup, overwrite = TRUE)

# Set up store path
PKG_VERSION <- tryCatch(

  utils::packageDescription("DeSurv", fields = "RemoteRef"),
  error = function(e) "HEAD"
)
if (is.null(PKG_VERSION) || is.na(PKG_VERSION)) PKG_VERSION <- "HEAD"

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

store_path <- paste0("store_PKG_VERSION=", PKG_VERSION, "_GIT_BRANCH=", GIT_BRANCH, "_quick")
tar_config_set(store = store_path)
message("Using store: ", store_path)

# Run the pipeline
result <- tryCatch({
  tar_make()
  "success"
}, error = function(e) {
  message("Pipeline error: ", e$message)
  "error"
})

# Restore original files
message("Restoring original configs...")
if (!is.null(backup_configs)) {
  writeLines(backup_configs, main_configs)
} else if (file.exists(main_configs)) {
  file.remove(main_configs)
}

if (!is.null(backup_setup)) {
  writeLines(backup_setup, main_setup)
} else if (file.exists(main_setup)) {
  file.remove(main_setup)
}

message("Quick mode run complete. Result: ", result)
