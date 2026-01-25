# submit_targets_sims.R - Local quick version
# Uses the quick simulation config with reduced settings for local desktop

library(targets)

# Use the same store as our main pipeline (naimedits0125 branch)
tar_config_set(store = "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125")

cat("Running simulation pipeline with quick settings...\n")
cat("Store:", tar_config_get("store"), "\n")
cat("Script: local_slurm/quick/_targets_sims.R\n\n")

tar_make(script = "local_slurm/quick/_targets_sims.R")
