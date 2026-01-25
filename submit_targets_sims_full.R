# submit_targets_sims_full.R - Run simulation pipeline with full configs
# Store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full
library(targets)

STORE_NAME <- "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"
tar_config_set(store = STORE_NAME)
message("Using targets store: ", STORE_NAME)
message("Running simulation pipeline with FULL configs...")
message("This includes: bo, bo_alpha0, bo_tune_ntop, bo_tune_ntop_alpha0, fixed, fixed_alpha0")
tar_make(script = "_targets_sims.R")
