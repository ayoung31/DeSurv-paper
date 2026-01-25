# submit_targets_full.R - Run main pipeline with full configs
# Store: store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full
library(targets)

STORE_NAME <- "store_PKG_VERSION=NA_GIT_BRANCH=naimedits0125_full"
tar_config_set(store = STORE_NAME)
message("Using targets store: ", STORE_NAME)
message("Running main pipeline with FULL configs...")
tar_make()
