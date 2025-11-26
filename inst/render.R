#' Render the manuscript via the targets pipeline.
library(targets)

tar_config_set(store = "store_PKG_VERSION=HEAD_GIT_BRANCH=new_init")
tar_make(names = "paper")
    
