library(targets)
PKG_VERSION        = utils::packageDescription("DeSurv", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()
tar_config_set(store=paste0("store_PKG_VERSION=",PKG_VERSION,"_GIT_BRANCH=",GIT_BRANCH,"_bladder"))

tar_make(data_val,script="_targets_bladder.R")
