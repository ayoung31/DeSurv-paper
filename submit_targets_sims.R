library(targets)
PKG_VERSION        = "HEAD" #utils::packageDescription("DeSurv", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()
tar_config_set(store=paste0("store_PKG_VERSION=",PKG_VERSION,"_GIT_BRANCH=",GIT_BRANCH,"_sims"))

tar_make(script="_targets_sims_full.R")

