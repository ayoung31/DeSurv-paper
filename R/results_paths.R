results_root_dir <- function(
    ngene,
    tol,
    maxit,
    pkg_version,
    git_branch,
    train_prefix,
    method_trans_train
) {
  folder <- sprintf("ng%d_tol%.0e_max%d", ngene, tol, maxit)
  root <- file.path(
    "results",
    sprintf("PKG_VERSION=%s_GIT_BRANCH=%s", pkg_version, git_branch),
    train_prefix,
    method_trans_train,
    folder
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}
