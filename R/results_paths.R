compact_path_parts <- function(parts) {
  keep <- vapply(
    parts,
    function(part) {
      if (is.null(part)) return(FALSE)
      if (!length(part)) return(FALSE)
      if (is.na(part)) return(FALSE)
      if (is.character(part)) return(nzchar(part))
      TRUE
    },
    logical(1)
  )
  parts[keep]
}

results_root_dir <- function(
    ngene,
    tol,
    maxit,
    pkg_version,
    git_branch,
    train_prefix,
    method_trans_train,
    config_tag = NULL,
    config_id = NULL
) {
  folder <- sprintf("ng%d_tol%.0e_max%d", ngene, tol, maxit)
  root_parts <- compact_path_parts(c(
    "results",
    sprintf("PKG_VERSION=%s_GIT_BRANCH=%s", pkg_version, git_branch),
    config_tag,
    config_id,
    train_prefix,
    method_trans_train,
    folder
  ))
  root <- do.call(file.path, as.list(root_parts))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}

get_bo_results_dir <- function(
    pkg_version,
    git_branch,
    train_prefix,
    method_trans_train,
    bo_config_tag = NULL,
    bo_config_id = NULL
) {
  root_parts <- compact_path_parts(c(
    "results",
    sprintf("PKG_VERSION=%s_GIT_BRANCH=%s", pkg_version, git_branch),
    "bo_runs",
    bo_config_tag,
    bo_config_id,
    train_prefix,
    method_trans_train
  ))
  root <- do.call(file.path, as.list(root_parts))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}
