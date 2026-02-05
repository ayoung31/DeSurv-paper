# Stable path per (combo, seed [, knobs])
create_filepath_clustering_output =
  function(ngene, tol, maxit, pkg_version,
           git_branch, train_prefix, method_trans_train,
           ntop, val_dataset, method) {
    
    stopifnot(
      length(ngene) == 1, !is.na(ngene),
      length(tol) == 1, !is.na(tol),
      length(maxit) == 1, !is.na(maxit),
      length(ntop) == 1, !is.na(ntop)
    )
    
    folder <- sanitize_path_component(
      sprintf("ng%d_tol%.0e_max%d", ngene, tol, maxit),
      "training folder"
    )
    pkg_branch <- sanitize_path_component(
      sprintf("PKG_VERSION=%s_GIT_BRANCH=%s", pkg_version, git_branch),
      "package/branch directory"
    )
    train_prefix <- sanitize_path_component(train_prefix, "training prefix")
    method_trans_train <- sanitize_path_component(
      method_trans_train,
      "training transform"
    )
    val_dataset <- sanitize_path_component(val_dataset, "validation dataset")
    ntop_dir <- sanitize_path_component(paste0("ntop=", ntop), "ntop directory")
    method <- sanitize_path_component(method, "clustering method")
    
    dir <- file.path(
      "results",
      pkg_branch,
      train_prefix,
      method_trans_train,
      folder,
      "clustering",
      val_dataset,
      ntop_dir,
      method
    )
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir
  }
