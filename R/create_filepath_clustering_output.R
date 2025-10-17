# Stable path per (combo, seed [, knobs])
create_filepath_clustering_output = 
  function(ngene, tol, maxit, pkg_version, 
           git_branch, train_prefix, method_trans_train, 
           ntop, val_dataset, method) {
  
  folder = sprintf("ng%d_tol%.0e_max%d",ngene,tol,maxit)
  
  dir=file.path("results", paste0("PKG_VERSION=",pkg_version,"_GIT_BRANCH=",git_branch), 
                train_prefix, method_trans_train, folder, 
                "clustering", val_dataset, paste0("ntop=",ntop), method)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}