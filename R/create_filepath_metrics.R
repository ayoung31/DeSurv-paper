# Stable path per (combo, seed [, knobs])
create_filepath_metrics = 
  function(params, ngene, tol, maxit, nfold, pkg_version, 
           git_branch, train_prefix, method_trans_train) {
  
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  
  folder = sprintf("ng%d_tol%.0e_max%d",ngene,tol,maxit)
  
  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e",
                     k, lambda, eta, lambdaW, lambdaH)
  
  dir=file.path("results", paste0("PKG_VERSION=",pkg_version,"_GIT_BRANCH=",git_branch), 
                train_prefix, method_trans_train, folder, 
                "metrics")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".RData"))
}