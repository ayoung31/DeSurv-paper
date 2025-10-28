# Stable path per (combo, seed [, knobs])
create_filepath_warmstart_runs_CV =
  function(params, ngene, tol, maxit, nfold,
           pkg_version, git_branch, train_prefix,
           method_trans_train) {

  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  fold               = params$fold

  folder = sprintf("ng%d_tol%.0e_max%d",ngene,tol,maxit)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_f%d_of_%d",
                     k, lambda, eta, lambdaW, lambdaH, fold, nfold)

  dir=file.path("results", paste0("PKG_VERSION=",pkg_version,"_GIT_BRANCH=",git_branch), 
                train_prefix, method_trans_train, folder, 
                "model_runs_cv")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}