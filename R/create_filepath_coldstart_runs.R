# Stable path per (combo, seed [, knobs])
create_filepath_coldstart_runs <- function(params) {
  
  train_prefix       = TRAIN_PREFIX
  method_trans_train = METHOD_TRANS_TRAIN
  ngene              = NGENE
  maxit              = MAXIT
  tol                = TOL
  imaxit             = IMAXIT
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  alpha              = params$alpha
  seed               = params$seed
  method_select      = params$method_select
  
  folder = sprintf("ng%d_imax%d_tol%.0e_max%d",ngene,imaxit,tol,maxit)
  
  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_a%.4f_i%d",
                     k, lambda, eta, lambdaW, lambdaH, alpha, seed)
  
  dir=file.path("results", paste0("PKG_VERSION=",PKG_VERSION), paste0("GIT_BRANCH=",GIT_BRANCH),train_prefix, method_trans_train, folder, 
                "model_runs", method_select)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}
