# Stable path per (combo, seed [, knobs])
create_filepath_warmstart_runs_CV <- function(params) {

  train_prefix       = params$train_prefix
  method_trans_train = params$method_trans_train
  ngene              = params$ngene
  maxit              = params$maxit
  tol                = params$tol
  imaxit             = params$imaxit
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  seed               = params$seed
  method_select      = params$method_select
  fold               = params$fold

  folder = sprintf("ng%d_imax%d_tol%.0e_max%d",ngene,imaxit,tol,maxit)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_i%d_f%d",
                     k, lambda, eta, lambdaW, lambdaH, seed, fold)

  dir=file.path("results", train_prefix, method_trans_train, folder, 
                "model_runs_cv", method_select)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}