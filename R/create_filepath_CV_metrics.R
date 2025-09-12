# Stable path per (combo, seed [, knobs])
create_filepath_CV_metrics <- function(params,train_prefix,method_trans_train, method_select) {

  train_prefix       = train_prefix
  method_trans_train = method_trans_train
  ngene              = params$ngene
  maxit              = params$maxit
  tol                = params$tol
  imaxit             = params$imaxit
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  method_select      = method_select
  fold               = params$fold

  folder = sprintf("ng%d_imax%d_tol%.0e_max%d",ngene,imaxit,tol,maxit)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_f%d",
                     k, lambda, eta, lambdaW, lambdaH, fold)

  dir=file.path("results", train_prefix, method_trans_train, folder, 
                "cv_metrics", method_select)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}