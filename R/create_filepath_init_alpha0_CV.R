# Stable path per (combo, seed [, knobs])
create_filepath_init_alpha0_CV <- function(param_grid) {
  
  p=param_grid

  train_prefix = p$train_prefix
  method_trans_train = p$method_trans_train
  ngene = p$ngene
  maxit = p$maxit
  tol = p$tol
  imaxit = p$imaxit
  k       = p$k
  lambda  = p$lambda
  eta     = p$eta
  lambdaW = p$lambdaW
  lambdaH = p$lambdaH
  fold    = p$fold

  folder = sprintf("ng%d_imax%d_tol%.0e_max%d",ngene,imaxit,tol,maxit)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_f%d",
                     k, lambda, eta, lambdaW, lambdaH, fold)
  
  dir=file.path("results", train_prefix, method_trans_train, folder, "alpha0_inits_CV")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}
