# Stable path per (combo, seed [, knobs])
create_filepath_init_alpha0 <- function(param_grid) {
  
  p=param_grid

  train_prefix = TRAIN_PREFIX
  method_trans_train = METHOD_TRANS_TRAIN
  ngene = NGENE
  maxit = MAXIT
  tol = TOL
  imaxit = IMAXIT
  k       = p$k
  lambda  = p$lambda
  eta     = p$eta
  lambdaW = p$lambdaW
  lambdaH = p$lambdaH

  folder = sprintf("ng%d_imax%d_tol%.0e_max%d",ngene,imaxit,tol,maxit)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e",
                     k, lambda, eta, lambdaW, lambdaH)
  
  dir=file.path("results", paste0("PKG_VERSION=",PKG_VERSION), paste0("GIT_BRANCH=",GIT_BRANCH),train_prefix, method_trans_train, folder, "alpha0_inits")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}
