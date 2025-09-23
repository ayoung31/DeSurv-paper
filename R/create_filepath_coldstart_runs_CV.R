# Stable path per (combo, seed [, knobs])
create_filepath_coldstart_runs_CV <- function(params) {

  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  fold               = params$fold

  folder = sprintf("ng%d_tol%.0e_max%d",NGENE,TOL,MAXIT)

  filename = sprintf("k%d_l%.0e_e%.0e_w%.0e_h%.0e_f%d_of_%d",
                     k, lambda, eta, lambdaW, lambdaH, fold, NFOLD)

  dir=file.path("results", paste0("PKG_VERSION=",PKG_VERSION,"_GIT_BRANCH=",GIT_BRANCH), TRAIN_PREFIX, METHOD_TRANS_TRAIN, folder, 
                "model_runs_cv_cold")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(filename, ".rds"))
}