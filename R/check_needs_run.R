check_needs_run = function(param_grid,ngene,alpha,ninit) {
  path_fits = create_filepath_warmstart_runs_CV(
    params   = param_grid,
    ngene    = ngene,
    tol      = TOL, 
    maxit    = MAXIT,
    nfold    = NFOLD,
    pkg_version      = PKG_VERSION, 
    git_branch       = GIT_BRANCH, 
    train_prefix     = TRAIN_PREFIX,
    method_trans_train = METHOD_TRANS_TRAIN
  )
  
  needs_run = TRUE
  needs_run = !is_valid_fit(path=path_fits,alpha=alpha,ninit=ninit)
  
  needs_run
}
