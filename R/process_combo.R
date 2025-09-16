process_combo = function(X,y,delta,p,paths){
  
  init_metrics <- init_alpha0(X, y, delta, p, paths$alpha0_inits, NINIT, IMAXIT, MAXIT, TOL)

  best   <- select_best_init(init_metrics,METHOD_TRANS_TRAIN)
  
  if (!file.exists(paths$warm_bundle)) {
    best <- readRDS(paths$best_init_meta)
    wb   <- run_warmstarts_over_alpha(X, y, delta, p, ALPHA_VEC, best$best_seed, IMAXIT, MAXIT, TOL)
    write_atomic_rds(wb, paths$warm_bundle)
  }
  if (!file.exists(paths$metrics)) {
    wb <- readRDS(paths$warm_bundle)
    mt <- compute_metrics(wb, X, y, delta, p)
    write_atomic_rds(mt, paths$metrics)
  }
  paths$metrics
}