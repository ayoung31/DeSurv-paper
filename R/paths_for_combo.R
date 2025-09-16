paths_for_combo <- function(version, method_norm, ngene, imaxit, tol, maxit,
                            k, lambda, eta, lambdaW, lambdaH) {
  mid = digest::digest(list(ng=ngene,imax=imaxit,tol=tol,max=maxit), algo = "md5")
  cid <- digest::digest(list(k = k, lambda = lambda, eta = eta, 
                             lambdaW = lambdaW, lambdaH = lambdaH), 
                        algo = "md5")
  root <- file.path("results", version, method_norm, mid, cid)
  list(
    root = root,
    # single-fit artifacts (kept for continuity)
    alpha0_inits = file.path(root, "alpha0", "inits_bundle.rds"),
    best_init_meta = file.path(root, "alpha0", "best_init_meta.rds"),
    warm_bundle = file.path(root, "warm", "warmstart_bundle.rds"),
    metrics = file.path(root, "metrics", "metrics.rds")
  )
}