collect_W_H = function(best_params,data_folds,ngene,tol,maxit,nfold,pkg_version,
                       git_branch,train_prefix,method_trans_train,ninit){
  Ws=list()
  Hs=list()
  k=1
  for(i in 1:nfold){
    best_params$fold=i
    path_fits = create_filepath_warmstart_runs_CV(params=best_params,
                                                  ngene=ngene, tol=TOL, 
                                                  maxit=MAXIT, nfold=NFOLD,
                                                  pkg_version=PKG_VERSION, 
                                                  git_branch=GIT_BRANCH, 
                                                  train_prefix=TRAIN_PREFIX,
                                                  method_trans_train=METHOD_TRANS_TRAIN)
    bundle=readRDS(path_fits)
    for(j in 1:ninit){
      fit=bundle[[j]]$fits[[as.character(best_params$alpha)]]
      Ws[[k]] = fit$W
      Hs[[k]] = t(fit$H)
      rownames(Hs[[k]]) = colnames(data_folds$data_train[[i]]$ex)
      k=k+1
    }
  }
  return(list(Ws=Ws,Hs=Hs))
}

