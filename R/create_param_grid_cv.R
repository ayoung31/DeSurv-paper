create_param_grid_cv = function(k,lambda,eta,lambdaW,lambdaH,nfold){
  # note alpha is not included here because of the warm start

  
  param_tbl <- tidyr::expand_grid(k       = k, 
                                  lambda  = lambda,
                                  eta     = eta,
                                  lambdaW = lambdaW,
                                  lambdaH = lambdaH,
                                  fold    = 1:nfold) #%>%
    # mutate(id = sprintf("%s_m=%s_ng=%d_mi=%d_t=%.0e_imi=%d_i=%d_k=%d_l=%.0e_e=%.0e_lW=%.0e_lH=%.0e", 
    #                     train_prefix, method_trans_train, ngene, maxit, 
    #                     tol, imaxit, init, k, lambda, eta, lambdaW, lambdaH))
  
  # note by splitting over named list, targets should only run new combos
  # param_list <- split(param_tbl, param_tbl$id)  # each element is a 1-row tibble
  
  return(param_tbl)
}
