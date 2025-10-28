create_param_grid_coldstarts = function(){
  # note alpha is not included here because of the warm start

  
  param_tbl <- tidyr::expand_grid(k       = K_VALS, 
                                  lambda  = LAMBDA_VALS,
                                  eta     = ETA_VALS,
                                  lambdaW = LAMBDAW_VALS,
                                  lambdaH = LAMBDAH_VALS,
                                  alpha = ALPHA_VALS) #%>%
    # mutate(id = sprintf("%s_m=%s_ng=%d_mi=%d_t=%.0e_imi=%d_i=%d_k=%d_l=%.0e_e=%.0e_lW=%.0e_lH=%.0e", 
    #                     train_prefix, method_trans_train, ngene, maxit, 
    #                     tol, imaxit, init, k, lambda, eta, lambdaW, lambdaH))
  
  # note by splitting over named list, targets should only run new combos
  # param_list <- split(param_tbl, param_tbl$id)  # each element is a 1-row tibble
  
  return(param_tbl)
}
