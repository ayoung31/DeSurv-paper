create_param_grid_cv = function(k,lambda,eta,lambdaW,lambdaH,nfold){
  # note alpha is not included here because of the warm start

  
  param_tbl <- tidyr::expand_grid(k       = k, 
                                  lambda  = lambda,
                                  eta     = eta,
                                  lambdaW = lambdaW,
                                  lambdaH = lambdaH,
                                  fold    = 1:nfold) %>%
    mutate(id = paste0("k", k, "_l", lambda, "_e", eta, "_f", fold))
  
  # note by splitting over named list, targets should only run new combos
  param_list <- split(param_tbl, param_tbl$id)  # each element is a 1-row tibble
  
  return(param_list)
}
