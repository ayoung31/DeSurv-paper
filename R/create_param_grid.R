create_param_grid = function(K_VALS, LAMBDA_VALS, ETA_VALS, 
                             LAMBDAW_VALS, LAMBDAH_VALS){
  # note alpha is not included here because of the warm start

  
  param_tbl <- tidyr::expand_grid(k       = K_VALS, 
                                  lambda  = LAMBDA_VALS,
                                  eta     = ETA_VALS,
                                  lambdaW = LAMBDAW_VALS,
                                  lambdaH = LAMBDAH_VALS) %>%
    distinct()

  return(param_tbl)
}
