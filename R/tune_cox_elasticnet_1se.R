tune_cox_elasticnet_1se <- function( 
  x, y_surv, 
  alpha_grid   = seq(0, 1, by = 0.2), 
  lambda_grid  = NULL,        # if NULL, coxnet picks its own decreasing sequence 
  nfolds       = 5, 
  type_measure = c("deviance", "C"), 
  parallel     = FALSE, 
  standardize  = TRUE, 
  seed         = 1, 
  foldid       = NULL         # optional pre-defined fold ids (1..nfolds) 
){ 
  
  type_measure <- match.arg(type_measure) 
  higher_is_better <- type_measure %in% c("C") 
  
  
  set.seed(seed) 
  
  n <- nrow(x) 
  if (is.null(foldid)) { 
    foldid <- sample(rep(1:nfolds, length.out = n)) 
  } 
  
  # store per-alpha cv objects and their 1SE choices 
  per_alpha_cv   <- list() 
  per_alpha_1se  <- list() 
  all_rows       <- list() 
  
  for (a in alpha_grid) { 
    
    cvfit <- cv.glmnet( 
      x, y_surv, 
      family       = "cox",
      alpha        = a, 
      lambda       = lambda_grid,   # can be NULL or a vector 
      nfolds       = nfolds, 
      foldid       = foldid,        # fixed folds across alphas
      type.measure = type_measure, 
      parallel     = parallel, 
      standardize  = standardize 
    ) 

    per_alpha_cv[[as.character(a)]] <- cvfit 
    
    df_alpha <- tibble( 
      alpha   = a, 
      lambda  = cvfit$lambda, 
      cv_mean = cvfit$cvm, 
      cv_se   = cvfit$cvsd 
    ) 
    
    all_rows[[length(all_rows) + 1]] <- df_alpha 
    
    
    # index at lambda.1se (match to grid safely) 
    idx_1se <- which.min(abs(df_alpha$lambda - cvfit$lambda.1se)) 
    row_1se <- df_alpha[idx_1se, ] 
    
    
    per_alpha_1se[[length(per_alpha_1se) + 1]] <- tibble( 
      alpha         = a, 
      lambda_1se    = df_alpha$lambda[idx_1se], 
      cv_mean_1se   = df_alpha$cv_mean[idx_1se], 
      cv_se_1se     = df_alpha$cv_se[idx_1se] 
    ) 
  } 
  
  
  cv_table_all  <- bind_rows(all_rows) 
  table_1se     <- bind_rows(per_alpha_1se) 
  
  
  # choose alpha by metric at its lambda.1se (use ties -> prefer larger lambda = more regularization) 
  
  if (higher_is_better) { 
    best_idx <- with(table_1se, order(-cv_mean_1se, -lambda_1se))[1] 
  } else { 
    best_idx <- with(table_1se, order(cv_mean_1se, -lambda_1se))[1] 
  } 
  
  
  
  best_alpha   <- table_1se$alpha[best_idx] 
  best_lambda  <- table_1se$lambda_1se[best_idx] 
  cvfit_best   <- per_alpha_cv[[as.character(best_alpha)]] 
  
  
  
  # coefficients and a handy predictor at exactly s = lambda.1se 
  
  coef_best <- as.matrix(coef(cvfit_best, s = best_lambda)) 
  predict_best <- function(newx, type = c("link", "response", "risk")) { 
    type <- match.arg(type) 
    predict(cvfit_best, newx = newx, s = best_lambda, type = type) 
  } 
  
  
  list( 
    cv_table_all   = cv_table_all,   # full CV grid across alphas & lambdas 
    per_alpha_1se  = table_1se,      # summary at lambda.1se for each alpha 
    best_alpha     = best_alpha, 
    best_lambda    = best_lambda,    # this is the chosen lambda.1se 
    cvfit_best     = cvfit_best,     # reuse this object; slice at s = best_lambda 
    coef_best      = coef_best, 
    predict_best   = predict_best 
  ) 
} 