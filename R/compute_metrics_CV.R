compute_metrics_CV = function(path,data_folds,ntop,score_bin=FALSE){
  
  b <- readRDS(path)
  stopifnot(is.list(b$fits), length(b$fits) > 0)
  pars <- b$meta
  
  data_test = data_folds$data_test[[pars$fold]]
  X = data_test$ex
  y = data_test$sampInfo$time
  delta = data_test$sampInfo$event
  res <- vector("list", length(b$fits))
  j = 1
  for (nm in sort(names(b$fits))){
    a <- as.numeric(nm)
    fit <- b$fits[[nm]]
    W <- fit$W
    if(!is.null(W) & sum(W)>0){
      tops = get_top_genes(W = W, ntop = ntop)
      score_data = compute_scores(tops,W,X,y,delta,score_bin)
      survfit = survival::coxph(Surv(time,event)~., data=score_data)
      
      bic = stats::BIC(survfit)
      
      
      
    }else{
      bic = NA
    }
    
    #create a dataframe of params and bic results
    res_df = pars
    res_df$alpha = a
    res_df$bic = bic
    res[[j]] = res_df
    
    j = j+1
  }
  
  metrics = dplyr::bind_rows(res)
  
  return(metrics)
}