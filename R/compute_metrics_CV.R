compute_metrics_CV = function(path_out,path,data_folds,ntop,score_bin=FALSE){
  
  b <- readRDS(path)
  stopifnot(is.list(b$fits), length(b$fits) > 0)
  pars <- b$meta
  
  data_test = data_folds$data_test[[pars$fold]]
  X = data_test$ex
  y = data_test$sampInfo$time
  delta = data_test$sampInfo$event
  res <- vector("list", length(b$fits))
  res_df = as.data.frame(pars)
  res_df$alpha = res_df$alpha_vec
  res_df$alpha_vec = NULL
  
  j = 1
  for (nm in sort(names(b$fits))){
    a <- as.numeric(nm)
    fit <- b$fits[[nm]]
    W <- fit$W
    if(!is.null(W) & sum(W)>0){
      tops = get_top_genes(W = W, ntop = ntop)
      score_data = compute_scores(tops,W,X,y,delta,score_bin)
      bic = tryCatch(
        {
          survfit = survival::coxph(Surv(time,event)~., data=score_data)
          stats::BIC(survfit)
        },
        error = function(e){
          NA
        }
      )
      
      
      
      
    }else{
      bic = NA
    }
    
    #create a dataframe of params and bic results
    cur_df = res_df[j,]
    cur_df$bic = bic
    res[[j]] = cur_df
    
    j = j+1
  }
  
  metrics = dplyr::bind_rows(res)
  
  tmp <- paste0(path_out, ".tmp")
  saveRDS(metrics, tmp, compress = "xz")
  file.rename(tmp, path_out)

  
  return(path_out)
}
