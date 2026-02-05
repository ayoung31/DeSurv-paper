compute_metrics_bundle = function(data_folds,bundle,fold){
  
  X = data_folds$data_train[[fold]]$ex
  y = data_folds$data_train[[fold]]$sampInfo$time
  d = data_folds$data_train[[fold]]$sampInfo$event
  
  Xtest = data_folds$data_test[[fold]]$ex
  ytest = data_folds$data_test[[fold]]$sampInfo$time
  dtest = data_folds$data_test[[fold]]$sampInfo$event
  
  mets_tr=list()
  mets_te=list()
  j=1
  for(i in 1:length(bundle)){
    b = bundle[[i]]
    meta = b$meta
    for(a in names(b$fits)){
      alpha = as.numeric(a)
      fit = b$fits[[a]]
      
      mets_tr[[j]] = tryCatch({
        m = compute_metrics(fit,X,y,d,test=FALSE)
        m$fold = fold
        m$alpha = alpha
        m$k = meta$k
        m$lambda = meta$lambda
        m$eta = meta$eta
        m$lambdaW = meta$lambdaW
        m$lambdaH = meta$lambdaH
        m$seed = meta$seed
        m
      },error=function(e) NULL)
      
      
      mets_te[[j]] = tryCatch({
        m = compute_metrics(fit,Xtest,ytest,dtest,test=TRUE)
        m$fold = fold
        m$alpha = alpha
        m$k = meta$k
        m$lambda = meta$lambda
        m$eta = meta$eta
        m$lambdaW = meta$lambdaW
        m$lambdaH = meta$lambdaH
        m$seed = meta$seed
        m
      },error=function(e) NULL)
      
      j=j+1
    }
  }
  
  mets_train = dplyr::bind_rows(mets_tr)
  mets_test = dplyr::bind_rows(mets_te)
  
  setNames(list(list(mets_train=mets_train, mets_test=mets_test)), 
           paste0(meta$k,meta$lambda,meta$eta,
                  meta$lambdaW,meta$lambdaH,fold))
}