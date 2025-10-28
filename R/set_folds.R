set_folds = function(data, nfold = 5, seed = 123, ngene=1000, method_trans_train="rank"){
  set.seed(seed)
  
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  strata = interaction(data$sampInfo$event,data$sampInfo$dataset,drop = FALSE)
  folds=caret::createFolds(strata,nfold,list=FALSE)
  
  data_train = list()
  data_test = list()
  for(i in 1:nfold){
    
    data_train[[i]] = list(ex = data$ex[,folds!=i,drop=FALSE], 
                           sampInfo = data$sampInfo[folds != i,,drop=FALSE],
                           featInfo = data$featInfo,
                           dataname = data$dataname,
                           samp_keeps = 1:sum(folds!=i))
    
    #fiter genes
    data_train[[i]] = preprocess_data(data_train[[i]], ngene=ngene, 
                                      method_trans_train=method_trans_train)

    
    data_test[[i]] = list(ex = data$ex[,folds==i,drop=FALSE],
                          sampInfo = data$sampInfo[folds == i,,drop=FALSE],
                          featInfo = data$featInfo,
                          dataname = data$dataname)
    
    #restrict to genes in training set
    data_test[[i]]$ex=data_test[[i]]$ex[rownames(data_train[[i]]$ex),]
    
    #normalize test data
    if(method_trans_train=='quant'){
      data_test[[i]]$ex = preprocessCore::normalize.quantiles.use.target(data_test[[i]]$ex,target = targets,keep.names=TRUE)
    }else if(method_trans_train=='rank'){
      data_test[[i]]$ex=apply(data_test[[i]]$ex,2,rank,ties.method="average")
    }

  }#end for loop
  
  return(list(data_train = data_train, data_test = data_test, folds=folds))
}
