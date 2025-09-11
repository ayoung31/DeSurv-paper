set_folds = function(data, nfold = 5){
  
  folds=caret::createFolds(data$sampInfo$time,k=nfold,list=FALSE)
  
  data_train = list()
  data_test = list()
  for(i in 1:nfold){

    data_train[[i]] = list(ex = data$ex[,folds!=i,drop=FALSE], 
                           sampInfo = data$sampInfo[folds != i,,drop=FALSE],
                           featInfo = data$featInfo,
                           dataname = data$dataname)

    data_test[[i]] = list(ex = data$ex[,folds==i,drop=FALSE],
                          sampInfo = data$sampInfo[folds == i,,drop=FALSE],
                          featInfo = data$featInfo,
                          dataname = data$dataname)

  }#end for loop
  
  return(list(data_train = data_train, data_test = data_test))
}
