# X_train, X_val: genes x samples (rows = genes), numeric
# Returns: list of per-gene maps from training and transformed matrices on [0,1]
fit_rank_maps <- function(X_train) {
  apply(X_train, 1, function(v) {
    o <- order(v); sv <- v[o]; n <- length(sv)
    rle_v <- rle(sv)
    ends   <- cumsum(rle_v$lengths)
    starts <- ends - rle_v$lengths + 1
    midr   <- (starts + ends) / 2              # midrank indices
    list(vals = rle_v$values, midr = midr / n, n = n)
  })
}

apply_rank_maps <- function(X, maps) {
  stopifnot(nrow(X) == length(maps))
  Y <- X
  for (g in seq_len(nrow(X))) {
    m <- maps[[g]]
    # findInterval gives the index of the last training value <= x
    idx <- findInterval(X[g, ], m$vals, rightmost.closed = TRUE)
    idx[idx < 1] <- 1                          # below min -> lowest bin
    idx[idx > length(m$midr)] <- length(m$midr)# above max -> highest bin
    Y[g, ] <- m$midr[idx]                      # map to training midrank in [0,1]
  }
  dimnames(Y) <- dimnames(X)
  Y
}

## usage

make_stratified_folds <- function(time, status, K = 5, seed = 1) {
  set.seed(seed)
  # bin event times into quantiles (events only)
  q <- quantile(time[status == 1], probs = seq(0, 1, length.out = 6), na.rm = TRUE)
  time_bin <- cut(time, breaks = unique(q), include.lowest = TRUE)
  # label censored separately so their distribution is also balanced
  strat <- interaction(factor(status, levels = c(0,1), labels = c("C","E")), addNA(time_bin))
  # create folds by stratified sampling
  levs <- levels(strat)
  idx <- split(seq_along(time), strat)
  folds <- rep(NA_integer_, length(time))
  for (lv in names(idx)) {
    ids <- idx[[lv]]
    if (length(ids) == 0) next
    krep <- rep(1:K, length.out = length(ids))
    krep <- sample(krep)  # shuffle
    folds[ids] <- krep
  }
  folds
}


set_folds = function(data, nfold = 5, seed = 123, ngene=1000, method_trans_train="rank"){
  set.seed(seed)
  
  data$ex = data$ex[,data$samp_keeps]
  data$sampInfo = data$sampInfo[data$samp_keeps,]
  strata = interaction(data$sampInfo$event,data$sampInfo$dataset,drop = FALSE)
  folds=caret::createFolds(strata,nfold,list=FALSE)
  # make_stratified_folds(data$sampInfo$time,data$sampInfo$event,nfold,seed=seed)
  
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
    # #normalize
    # if(METHOD_TRANS_TRAIN=="quant"){
      # rns = rownames(data_train[[i]]$ex)
      # targets=preprocessCore::normalize.quantiles.determine.target(data_train[[i]]$ex)
      # data_train[[i]]$ex = preprocessCore::normalize.quantiles.use.target(data_train[[i]]$ex,
      #                                                                     target=targets)
      # rownames(data_train[[i]]$ex)=rns
      
    # }else if(METHOD_TRANS_TRAIN=="rank"){
    #   data_train[[i]]$ex=apply(data_train[[i]]$ex,2,rank,ties.method="average")
    # }
    
    
    
    data_test[[i]] = list(ex = data$ex[,folds==i,drop=FALSE],
                          sampInfo = data$sampInfo[folds == i,,drop=FALSE],
                          featInfo = data$featInfo,
                          dataname = data$dataname)
    
    data_test[[i]]$ex=data_test[[i]]$ex[rownames(data_train[[i]]$ex),]
    if(method_trans_train=='quant'){
      # rns = rownames(data_test[[i]]$ex)
      # data_test[[i]]$ex = preprocessCore::normalize.quantiles.use.target(data_test[[i]]$ex,target = targets)
      # rownames(data_test[[i]]$ex) = rns
      stop("need to return quantiles of training data from preprocess_data function")
    }else if(method_trans_train=='rank'){
      data_test[[i]]$ex=apply(data_test[[i]]$ex,2,rank,ties.method="average")
    }
    # 
    
    
    # maps       <- fit_rank_maps(data_train[[i]]$ex)
    # data_train[[i]]$ex  <- apply_rank_maps(data_train[[i]]$ex, maps)   # train-in, train-out (idempotent up to ties)
    # data_test[[i]]$ex    <- apply_rank_maps(data_test[[i]]$ex,   maps)   # apply train map to validation (no leakage)
    
    
  }#end for loop
  
  return(list(data_train = data_train, data_test = data_test, folds=folds))
}
