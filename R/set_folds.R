set_folds = function(data, nfold = 5, seed = 123, ngene = 1000, method_trans_train = "log2") {
  set.seed(seed)

  data$ex = data$ex[, data$samp_keeps, drop = FALSE]
  data$sampInfo = data$sampInfo[data$samp_keeps, , drop = FALSE]
  strata = interaction(data$sampInfo$event, data$sampInfo$dataset, drop = FALSE)
  folds = caret::createFolds(strata, nfold, list = FALSE)

  data_train = vector("list", nfold)
  data_test = vector("list", nfold)

  for (i in seq_len(nfold)) {
    keep_train <- folds != i
    keep_test <- folds == i

    train_proc <- DeSurv::preprocess_data(
      X = data$ex[, keep_train, drop = FALSE],
      y = data$sampInfo$time[keep_train],
      d = data$sampInfo$event[keep_train],
      dataset = data$sampInfo$dataset[keep_train],
      samp_keeps = seq_len(sum(keep_train)),
      ngene = ngene,
      method_trans_train = method_trans_train,
      verbose = FALSE
    )
    train_proc$dataname <- data$dataname
    data_train[[i]] <- train_proc

    test_proc <- DeSurv::preprocess_data(
      X = data$ex[, keep_test, drop = FALSE],
      y = data$sampInfo$time[keep_test],
      d = data$sampInfo$event[keep_test],
      dataset = data$sampInfo$dataset[keep_test],
      samp_keeps = seq_len(sum(keep_test)),
      genes = train_proc$featInfo,
      method_trans_train = method_trans_train,
      verbose = FALSE
    )
    test_proc$dataname <- data$dataname
    data_test[[i]] <- test_proc
  }

  list(data_train = data_train, data_test = data_test, folds = folds)
}
