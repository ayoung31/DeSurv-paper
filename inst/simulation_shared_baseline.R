purrr::walk(list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$"), source)
library(NMF)
library(DeSurv)

sim <- simulate_desurv_data_shared_baseline()  
X <- sim$counts              # or normalized version if deSurv expects normalized
surv_all <- sim$surv         # data.frame(patient,time,status)
rownames(surv_all) <- sim$surv$patient

X = apply(X,2,rank)

set.seed(1)
idx <- split_train_test(ncol(X), train_frac = 0.7)

X_train <- X[, idx$train, drop = FALSE]
X_test  <- X[, idx$test,  drop = FALSE]

surv_all$time <- pmax(surv_all$time, 1e-3)
surv_train <- surv_all[idx$train, , drop = FALSE]
surv_test  <- surv_all[idx$test,  , drop = FALSE]

nmf_fit <- NMF::nmf(X,4,method="lee",nrun=20)
nmf_fit2 = list(W=nmf_fit@fit@W,H=nmf_fit@fit@H)
W_nmf   <- nmf_fit2$W  # G x K

Z_train_nmf <- project_scores(X_train, W_nmf)  # Ntrain x K
Z_test_nmf  <- project_scores(X_test,  W_nmf)  # Ntest x K

cox_nmf <- fit_cox_on_Z(Z_train_nmf, surv_train)

# Get risk predictions on test using the same Cox model
# predict.coxph() with newdata:
predict_test_nmf <- {
  df_test <- data.frame(time = surv_test$time,
                        status = surv_test$status,
                        Z_test_nmf)
  colnames(df_test)[-(1:2)] <- paste0("X", seq_len(ncol(Z_test_nmf)))
  predict(cox_nmf$fit, newdata = df_test, type = "lp")
}

cindex_nmf <- compute_cindex(predict_test_nmf, surv_test)

desurv_fit <- DeSurv::desurv_fit(
  X = X_train,
  y = surv_train$time,
  d = surv_train$status,
  k = 4,
  alpha = 0.35,
  lambda = 0.1,
  nu = 0.9,
  lambdaW = 1e-3,
  lambdaH = 1e-3,
  tol = 1e-5,
  imaxit = 50,
  ninit = 5,
  maxit = 6000,
  verbose = TRUE
)

W_desurv    <- desurv_fit$W        # G x K
beta_desurv <- as.numeric(desurv_fit$beta)     # length K (optional but nice)
cvwrapr::getCindex(t(X_train)%*%W_desurv%*%beta_desurv,Surv(surv_train$time,surv_train$status))
Z_train_desurv <- project_scores(X_train, W_desurv)  # Ntrain x K
Z_test_desurv  <- project_scores(X_test,  W_desurv)  # Ntest x K

cox_desurv <- fit_cox_on_Z(Z_train_desurv, surv_train)

predict_test_desurv <- {
  df_test <- data.frame(time = surv_test$time,
                        status = surv_test$status,
                        Z_test_desurv)
  colnames(df_test)[-(1:2)] <- paste0("X", seq_len(ncol(Z_test_desurv)))
  predict(cox_desurv$fit, newdata = df_test, type = "lp")
}

cindex_desurv <- compute_cindex(predict_test_desurv, surv_test)

risk_test_desurv_direct <- as.vector(Z_test_desurv %*% beta_desurv)
cindex_desurv_direct <- compute_cindex(risk_test_desurv_direct, surv_test)
