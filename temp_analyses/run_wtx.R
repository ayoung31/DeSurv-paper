library(dplyr)
library(survival)
library(ggplot2)
library(parallel)
library(foreach)
source("R/load_data.R")
source("R/load_data_internal.R")
source("R/preprocess_data.R")
source("R/get_top_genes.R")
# purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)
setwd("/work/users/a/y/ayoung31/DeSurv-paper-dev")
source("code/misc.R")

args <- commandArgs(trailingOnly = TRUE)
ntop <- as.integer(args[[1]])
validation_dataset = args[[2]]
score = args[[3]]
method = args[[4]]

print(ntop)
print(validation_dataset)
print(score)
print(method)

ncore=11

model.params = list(
  dataset="TCGA_PAAD", # dataset we trained the model on
  version="TCGA_PAAD",
  ngene=1000,
  ninit=100,
  imaxit=6000,
  tol=1e-6,
  maxit=6000
)

select_model_v3_foreach_outer <- function(model.params, ntop, data, score_bin = TRUE,
                                          method = "none") {
  stopifnot(requireNamespace("Matrix", quietly = TRUE),
            requireNamespace("survival", quietly = TRUE),
            requireNamespace("foreach", quietly = TRUE))
  
  version <- model.params$version
  ngene   <- model.params$ngene
  ninit   <- model.params$ninit
  imaxit  <- model.params$imaxit
  tol     <- model.params$tol
  maxit   <- model.params$maxit
  
  ks <- 2:12
  
  res <- foreach::foreach(
    k = ks,
    .combine  = rbind,                      # combine per-k data.frames
    .inorder  = FALSE,                      # allow out-of-order completion
    .packages = c("Matrix", "survival"),    # load on workers
    .export = c("get_top_genes",
                "fit_val_model",
                "build_metrics_path",
                "build_model_filename",
                "build_model_path",
                "preprocess_data",
                "subset_with_zeros")
  ) %dopar% {
    
    message(sprintf(">>> Worker starting k = %d", k))
    
    metrics_path <- build_metrics_path(version, ngene, k, ninit, imaxit, tol, maxit)
    if (!file.exists(metrics_path)) {
      warning("Metrics file not found: ", metrics_path)
      return(data.frame(k = k, alpha = NA_real_, lambda = NA_real_, eta = NA_real_,
                        lambdaW = NA_real_, lambdaH = NA_real_, bic = NA_real_))
    }
    load(metrics_path)  # -> metrics
    nJ <- nrow(metrics)
    message(sprintf("[k=%d] Loaded metrics: %d rows", k, nJ))
    
    bics_k <- list(); bix <- 0L
    # pb <- txtProgressBar(min = 0, max = nJ, style = 3)
    
    for (j in seq_len(nJ)) {
      print(sprintf("k=%d j=%d",k,j))
      # setTxtProgressBar(pb, j)
      best_row <- metrics[j, ]
      
      model_filename <- build_model_filename(
        k, best_row$alpha, best_row$lambda, best_row$eta,
        best_row$lambdaW, best_row$lambdaH, ninit, imaxit, tol, maxit
      )
      model_path <- build_model_path(version, ngene, model_filename)
      if (!file.exists(model_path)) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_,
          c=NA_real_,k_eff=NA_real_
        )
        next
      }
      
      # Load model into a local env to avoid polluting the worker's global env
      e <- new.env(parent = emptyenv())
      load(model_path, envir = e)  # expects e$fit_cox
      if (!exists("fit_cox", envir = e, inherits = FALSE)) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_,
          c=NA_real_,k_eff = NA_real_
        )
        next
      }
      
      tops <- get_top_genes(e$fit_cox$W,ntop)
      
      if (is.null(tops)) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_,
          c=NA_real_, k_eff=NA_real_
        )
        next
      }
      
      # ---- Vectorized scoring (same as v3) ----
      W <- e$fit_cox$W
      
      data_filtered = preprocess_data(data,genes=rownames(W),method_trans_train=method)
      X=data_filtered$ex
      
      fit = fit_val_model(X,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops,W,score_bin=score_bin)
      
      bic_val <- stats::BIC(fit)
      c = fit$concordance[6]
      
      cols = which(unlist(lapply(tops, function(x) !is.null(x))))
      k_eff=length(cols)
      
      bix <- bix + 1L
      bics_k[[bix]] <- data.frame(
        k = k,
        alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
        lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH,
        bic = bic_val,c = c,k_eff=k_eff
      )
    } # j loop
    
    # close(pb)
    
    if (length(bics_k)) do.call(rbind, bics_k) else {
      data.frame(k = k, alpha = NA_real_, lambda = NA_real_, eta = NA_real_,
                 lambdaW = NA_real_, lambdaH = NA_real_, 
                 bic = NA_real_,c=NA_real_,k_eff=NA_real_)
    }
  } # foreach over k
  
  rownames(res) <- NULL
  res
}



#c("CPTAC","Dijk","Linehan","Moffitt_GEO_array","PACA_AU_array",
                        #"PACA_AU_seq","Puleo_array")

if(score=="bin"){
  score_bin=TRUE
}else{
  score_bin=FALSE
}
  

data=load_data(datasets=validation_dataset)

cl = parallel::makeCluster(ncore,outfile="")
doParallel::registerDoParallel(cl)
parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())

bics = select_model_v3_foreach_outer(model.params,ntop=ntop,data,
                                     score_bin=score_bin,method=method)

save(bics,file=paste0("wtx_bics_top",ntop,"_val_",validation_dataset,
                      "_",score,"_",method,".RData"))

stopCluster(cl)
