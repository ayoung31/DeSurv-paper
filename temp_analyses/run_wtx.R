library(dplyr)
library(survival)
library(ggplot2)
library(parallel)
library(foreach)
library(coxNMF)
purrr::walk(list.files("/work/users/a/y/ayoung31/DeSurv-paper/R", full.names = TRUE, pattern = "[.]R$"), source)
    
setwd("/work/users/a/y/ayoung31/DeSurv-paper-dev")


args <- commandArgs(trailingOnly = TRUE)
ntop <- as.integer(args[[1]])

ncore=16

model.params = list(
  dataset="TCGA_PAAD", # dataset we trained the model on
  version="TCGA_PAAD",
  ngene=1000,
  ninit=100,
  imaxit=6000,
  tol=1e-6,
  maxit=6000
)

select_model_v3_foreach_outer <- function(model.params, ntop, data, 
                                          score_bin = FALSE, rank = FALSE) {
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
    .packages = c("Matrix", "survival", "dplyr", "coxNMF")     # load on workers
  ) %dopar% {
    
    purrr::walk(list.files("/work/users/a/y/ayoung31/DeSurv-paper/R", full.names = TRUE, pattern = "[.]R$"), source)
    
    
    message(sprintf(">>> Worker starting k = %d", k))
    
    build_metrics_path <- function(version, ngene, k, ninit, imaxit, tol, maxit) {
      file.path('results', version, 'full', paste0('ngene', ngene), 'summary',
                paste0('metrics_ngene', ngene, '_k=', k,
                       '_ninit=', ninit, '_imaxit=', imaxit,
                       '_tol=', tol, '_maxit=', maxit, '.RData'))
    }
    
    build_model_filename <- function(k, alpha, lambda, eta, lambdaW, lambdaH,
                                     ninit, imaxit, tol, maxit) {
      paste0("k=", k, "_alpha", alpha,
             "_lambda", lambda, "_eta", eta,
             "_lambdaW", lambdaW, "_lambdaH", lambdaH,
             "_full_ninit", ninit, "_imaxit", imaxit,
             "_tol", tol, "_maxit", maxit, ".RData")
    }
    
    #' Generate full path to the raw model
    build_model_path <- function(version, ngene, filename) {
      file.path('results', version, 'full', paste0('ngene', ngene), 'raw', filename)
    }
    
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
      # print(sprintf("k=%d j=%d",k,j))
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
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_
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
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_
        )
        next
      }
      
      tops <- get_top_genes(e$fit_cox$W,ntop)
      
      if (is.null(tops)) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_
        )
        next
      }
      
      # ---- Vectorized scoring (same as v3) ----
      W <- e$fit_cox$W
      X <- data$ex
      
      
      g_common <- intersect(rownames(W), rownames(X))
      if (length(g_common) == 0L) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_
        )
        next
      }
      
      Wc <- W[g_common, , drop = FALSE]
      Xc <- X[g_common, , drop = FALSE]
      
      p <- ncol(tops)
      idx_list <- lapply(seq_len(p), function(i) {
        m <- match(tops[, i], g_common)
        m[!is.na(m)]
      })
      rows <- unlist(idx_list, use.names = FALSE)
      cols <- rep.int(seq_len(p), vapply(idx_list, length, integer(1)))
      if (length(rows) == 0L) {
        bix <- bix + 1L
        bics_k[[bix]] <- data.frame(
          k = k,
          alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
          lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH, bic = NA_real_
        )
        next
      }
      
      M <- Matrix::sparseMatrix(
        i = rows, j = cols, x = 1,
        dims = dim(Wc), dimnames = dimnames(Wc)
      )
      
      if (score_bin) {
        W_eff <- Matrix::drop0(Matrix::Matrix((Wc > 0) * (M != 0), sparse = TRUE))
      } else {
        W_eff <- Matrix::drop0(Matrix::Matrix(Wc, sparse = TRUE) * M)
      }
      
      scores <- as.matrix(Matrix::crossprod(as.matrix(Xc), W_eff))  # samples x factors
      
      score_data <- data.frame(scores,
                               time = data$sampInfo$time,
                               event = data$sampInfo$event,
                               check.names = FALSE)
      
      fit <- survival::coxph(
        survival::Surv(time, event) ~ .,
        data = score_data,
        model = FALSE, x = FALSE, y = FALSE
      )
      bic_val <- stats::BIC(fit)
      
      bix <- bix + 1L
      bics_k[[bix]] <- data.frame(
        k = k,
        alpha = best_row$alpha, lambda = best_row$lambda, eta = best_row$eta,
        lambdaW = best_row$lambdaW, lambdaH = best_row$lambdaH,
        bic = bic_val
      )
    } # j loop
    
    # close(pb)
    
    if (length(bics_k)) do.call(rbind, bics_k) else {
      data.frame(k = k, alpha = NA_real_, lambda = NA_real_, eta = NA_real_,
                 lambdaW = NA_real_, lambdaH = NA_real_, bic = NA_real_)
    }
  } # foreach over k
  
  rownames(res) <- NULL
  res
}



validation_datasets = "CPTAC"#c("CPTAC","Dijk","Linehan","Moffitt_GEO_array","PACA_AU_array",
                        #"PACA_AU_seq","Puleo_array")



cl = parallel::makeCluster(ncore,outfile="")
doParallel::registerDoParallel(cl)
parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())

data=load_data(datasets=validation_datasets)
data = preprocess_data(data,method_trans_train="quant")

bics = select_model_v3_foreach_outer(model.params,ntop=ntop,data,
                                     score_bin=FALSE, rank=FALSE)


save(bics,file=paste0("wtx_bics_top",ntop,"_val_cptac_quant_nobin.RData"))


bics = select_model_v3_foreach_outer(model.params,ntop=ntop,data,
                                     score_bin=TRUE, rank=FALSE)


save(bics,file=paste0("wtx_bics_top",ntop,"_val_cptac_quant_bin.RData"))

stopCluster(cl)
