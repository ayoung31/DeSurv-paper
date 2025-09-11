load("wtx_bics_top25_val_cptac_rank.RData")
load("data/derv/cmbSubtypes_formatted.RData")

temp = bics %>% group_by(k) %>% slice_min(order_by = bic,n=1,with_ties = FALSE)

plot(temp$k,temp$bic,type='l')



best = temp[temp$k==5,]

bics_a0 = bics %>% filter(alpha==0)
best_a0=bics_a0[which.min(bics_a0$bic),]

ntop = 25

model.params = list(
  dataset="TCGA_PAAD", # dataset we trained the model on
  version="TCGA_PAAD",
  ngene=1000,
  ninit=100,
  imaxit=6000,
  tol=1e-6,
  maxit=6000
)


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


model_filename = build_model_filename(best$k,best$alpha,best$lambda,best$eta,
                                      best$lambdaW,best$lambdaH,model.params$ninit,
                                      model.params$imaxit,model.params$tol,model.params$maxit)

model_path = build_model_path(model.params$version,model.params$ngene,model_filename)

load(model_path)

results <- list(ntop = ntop, fit_cox = fit_cox, model.params = model.params)

tops =get_top_genes_v2(results)
results$tops = tops$tops
create_table_v2(results,top_genes,"DECODER",colors,save=FALSE)

validation_datasets = "Dijk"#c("CPTAC","Dijk","Linehan","Moffitt_GEO_array","PACA_AU_array",
#"PACA_AU_seq","Puleo_array")

data=load_data(datasets=validation_datasets,save=TRUE,replace=FALSE)
dim(data$ex)

X=data$ex
W=fit_cox$W

g_common <- intersect(rownames(W), rownames(X))

Wc <- W[g_common, , drop = FALSE]
Xc <- X[g_common, , drop = FALSE]

p <- ncol(results$tops)
idx_list <- lapply(seq_len(p), function(i) {
  m <- match(results$tops[, i], g_common)
  m[!is.na(m)]
})
rows <- unlist(idx_list, use.names = FALSE)
cols <- rep.int(seq_len(p), vapply(idx_list, length, integer(1)))


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
