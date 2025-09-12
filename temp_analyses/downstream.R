library(targets)
library(NMF)
library(coxNMF)

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


purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

model.params = list(
  dataset="TCGA_PAAD", # dataset we trained the model on
  version="TCGA_PAAD",
  ngene=1000,
  ninit=100,
  imaxit=6000,
  tol=1e-6,
  maxit=6000
)

version <- model.params$version
ngene   <- model.params$ngene
ninit   <- model.params$ninit
imaxit  <- model.params$imaxit
tol     <- model.params$tol
maxit   <- model.params$maxit

# load processed data
tar_load(data_filtered)

# run final model
load("/work/users/a/y/ayoung31/DeSurv-paper-dev/wtx_bics_top25_val_cptac_quant_nobin.RData")
best_row=bics %>% slice_min(order_by=bic,n=1,with_ties=FALSE)
model_filename <- build_model_filename(
  best_row$k, best_row$alpha, best_row$lambda, best_row$eta,
  best_row$lambdaW, best_row$lambdaH, ninit, imaxit, tol, maxit
)
model_path <- build_model_path(version, ngene, model_filename)
load(paste0("/work/users/a/y/ayoung31/DeSurv-paper-dev/",model_path))

## now standard NMF
# select k from visual inspection
# fit_std = NMF::nmf(data_filtered$ex,rank=2:12,method="lee", seed=123, nrun = 30, .options="v4")
# dir.create("results/TCGA_PAAD/std_NMF/")
# save(fit_std,file="results/TCGA_PAAD/std_NMF/select_k_fits.RData")
load("results/TCGA_PAAD/std_NMF/select_k_fits.RData")
dir.create("figures/TCGA_PAAD/std_NMF/",recursive=TRUE)
png("figures/TCGA_PAAD/std_NMF/select_k.png")
plot(fit_std)
dev.off()
# run with selected k for 100 initializations
# fit_std_k5 = NMF::nmf(data_filtered$ex, rank=5, seed=123, method = "lee", nrun=100)
# save(fit_std_k5,file="results/TCGA_PAAD/std_NMF/fit_k=5.RData")
load("results/TCGA_PAAD/std_NMF/fit_k=5.RData")
# load and eval on validation datasets
score_bin=FALSE
val_datasets = c("Dijk","Linehan","Moffitt_GEO_array",
                   "PACA_AU_array","PACA_AU_seq","Puleo_array")
bic_ind = numeric(length=length(val_datasets))
bic_ind_std = numeric(length=length(val_datasets))
scores = list()
scores_std = list()
for(i in 1:length(val_datasets)){
  print(val_datasets[i])
  raw = load_data(val_datasets)
  data = preprocess_data(data = raw,
                         ngene = 1000,
                         method_trans_train = "rank")
  
  # DeSurv
  tops = get_top_genes(W = fit_cox$W, ntop = 25)
  score_data = compute_scores(tops, fit_cox$W,data$ex,data$sampInfo$time,data$sampInfo$event,score_bin)
  survfit = survival::coxph(Surv(time,event)~., data=score_data)
  bic_ind[i] = stats::BIC(survfit)
  score_data$dataset = val_datasets[i]
  scores[[i]] = score_data
  
  #std NMF
  tops_std = get_top_genes(W = fit_std_k5@fit@W, ntop = 25)
  score_data_std = compute_scores(tops_std, fit_std_k5@fit@W, data$ex, data$sampInfo$time, data$sampInfo$event, score_bin)
  survfit_std = survival::coxph(Surv(time,event)~., data=score_data_std)
  bic_ind_std[i] = stats::BIC(survfit_std)
  score_data_std$dataset = val_datasets[i]
  scores_std[[i]] = score_data_std
  
  rm(raw)
  rm(data)
  rm(score_data)
  rm(score_data_std)
}

score_data = dplyr::bind_rows(scores)
score_data_std = dplyr::bind_rows(scores_std)
dir.create("results/comps/")
save(bic_ind,bic_ind_std,score_data,score_data_std,file=
       "results/comps/bics_and_scores_rank.RData")

#### survival curves deSurv

fit_surv = survival::coxph(Surv(time,event)~.+strata(dataset),data=score_data)
summary(fit_surv)
cox.zph(fit_surv)

#plot survival by quantiles of scores
score_data$lp = predict(fit_surv,type="lp")
score_data = score_data %>% 
  mutate(lp_group = cut(lp,
                        breaks = quantile(lp, probs = seq(0,1,by=.25), na.rm=TRUE),
                        include.lowest = TRUE,
                        labels = paste0("Q",1:4)))

km_survq = survfit(Surv(time, event) ~lp_group, data=score_data)

ggsurvplot(
  km_survq,
  data = df,
  risk.table = TRUE,
  pval = TRUE,
  legend.title = "Risk (Cox LP) quantile",
  xlab = "Time",
  ylab = "Survival probability"
)

# should also break each factors scores into quantiles and plot by this


# get the basal/classical factors
bas_fac = fill
clas_fac = fill

#cluster on these factors?


#####survival curves std nmf





