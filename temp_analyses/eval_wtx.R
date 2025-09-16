library(targets)
tar_load_globals()

load("../DeSurv-paper-dev/data/derv/cmbSubtypes_formatted.RData")
library(dplyr)
library(ggplot2)

ntop = 25
dataset = "CPTAC"
score="bin"
method="rank"

load(paste0("../DeSurv-paper-dev/wtx_bics_top",ntop,"_val_",dataset,"_",score,"_",method,".RData"))


temp = bics%>%group_by(alpha,k) %>% 
  slice_min(order_by = bic,n=1,with_ties=FALSE) %>%
  ungroup()

ggplot(temp,aes(x=k,y=alpha,fill=bic))+
  geom_tile()+
  scale_fill_gradient(low = "red", high = "white")


best = temp %>% slice_min(order_by = bic,n=1,with_ties=FALSE)

bics_a0 = bics%>% filter(alpha==0)
best_a0=bics_a0[which.min(bics_a0$bic),]



test_datasets = "Puleo_array"#c("Dijk","Linehan","Moffitt_GEO_array","PACA_AU_array",
#"PACA_AU_seq","Puleo_array")
data=load_data(datasets=test_datasets)

## best model
fit_cox=load_model(best)
tops =get_top_genes(fit_cox$W,ntop)
W=fit_cox$W

data_filtered = preprocess_data(data,genes=rownames(W))
X=data_filtered$ex

fit=fit_val_model(X,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops,W)
bic_val <- stats::BIC(fit)
summary(fit)
cols = which(unlist(lapply(tops, function(x) !is.null(x))))
tops = as.data.frame(tops[cols])
create_table(as.data.frame(tops),top_genes,"DECODER",colors)

# alpha 0 best model
fit_cox_a0=load_model(best_a0)
tops_a0 =get_top_genes(fit_cox_a0$W,ntop)
W_a0=fit_cox_a0$W

data_filtered_a0 = preprocess_data(data,genes=rownames(W_a0))
X_a0=data_filtered_a0$ex

fit_a0=fit_val_model(X_a0,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops_a0,W_a0)
bic_val_a0 <- stats::BIC(fit_a0)
summary(fit_a0)
cols = which(unlist(lapply(tops_a0, function(x) !is.null(x))))
tops_a0 = as.data.frame(tops_a0[cols])
create_table(tops_a0,top_genes,"DECODER",colors,)


# std NMF
tar_load(data_filtered)
nmf_fit=NMF::nmf(as.matrix(data_filtered$ex),rank=5,method="lee",nrun=50)
tops_std =get_top_genes(nmf_fit@fit@W,ntop)
W_std=nmf_fit@fit@W

data_filtered_std = preprocess_data(data,genes=rownames(W_std))
X_std=data_filtered_std$ex

fit_std=fit_val_model(X_std,data_filtered_std$sampInfo$time,data_filtered_std$sampInfo$event,tops_std,W_std)
bic_val_std <- stats::BIC(fit_std)
summary(fit_std)
create_table(tops_std,top_genes,"DECODER",colors,)
