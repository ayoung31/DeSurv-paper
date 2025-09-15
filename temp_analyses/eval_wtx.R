library(targets)
tar_load_globals()

load("../DeSurv-paper-dev/data/derv/cmbSubtypes_formatted.RData")
library(dplyr)
library(ggplot2)

ntop = 25

load("../DeSurv-paper-dev/wtx_bics_top25_val_cptac_quant_bin.RData")
bics_quant_bin=bics

temp_quant_bin = bics_quant_bin %>% group_by(alpha,k) %>% 
  slice_min(order_by = bic,n=1,with_ties=FALSE) %>%
  ungroup()

ggplot(temp_quant_bin,aes(x=k,y=alpha,fill=bic))+
  geom_tile()


best = temp_quant_bin %>% filter(k==4)%>% slice_min(order_by = bic,n=1,with_ties=FALSE)

bics_a0_quant_bin = bics_quant_bin %>% filter(alpha==0)
best_a0=bics_a0_quant_bin[which.min(bics_a0_quant_bin$bic),]



validation_datasets = "TCGA_PAAD"#c("Dijk","Linehan","Moffitt_GEO_array","PACA_AU_array",
#"PACA_AU_seq","Puleo_array")
data=load_data(datasets=validation_datasets)
table(rownames(data$ex) %in% rownames(W))

## best model at k=4
fit_cox=load_model(best)
tops =get_top_genes(fit_cox$W,ntop)
W=fit_cox$W

data_filtered = preprocess_data(data,genes=rownames(W))
X=data_filtered$ex

fit=fit_val_model(X,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops,W)
bic_val <- stats::BIC(fit)
summary(fit)
# create_table(tops,top_genes,"DECODER",colors)

# alpha 0 at k=4
fit_cox_a0=load_model(best_a0)
tops_a0 =get_top_genes(fit_cox_a0$W,ntop)
W_a0=fit_cox_a0$W

data_filtered_a0 = preprocess_data(data,genes=rownames(W_a0))
X_a0=data_filtered_a0$ex

fit_a0=fit_val_model(X_a0,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops_a0,W_a0)
bic_val_a0 <- stats::BIC(fit_a0)
summary(fit_a0)
# create_table(tops_a0,top_genes,"DECODER",colors,)


# std NMF at k=4
# tar_load(data_filtered)
# nmf_fit=NMF::nmf(as.matrix(data_filtered$ex),rank=4,method="lee",nrun=50)
tops_std =get_top_genes(nmf_fit@fit@W,ntop)
W_std=nmf_fit@fit@W

data_filtered_std = preprocess_data(data,genes=rownames(W_std))
X_std=data_filtered_std$ex

fit_std=fit_val_model(X_std,data_filtered$sampInfo$time,data_filtered$sampInfo$event,tops_std,W_std)
bic_val_std <- stats::BIC(fit_std)
summary(fit_std)
# create_table(tops_std,top_genes,"DECODER",colors,)
