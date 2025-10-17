data_val_filtered = lapply(data_val,preprocess_data_val,ngene=2000,method_trans_train="rank")
dataset=1
data_val[[6]]$ex
X=gene_filter_test(data_val[[6]]$ex,1000)
dataname = clusters_best_1000[[dataset]]$data$dataname
dataname
X = data_val[[dataname]]$ex
plot_heatmap.clustering_result(clusters_best_1000[[dataset]]$clus_res,k=3,Xtemp=data_val_filtered$Puleo_array$ex)
# 
# gcommon = intersect(rownames(data_val_filtered$Puleo_array$ex),unlist(tops_best_50_1000$top_genes[,c(1,3)]))
# X=t(apply(data_val_filtered$Puleo_array$ex[gcommon,],1,scale))
# X=scale(X[gcommon,])
# X[X>2]=2
# X[X< -2]=-2
# pheatmap::pheatmap(X)

# rownames(X)
# hgnc_fix=HGNChelper::checkGeneSymbols(rownames(expr_clean))
# genes_fixed <- hgnc_fix$Suggested.Symbol
# genes_fixed[is.na(genes_fixed)] <- hgnc_fix$x[is.na(genes_fixed)]  # fallback to original if no suggestion
# rownames(expr_clean) <- genes_fixed
# matches <- grep(paste(unlist(tops_best_1000$top_genes[,c(1,3)]), collapse = "|"), rownames(expr_clean), value = TRUE)
# matches
# 

tar_load(data_val_filtered_1000)
tar_load(tops_best_75_1000)
data=data_val_filtered_1000[[6]]
run_clustering(tops_best_75_1000$top_genes,data,top_genes,colors,
               facs=c(1,3),plot=FALSE,maxKcol = 5, maxKrow = 5)

tar_load(clusters_best_75_1000)
clusters_best_1000=clusters_best_75_1000
dataset=6

nclus_gene = rep(4,6)#c(4,4,4,4,4,4)
# nclus_gene = rep(2,6)

groups=list()
genes_per=list()
for(dataset in 1:6){
  dataname=clusters_best_1000[[dataset]]$data$dataname
  X=clusters_best_1000[[dataset]]$data$ex
  plot_heatmap.clustering_result(clusters_best_1000[[dataset]]$clus_res,k=4,Xtemp=X)
  
  #create dataframe with cols subject,gene,sub clus, gene clus,expression
  genes = rownames(X)
  expr = as.data.frame(X) %>% pivot_longer(everything(),names_to="subject",values_to ="expression")
  expr$gene = rep(genes,each=ncol(X))
  col_class = clusters_best_1000[[dataset]]$clus_res$clusCol[[4]]$consensusClass
  row_class = clusters_best_1000[[dataset]]$clus_res$clusRow[[nclus_gene[dataset]]]$consensusClass
  row_lab = data.frame(gene_clus = row_class,gene=names(row_class))
  col_lab = data.frame(samp_clus = col_class,subject=names(col_class))
  expr = expr %>% left_join(row_lab) %>% left_join(col_lab)
  expr=expr[!is.na(expr$gene_clus),]
  
  for(g in 1:nclus_gene[dataset]){
    genes_per[[dataname]][[g]] = row_lab$gene[row_lab$gene_clus==g]
  }
  # genes_per[[dataname]] = row_lab %>% group_by(gene_clus) %>% summarise(genes = paste(gene,collapse=", "))
  # genes_per[[dataset]]$dataset = dataname
  temp = expr %>% 
    group_by(gene_clus,samp_clus) %>% 
    summarise(mean_expr = mean(expression)) %>%
    group_by(samp_clus) %>%
    summarise(gene_clus=gene_clus[which.max(mean_expr)])
    # pivot_wider(names_from = gene_clus, values_from = mean_expr)%>%
    # ungroup()%>%
    # select(-samp_clus)
  # assign = clue::solve_LSAP(as.matrix(temp),maximum=TRUE)
  groups[[dataset]] = temp#data.frame(samp_clus=rownames(temp),gene_clus=colnames(temp)[assign])
}

l=1
inter=list()
for(dataset1 in 1){
  for(dataset2 in (dataset1 + 1):6){
    dataname1 = clusters_best_1000[[dataset1]]$data$dataname
    dataname2 = clusters_best_1000[[dataset2]]$data$dataname
    mat = matrix(nrow=nclus_gene[dataset1],ncol=nclus_gene[dataset2])
    for(i in 1:nclus_gene[dataset1]){
      for(j in 1:nclus_gene[dataset2]){
        mat[i,j] = length(intersect(genes_per[[dataname1]][[i]],genes_per[[dataname2]][[j]]))
        
      }
    }
    inter[[paste0(dataname1," x ",dataname2)]] = mat
    l=l+1
  }
}


names(clusters_best_1000[[1]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[1]]$sampInfo$sampID
data_val_filtered_1000[[1]]$sampInfo$clus = clusters_best_1000[[1]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[1]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[1]]$sampInfo$clus %in% c(1,3),2,
         ifelse(data_val_filtered_1000[[1]]$sampInfo$clus==2,1,
                ifelse(data_val_filtered_1000[[1]]$sampInfo$clus==4,2,NA)))


names(clusters_best_1000[[2]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[2]]$sampInfo$sampID
data_val_filtered_1000[[2]]$sampInfo$clus = clusters_best_1000[[2]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[2]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[2]]$sampInfo$clus %in% c(1,3),4,
         ifelse(data_val_filtered_1000[[2]]$sampInfo$clus==2,1,
                ifelse(data_val_filtered_1000[[2]]$sampInfo$clus==4,2,NA)))


names(clusters_best_1000[[3]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[3]]$sampInfo$sampID
data_val_filtered_1000[[3]]$sampInfo$clus = clusters_best_1000[[3]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[3]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[3]]$sampInfo$clus %in% c(1,2,4),4,
         ifelse(data_val_filtered_1000[[3]]$sampInfo$clus==3,1,NA))

names(clusters_best_1000[[4]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[4]]$sampInfo$sampID
data_val_filtered_1000[[4]]$sampInfo$clus = clusters_best_1000[[4]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[4]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[4]]$sampInfo$clus %in% c(1,2,4),3,
         ifelse(data_val_filtered_1000[[4]]$sampInfo$clus==3,1,NA))


names(clusters_best_1000[[5]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[5]]$sampInfo$sampID
data_val_filtered_1000[[5]]$sampInfo$clus = clusters_best_1000[[5]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[5]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[5]]$sampInfo$clus %in% c(1,2,4),3,
         ifelse(data_val_filtered_1000[[5]]$sampInfo$clus==3,1,NA))


names(clusters_best_1000[[6]]$clus_res$clusCol[[4]]$consensusClass)==data_val_filtered_1000[[6]]$sampInfo$sampID
data_val_filtered_1000[[6]]$sampInfo$clus = clusters_best_1000[[6]]$clus_res$clusCol[[4]]$consensusClass
data_val_filtered_1000[[6]]$sampInfo$assign = 
  ifelse(data_val_filtered_1000[[6]]$sampInfo$clus %in% c(1,2),2,
         ifelse(data_val_filtered_1000[[6]]$sampInfo$clus==3,1,
                ifelse(data_val_filtered_1000[[6]]$sampInfo$clus==4,4,NA)))


samps = list()
for(i in 1:6){
  samps[[i]] = data_val_filtered_1000[[i]]$sampInfo
}
library(survival)

sampInfo = dplyr::bind_rows(samps)
table(sampInfo$assign,sampInfo$DeCAF,sampInfo$PurIST)
fit=survfit(Surv(time,event)~assign,data=sampInfo)
ggsurvplot(fit,data=sampInfo)

clusters_best_1000[[dataset]]$data$dataname
ncol(clusters_best_1000[[dataset]]$data$ex)
clusters_best_1000[[dataset]]$data$sampInfo$subtype = 
  ifelse(clusters_best_1000[[dataset]]$data$sampInfo$Elyada_CAF=="myCAF","myCAF",
         ifelse(clusters_best_1000[[dataset]]$data$sampInfo$PurIST=="Classical","Classical","restCAF"))
nclus=4

table(clusters_best_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_best_1000[[dataset]]$data$sampInfo$PurIST)
table(clusters_best_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_best_1000[[dataset]]$data$sampInfo$DeCAF)
table(clusters_best_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_best_1000[[dataset]]$data$sampInfo$subtype)

table(clusters_best_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      paste0(clusters_best_1000[[dataset]]$data$sampInfo$DeCAF," * ",
             clusters_best_1000[[dataset]]$data$sampInfo$PurIST))



dataset=2
clusters_std_1000[[dataset]]$data$dataname
ncol(clusters_std_1000[[dataset]]$data$ex)
clusters_std_1000[[dataset]]$data$sampInfo$subtype = 
  ifelse(clusters_std_1000[[dataset]]$data$sampInfo$DeCAF=="permCAF","permCAF",
         ifelse(clusters_std_1000[[dataset]]$data$sampInfo$PurIST=="Classical","Classical","restCAF"))
nclus=2
table(clusters_std_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_std_1000[[dataset]]$data$sampInfo$DeCAF)
table(clusters_std_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_std_1000[[dataset]]$data$sampInfo$PurIST)
table(clusters_std_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      clusters_std_1000[[dataset]]$data$sampInfo$subtype)

table(clusters_std_1000[[dataset]]$clus_res$clusCol[[nclus]]$consensusClass,
      paste0(clusters_std_1000[[dataset]]$data$sampInfo$DeCAF," * ",
             clusters_std_1000[[dataset]]$data$sampInfo$PurIST))



roc_res_std <- timeROC(roc_res_std <- timeROC(SCISSORS_CAF
  T = lp_std$y,
  delta = lp_std$d,
  marker = lp_std$lp,
  cause = 1,
  times = seq(2,80,by=2),     # in months, for example
  iid = FALSE
)
roc_res_std$AUC

roc_res <- timeROC(
  T = lp_desurv$y,
  delta = lp_desurv$d,
  marker = lp_desurv$lp,
  cause = 1,
  times = seq(2,80,by=2),     # in months, for example
  iid = FALSE
)
roc_res$AUC
plot(roc_res, time = 36, col = "red")
plot(roc_res_std,time=36,col="blue")

lp_desurv$lp_bin = lp_desurv$lp>median(lp_desurv$lp)

fitsurv_de = survfit(Surv(y,d)~lp_bin,data=lp_desurv)
ggsurvplot(fitsurv_de,data=lp_desurv)

lp_std$lp_bin = lp_std$lp > median(lp_std$lp)
fitsurv_std = survfit(Surv(y,d)~lp_bin,data=lp_std)
ggsurvplot(fitsurv_std,data=lp_desurv)


### clustering
ntop=25
tops_std_1000 = get_top_genes(fit_std_beta$W,ntop)$top_genes
tops_best_1000 = get_top_genes(fit_consensus_1000$W,ntop)$top_genes

data_val = readRDS("data/derv/puleo_formatted.rds")
train_genes = rownames(data_filtered_1000$ex)
data_val_filtered = preprocess_data_val(data = data_val,
                                    genes = train_genes,
                                    method_trans_train = METHOD_TRANS_TRAIN)

clus_std = run_clustering(tops_std_1000,data_val_filtered,top_genes,colors,facs=c(2,3))

clus_de = run_clustering(tops_best_1000,data_val_filtered,top_genes,colors,facs=c(1,3))
nclus=4
data_val_filtered$sampInfo[[paste0("de_",nclus,"_",ntop)]]=as.factor(clus_de$clus_res$clusCol[[nclus]]$consensusClass)
data_val_filtered$sampInfo[[paste0("std_",nclus,"_",ntop)]]=as.factor(clus_std$clus_res$clusCol[[nclus]]$consensusClass)

fit1=coxph(Surv(time,event)~de_4_10,data=data_val_filtered$sampInfo)
fit2=coxph(Surv(time,event)~std_4_10,data=data_val_filtered$sampInfo)
summary(fit1)
summary(fit2)

cvwrapr::getCindex(clus_de$clus_res$clusCol[[nclus]]$consensusClass,
                   Surv(data_val_filtered$sampInfo$time,data_val_filtered$sampInfo$event))

cvwrapr::getCindex(clus_std$clus_res$clusCol[[nclus]]$consensusClass,
                   Surv(data_val_filtered$sampInfo$time,data_val_filtered$sampInfo$event))

set = "DeCAF"
nclus=2
table(clus_de$clus_res$clusCol[[nclus]]$consensusClass , 
      data_val_filtered$sampInfo[[set]])
table(clus_std$clus_res$clusCol[[nclus]]$consensusClass , 
      data_val_filtered$sampInfo[[set]])




data_val_filtered$sampInfo$subtype=paste0(data_val_filtered$sampInfo$PurIST," * ",
                                          data_val_filtered$sampInfo$DeCAF)
nclus=2
table(clus_de$clus_res$clusCol[[nclus]]$consensusClass , 
      data_val_filtered$sampInfo$subtype)
table(clus_std$clus_res$clusCol[[nclus]]$consensusClass , 
      data_val_filtered$sampInfo$subtype)


#scores
data_val = readRDS("data/derv/Dijk.Linehan.Moffitt_GEO_array.PACA_AU_array.PACA_AU_seq.Puleo_array_formatted.rds")
train_genes = rownames(data_filtered_1000$ex)#unlist(tops_best_1000[c(3)])
data_val_filtered = preprocess_data(data = data_val,
                                    genes = train_genes,
                                    method_trans_train = METHOD_TRANS_TRAIN)
ntop=20
tops_std = get_top_genes(fit_std_beta$W,ntop)$top_genes
tops_best = get_top_genes(fit_consensus_1000$W,ntop)$top_genes
scores_std = compute_scores(tops_std,fit_std_beta$W,data_val_filtered$ex,y,delta,score_bin = TRUE)
scores_de = compute_scores(tops_best,fit_consensus_1000$W,data_val_filtered$ex,y,delta,score_bin = TRUE)

length(intersect(tops_std[,3],tops_best[,1]))
length(intersect(tops_std[,2],tops_best[,3]))
fit1=coxph(Surv(time,event)~`1`+`2`+`3`,data=scores_std)
fit2=coxph(Surv(time,event)~`1`+`2`+`3`,data=scores_de)
summary(fit1)
summary(fit2)

scores_std$dataset = data_val_filtered$sampInfo$dataset
scores_de$dataset = data_val_filtered$sampInfo$dataset
datcur="Moffitt_GEO_array"
scores_std_1 = scores_std %>% filter(dataset==datcur)
scores_de_1 = scores_de %>% filter(dataset==datcur)
fit1=coxph(Surv(time,event)~`1`+`2`+`3`,data=scores_std_1)
fit2=coxph(Surv(time,event)~`1`+`2`+`3`,data=scores_de_1)
summary(fit1)
summary(fit2)

cvwrapr::getCindex(scores_std$`3`,Surv(scores_std$time,scores_std$event))
cvwrapr::getCindex(scores_de$`3`,Surv(scores_de$time,scores_de$event))
