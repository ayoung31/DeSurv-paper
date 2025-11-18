# tar_load(data_val)
ntop=50
method="std"
tar_load_globals()
# tar_load(paste0("tops_",method,"_",ntop,"_1000"))
# # tar_load(paste0("clusters_",method,"_",ntop,"_1000"))
# # clus=get(paste0("clusters_",method,"_",ntop,"_1000"))
# tops = get(paste0("tops_",method,"_",ntop,"_1000"))
# tops = tops$top_genes
tar_load(data_val)
clus=list()
for(i in 1:5){
  dataset <- data_val[[i]]
  dataset <- DeSurv::preprocess_data(
    X = dataset$ex,
    y = dataset$sampInfo$time,
    d = dataset$sampInfo$event,
    dataset = dataset$sampInfo$dataset,
    samp_keeps = dataset$samp_keeps,
    ngene = 1000,
    method_trans_train = "quant",
    verbose = FALSE
  )
  dataname <- if (!is.null(names(data_val))) names(data_val)[i] else NULL
  if (is.null(dataname) || !nzchar(dataname)) {
    dataname <- dataset$dataname
  }
  if (is.null(dataname) || !nzchar(dataname)) {
    dataname <- paste0("dataset_", i)
  }
  dataset$dataname <- dataname
  data_val[[i]] = dataset
  clus[[i]]=run_clustering_v2(tops,data_val[[i]],top_genes,facs=c(2,6),
                              plot=FALSE,maxKcol = 5,maxKrow = 5)
}
clus_save=clus
# which factors to use?
create_table(tops,top_genes,"DECODER",colors)

tops = unlist(tops[,c(2,6)])

nclus=rep(3,5)#c(3,3,2,2,3)
plot_cluster_cor(clus,tops,nclus)

gene_cluster_overlap(clus,2)

#manually align clusters
#top50
# dijk =   c(1,2)
# moff =   c(1,2)
# paca_a = c(2,1)
# paca_s = c(1,2)
# pul =    c(1,2)

#top25 
dijk =   c(1,2)
moff =   c(1,2)
paca_a = c(2,1)
paca_s = c(2,1)
pul =    c(1,2)




gene_clus=data.frame(dijk=dijk,moff=moff,paca_a=paca_a,paca_s=paca_s,pul=pul)


for(i in 1:5){
  nms = names(clus[[i]]$clus_res$clusRow[[2]]$consensusClass)
  y = match(clus[[i]]$clus_res$clusRow[[2]]$consensusClass,gene_clus[,i])
  names(y) = nms
  clus[[i]]$clus_res$clusRow[[2]]$consensusClass=y
}


#subject cluster assignments
#top50
# dijk =   c(1,3,2)
# moff =   c(1,3,2)
# paca_a = c(1,2,3)
# paca_s = c(2,1,3)
# pul =    c(3,2,1)

#top25 2clus moff pul
# dijk =   c(2,3,1)
# moff =   c(2,NA,1)
# paca_a = c(2,3,1)
# paca_s = c(3,1,2)
# pul =    c(2,NA,1)

#top 25 3 clus all
# dijk =   c(1,2,3)
# moff =   c(1,3,2)
# paca_a = c(3,1,2)
# paca_s = c(2,1,3)
# pul =    c(1,3,2)
# dijk =   c(1,2,3,NA)
# moff =   c(1,2,3,NA)
# paca_a = c(2,1,1,3)
# paca_s = c(2,1,1,3)
# pul =    c(3,2,1,NA)
dijk =   c(1,3,2)
moff =   c(1,3,2)
paca_a = c(2,3,1)
paca_s = c(1,2,3)
pul =    c(3,1,2)

samp_clus = data.frame(dijk=dijk,moff=moff,paca_a=paca_a,paca_s=paca_s,pul=pul)

for(i in 1:5){
  
  nms = names(clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass)
  # y = samp_clus[clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass,i]
  y = match(clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass,samp_clus[,i])
  names(y) = nms
  clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass=y
  
}

# clus[[2]]$clus_res$clusCol[[3]]$consensusClass[clus[[2]]$clus_res$clusCol[[3]]$consensusClass==3]

#double check with heatmap
for(i in 1:5){
  print(plot_heatmap(clus[[i]]$clus_res,k=nclus[i],k_row=2,Xtemp=clus[[i]]$data$ex))
}

for(i in 1:5){
  print(table(clus[[i]]$data$sampInfo$sampID==names(clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass)))
  clus[[i]]$data$sampInfo$samp_cluster = clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass
}

#save clustering alignment
# if(method=="best")method="DeSurv"
# dir.create(paste0("results/PKG_VERSION=lambda_scaling_GIT_BRANCH=fold_filtering/TCGA_PAAD.CPTAC/rank/ng1000_tol1e-04_max6000/clustering/aligned/ntop=",ntop,"/",method,"/"),recursive = TRUE)
# save(clus,file=paste0("results/PKG_VERSION=lambda_scaling_GIT_BRANCH=fold_filtering/TCGA_PAAD.CPTAC/rank/ng1000_tol1e-04_max6000/clustering/aligned/ntop=",ntop,"/",method,"/aligned_clusters.RData"))



i=5
fit = survfit(Surv(time,event)~samp_cluster,data=clus[[i]]$data$sampInfo)
ggsurvplot(fit,data=clus[[i]]$data$sampInfo,title=clus[[i]]$data$dataname)

sinfos=list()
for(i in 1:5){
  sinfos[[i]]=clus[[i]]$data$sampInfo
}
sampInfo=bind_rows(sinfos)

fit = survfit(Surv(time,event)~samp_cluster,data=sampInfo)
p=ggsurvplot(fit,data=sampInfo,pval = TRUE,risk.table = TRUE,
             legend.labs=c("1","2","3"),legend.title="Cluster",
             pval.coord=c(0,.05))
coxph(Surv(time,event)~as.factor(samp_cluster)+strata(dataset),data=sampInfo)

bottom=plot_grid(NULL,p$table,ncol=2,rel_widths = c(.2,8))
survplot=plot_grid(p$plot,bottom,nrow=2,rel_heights = c(2,1))
ggsave(filename="paper/figures/survival_clusters_stdNMF_bas_immune_3clus.tiff",
       plot=survplot,width=5,height=6)

sampInfo$subtype = ifelse(sampInfo$PurIST=="Basal-like","Basal",
                          ifelse(sampInfo$DeCAF=="restCAF","restCAF","proCAF"))
table(sampInfo$samp_cluster,sampInfo$subtype)
table(sampInfo$samp_cluster,sampInfo$PurIST)
table(sampInfo$samp_cluster,sampInfo$DeCAF)

tar_load(data_val)
for(i in 1:5){
  data_val[[i]]$ex = data_val[[i]]$ex[,data_val[[i]]$samp_keeps]
  data_val[[i]]$sampInfo = data_val[[i]]$sampInfo[data_val[[i]]$samp_keeps,]
}

ex=list()
genes = list()
samps = list()
gene_clus = list()
tops = unlist(tops_all[,c(2,6)])
for(i in 1:5){
  ex[[i]] = data_val[[i]]$ex
  genes[[i]] = intersect(rownames(ex[[i]]),tops)
  samps[[i]] = clus[[i]]$data$sampInfo
  gene_clus[[i]] = clus[[i]]$clus_res$clusRow[[2]]$consensusClass
}

# genes = Reduce(intersect,lapply(gene_clus,function(x) names(x)))
genes=Reduce(intersect,genes)
ex=lapply(ex,function(x) t(scale(t(x[genes,]))))
X = do.call('cbind',ex)
X = X[genes,]
sampInfo = do.call('rbind',samps)
col_anno = data.frame(NMF_cluster=as.factor(sampInfo$samp_cluster),
                      DeCAF = sampInfo$DeCAF,
                      PurIST = sampInfo$PurIST,
                      dataset=sampInfo$dataset
)
col_anno$DeCAF = ifelse(col_anno$DeCAF=="permCAF","proCAF","restCAF")
rownames(col_anno) = colnames(X)
X = X[,order(col_anno$NMF_cluster)]
col_anno = col_anno[order(col_anno$NMF_cluster),]
row_anno = data.frame(gene = genes)#,gene_cluster=as.factor(gene_clus[[1]][genes])
row_anno$`NMF factor` = as.factor(ifelse(genes %in% tops_all[,2],2,
                                         ifelse(genes %in% tops_all[,6],6,7)))
rownames(row_anno) = row_anno$gene
row_anno$gene = NULL
X = X[order(row_anno$`NMF factor`),]
row_anno = row_anno[order(row_anno$`NMF factor`),,drop=FALSE]

min_val = -2
max_val = 4
ncolors <- 200
my_colors <- colorRampPalette(c("blue", "white", "red"))(ncolors)

# Create symmetric breaks centered at 0
breaks_centered <- c(
  seq(min_val, 0, length.out = ceiling(ncolors / 2) + 1),
  seq(0, max_val, length.out = floor(ncolors / 2) + 1)[-1]
)
annotation_colors = list(
  `NMF factor`=c(
    # `3` = "black",
    `2` = "black",
    `6` = "lightgrey"
  ),
  NMF_cluster = c(
    `1` = "#F8766D",
    `2` = "#00BA38",
    `3` = "#619CFF"
  ),
  DeCAF = c(
    `proCAF` = "violetred2",
    `restCAF` = "cyan4"
  ),
  PurIST = c(
    `Basal-like` = "orange",
    `Classical` = "blue"
  ),
  dataset = c(
    `Dijk` = "slateblue1",
    `Moffitt_GEO_array` = "springgreen4",
    `PACA_AU_array` = "yellow3",
    `PACA_AU_seq` = "coral",
    `Puleo_array` = "dodgerblue3"
  )
)
png("paper/figures/heatmap_stdNMF_bas_immune_3clus.png",width=8,height=6,units="in",res=600)
pheatmap(X,
         annotation_col = col_anno,
         annotation_row = row_anno,
         annotation_colors = annotation_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = my_colors,
         breaks = breaks_centered,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         show_rownames = FALSE)
dev.off()
