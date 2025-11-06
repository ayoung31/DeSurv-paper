tar_load(clusters_desurv_50_2000)
clus=clusters_desurv_50_2000
tops=tops_desurv$top_genes
for(i in 1:length(clus)){
  dat=data.frame(matrix(ncol=3,nrow=0))
  colnames(dat) = c("cluster","immune","basal")
  nc=nclus_tbl$nclus[i]
  classes = clus[[i]]$clus_res$clusCol[[4]]$consensusClass
  nc=unique(classes)
  for(j in 1:4){
    
    subs = names(classes)[classes==j]
    immune=mean(as.numeric(clus[[i]]$data$ex[rownames(clus[[i]]$data$ex) %in% tops[,1],subs]))
    basal=mean(as.numeric(clus[[i]]$data$ex[rownames(clus[[i]]$data$ex) %in% tops[,3],subs]))
    dat[nrow(dat)+1,] = c(j,immune,basal)
  }
  clus[[i]]$data$dataname
  dat
}

subs = clus[[5]]$data$sampInfo$sampID[clus[[5]]$data$sampInfo$PurIST=="Basal-like"]
mean(as.numeric(clus[[5]]$data$ex[rownames(clus[[5]]$data$ex) %in% tops[,1],subs]))
mean(as.numeric(clus[[5]]$data$ex[rownames(clus[[5]]$data$ex) %in% tops[,3],subs]))

subs = clus[[5]]$data$sampInfo$sampID[clus[[5]]$data$sampInfo$PurIST=="Classical"]
mean(as.numeric(clus[[5]]$data$ex[rownames(clus[[5]]$data$ex) %in% tops[,1],subs]))
mean(as.numeric(clus[[5]]$data$ex[rownames(clus[[5]]$data$ex) %in% tops[,3],subs]))

     