for(i in 1:length(clus)){
  dat=data.frame(matrix(ncol=3,nrow=0))
  colnames(dat) = c("cluster","immune","basal")
  nc=nclus_tbl$nclus[i]
  genes_immune=rownames(clus[[i]]$data$ex)[rownames(clus[[i]]$data$ex) %in% tops[,1]]
  genes_basal=rownames(clus[[i]]$data$ex)[rownames(clus[[i]]$data$ex) %in% tops[,3]]
  dat=data.frame(gene=rep(c(genes_immune,genes_basal),nc),
                 class=rep(c(rep("immune",length(genes_immune)),rep("basal",length(genes_basal))),nc),
                 cluster=rep(1:nc,each=length(c(genes_immune,genes_basal))))
  expr=numeric()
  for(j in 1:nc){
    classes = clus[[i]]$clus_res$clusCol[[nc]]$consensusClass
    subs = names(classes)[classes==j]
    immune=rowMeans(clus[[i]]$data$ex[genes_immune,subs])
    basal=rowMeans(clus[[i]]$data$ex[genes_basal,subs])
    
    expr=c(expr,immune,basal)
    
  }
  dat$expr=expr
  dat
  ggplot(dat,aes(y=expr,color=class))+
    geom_boxplot()+
    facet_wrap(~cluster)
}