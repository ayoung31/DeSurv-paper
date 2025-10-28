gene_cluster_overlap = function(clus,nclus){

  genes_per=list()
  for(dataset in 1:5){
    dataname=clus[[dataset]]$data$dataname
    X=clus[[dataset]]$data$ex
    
    row_class = clus[[dataset]]$clus_res$clusRow[[nclus]]$consensusClass
    row_lab = data.frame(gene_clus = row_class,gene=names(row_class))
    
    for(g in 1:nclus){
      genes_per[[dataname]][[g]] = row_lab$gene[row_lab$gene_clus==g]
    }
  }
  
  inter=list()
  
  for(dataset1 in 1:4){
    for(dataset2 in (dataset1 + 1):5){
      dataname1 = clus[[dataset1]]$data$dataname
      dataname2 = clus[[dataset2]]$data$dataname
      mat = matrix(nrow=nclus,ncol=nclus)
      for(i in 1:nclus){
        for(j in 1:nclus){
          mat[i,j] = length(intersect(genes_per[[dataname1]][[i]],genes_per[[dataname2]][[j]]))
          
        }
      }
      inter[[paste0(dataname1," x ",dataname2)]] = mat
    }
  }
  return(inter)
}
