plot_cluster_cor = function(clus,tops,nclus){
  ex =list()
  for(i in 1:5){
    X = clus[[i]]$data$ex
    genes = intersect(unlist(tops),rownames(X))
    X=as.data.frame(t(X[genes,]))
    X$dataset = clus[[i]]$data$dataname
    X$subject = rownames(X)
    
    X$cluster = as.character(clus[[i]]$clus_res$clusCol[[nclus[i]]]$consensusClass)
    
    ex[[i]] = X
  }
  X = dplyr::bind_rows(ex)
  X = X %>% pivot_longer(cols = where(is.numeric),names_to="gene",values_to = "expr")
  X$expr[is.nan(X$expr)] = NA
  X = X %>% group_by(dataset,gene,cluster) %>% summarise(mean_expr = mean(expr,na.rm=TRUE))
  X=X%>% group_by(dataset,gene) %>% mutate(mean_expr = scale(mean_expr))
  X$mean_expr[is.nan(X$mean_expr)] = NA
  
  X = X %>% pivot_wider(names_from=c(cluster,dataset),values_from = mean_expr)
  
  palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  pheatmap(cor(X[,2:ncol(X)],use="pairwise.complete.obs"),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = palette,
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE,
           number_color = "black")
  
}
