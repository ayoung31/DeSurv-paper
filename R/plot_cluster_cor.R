plot_cluster_cor = function(clus,tops,nclus){
  dataset_names <- vapply(clus, function(x) x$data$dataname, character(1))
  if (is.null(names(nclus))) {
    if (length(nclus) != length(dataset_names)) {
      stop("Length of nclus must match number of datasets or include names.")
    }
    nclus <- setNames(nclus, dataset_names)
  }
  missing_sets <- setdiff(dataset_names, names(nclus))
  if (length(missing_sets)) {
    stop("nclus is missing entries for: ", paste(missing_sets, collapse = ", "))
  }

  ex =vector("list", length(clus))
  for(i in seq_along(clus)){
    current_clus <- clus[[i]]
    dataname <- dataset_names[[i]]
    nclus_i <- nclus[[dataname]]
    if (is.na(nclus_i) || nclus_i < 1) {
      stop("Invalid nclus value for dataset: ", dataname)
    }
    X = current_clus$data$ex
    genes = intersect(unlist(tops),rownames(X))
    X=as.data.frame(t(X[genes,]))
    X$dataset = dataname
    X$subject = rownames(X)
    
    X$cluster = as.character(current_clus$clus_res$clusCol[[nclus_i]]$consensusClass)
    
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
