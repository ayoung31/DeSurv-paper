
ORA_plots = function(ora_res_list,msigdbr,p.adj){
  
  for (i in seq_along(ora_res_list)) {
    ora_res=ora_res_list[[i]]
    ora_df = as.data.frame(ora_res)
    
    if (!is.null(ora_res) && nrow(ora_df) > 0){
      p1 <- clusterProfiler::dotplot(ora_res_list[[i]],
                                     title = paste0("factor ", i))
      p1 =p1+scale_fill_gradientn(
        name = "adjusted p-values",
        colours = c("red", "white", "blue"),
        limits = c(0, p.adj))
      print(p1)

    }
    
  }
  
  
}
