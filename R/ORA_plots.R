
ORA_plots = function(results,list_name,msigdbr,p.adj){
  save_dir = file.path(results$model_save_dir,"dotplots")
  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive = TRUE)
  }
  ora_res_list = results[[list_name]]
  
  if(!msigdbr){
    plot_files <- paste0(save_dir, paste0("/dotplot_factor", 1:results$model.params$k, ".png"))
  }else{
    plot_files <- paste0(save_dir, paste0("/dotplot_factor", 1:results$model.params$k, "_msigdbr.png"))
  }
  
  if(all(file.exists(plot_files)[!is.na(results$labels$label)])){
    print("dotplots already exist")
    return()
  }
  
  
  if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
  }
  
  for (i in seq_along(ora_res_list)) {
    ora_res=ora_res_list[[i]]
    ora_df = as.data.frame(ora_res)
    
    if (!is.null(ora_res) && nrow(ora_df) > 0){
      p1 <- clusterProfiler::dotplot(ora_res_list[[i]],
                                     title = paste0("alpha=",results$alpha," factor ", i))
      p1 =p1+scale_fill_gradientn(
        name = "adjusted p-values",
        colours = c("red", "white", "blue"),
        limits = c(0, p.adj))
      png(plot_files[i], width = 480, height = 480)
      if (inherits(p1, "gg")) print(p1)
      dev.off()
    }
    
  }
  
  
}
