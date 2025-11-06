build_cluster_alignment_table = 
  function(nclus_tbl,clus,path){
    tbl = data.frame(matrix(ncol=nrow(nclus_tbl)+1,nrow=max(nclus_tbl$nclus)))
    colnames(tbl) = c("dataset_cluster",nclus_tbl$dataset)
    tbl$dataset_cluster = 1:nrow(tbl)
    temp=as.data.frame(sapply(clus,function(x){
      dataname=x$data$dataname
      maxclus = max(nclus_tbl$nclus)
      nclus = nclus_tbl$nclus[nclus_tbl$dataset==dataname]
      classes = x$clus_res$clusCol[[nclus]]$consensusClass
      subs=character()
      for(i in 1:maxclus){
        if(i <= nclus){
          subs[i] = paste0(names(classes)[classes==i],collapse=",")
        }else{
          subs[i] = NA
        }
        
      }
      subs
    }))
    colnames(temp) = paste0(nclus_tbl$dataset,"_subjects")
    temp$dataset_cluster = 1:max(nclus_tbl$nclus)
    tbl_new = tbl %>% left_join(temp)
    if(file.exists(path)){
      tbl_old = read.csv(path)
      tbl_old_comp = tbl_old[,paste0(nclus_tbl$dataset,"_subjects")]
      tbl_new_comp = tbl_new[,paste0(nclus_tbl$dataset,"_subjects")]
      
      if(identical(tbl_old_comp,tbl_new_comp)){
        print("Existing cluster alignment table already corresponds to current runs. Not regenerating.")
        return(path)
      }
      
      # 1. make a timestamped backup of the old file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_path <- file.path(dirname(path),
                               paste0("cluster_alignment_", timestamp, ".csv")
      )
      file.copy(from = path, to = backup_path)
      
      message("Existing ncluster selection table backed up as: ", backup_path)
      
      write.csv(tbl_new, path, row.names = FALSE)
      
      message("New cluster alignment table created from runs. Please re-curate.")
    }else{
      write.csv(tbl_new, path, row.names = FALSE)
      message("Cluster alignment table created. Please curate.")
    }
    
    return(path)
  }