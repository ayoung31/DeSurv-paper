build_ncluster_selection_table = 
  function(clus,path){
    datanames = sapply(clus,function(x) x$data$dataname)
    maxclus = unique(sapply(clus,function(x) length(x$clus_res$clusCol)))
    if(length(maxclus)!=1){
      stop("All validation datasets should have the same number of max clusters")
    }
    tbl = data.frame(dataset=datanames,nclus=NA)
    colnames(tbl) = c("dataset","nclus")
    temp=as.data.frame(t(sapply(clus,function(x){
      sapply(x$clus_res$clusCol[2:maxclus],function(y) paste(y$consensusClass,collapse = ","))
    })))
    colnames(temp) = paste0("assignments_nclus",2:maxclus)
    temp$dataset=datanames
    tbl_new = tbl %>% left_join(temp)
    
    if(file.exists(path)){
      tbl_old = read.csv(path)
      tbl_old_comp = tbl_old[,paste0("assignments_nclus",2:maxclus)]
      tbl_new_comp = tbl_new[,paste0("assignments_nclus",2:maxclus)]
      
      if(identical(tbl_old_comp,tbl_new_comp)){
        print("Existing ncluster selection table already corresponds to current runs. Not regenerating.")
        return(path)
      }
      
      # 1. make a timestamped backup of the old file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_path <- file.path(dirname(path),
                               paste0("ncluster_selection_", timestamp, ".csv")
      )
      file.copy(from = path, to = backup_path)
      
      message("Existing ncluster selection table backed up as: ", backup_path)
      
      write.csv(tbl_new, path, row.names = FALSE)
      
      message("New ncluster selection table created from runs. Please re-curate.")
    }else{
      write.csv(tbl_new, path, row.names = FALSE)
      message("ncluster selection table created. Please curate and mark TRUE for selected ncluster for each dataset.")
    }
    
    return(path)
  }