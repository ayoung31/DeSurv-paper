build_factor_selection_table = 
  function(tops,path){
    genes = character()
    for(i in 1:ncol(tops)){
      genes[i] = paste(unlist(tops[,i]),collapse = ", ")
    }
    df_new = data.frame(
      factor = seq_len(ncol(tops)),
      selected = FALSE,
      genes = genes,
      stringsAsFactors = FALSE
    )
    
    if(file.exists(path)){
      df_old = read.csv(path)
      old_cols = setdiff(names(df_old), "selected")
      new_cols = setdiff(names(df_new), "selected")
      df_old_comp = df_old[, old_cols, drop = FALSE]
      df_new_comp = df_new[, new_cols, drop = FALSE]
      if(identical(df_old_comp,df_new_comp)){
        print("Existing factor selection table already corresponds to current runs. Not regenerating.")
        return(path)
      }
      # 1. make a timestamped backup of the old file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_path <- file.path(dirname(path),
                               paste0("factor_selection_", timestamp, ".csv")
      )
      file.copy(from = path, to = backup_path)
      
      message("Existing factor selection table backed up as: ", backup_path)
      
      # 2. overwrite cluster_review.csv with a clean table
      #    (i.e. you must re-curate for this new clustering)
      write.csv(df_new, path, row.names = FALSE)
      
      message("New factor selection table created from runs. Please re-curate.")
    } else {
      # first time: just create it
      write.csv(df_new, path, row.names = FALSE)
      message("factor selection table created. Please curate and mark TRUE for selected k.")
    }
    
    return(path)
}
