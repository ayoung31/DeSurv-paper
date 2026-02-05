build_std_nmf_k_selection_table = 
  function(k_vals, path, fit_std){
    df_new = fit_std$measures
    df_new$selected = FALSE
    
    if (file.exists(path)) {
      df_old = read.csv(path)
      df_old_comp = df_old %>% select(-selected)
      df_new_comp = df_new %>% select(-selected)
      if(identical(df_old_comp,df_new_comp)){
        print("Existing std_nmf_k_selection_table.csv already corresponds to current runs. Not regenerating.")
        return(path)
      }
      # 1. make a timestamped backup of the old file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_path <- file.path(dirname(path),
                               paste0("std_nmf_k_selection_table_", timestamp, ".csv")
                              )
      file.copy(from = path, to = backup_path)
      
      message("Existing std_nmf_k_selection_table.csv backed up as: ", backup_path)
      
      # 2. overwrite cluster_review.csv with a clean table
      #    (i.e. you must re-curate for this new clustering)
      write.csv(df_new, path, row.names = FALSE)
      
      message("New std_nmf_k_selection_table.csv created from runs. Please re-curate.")
    } else {
      # first time: just create it
      write.csv(df_new, path, row.names = FALSE)
      message("std_nmf_k_selection_table.csv created. Please curate and mark TRUE for selected k.")
    }
    return(path)
}