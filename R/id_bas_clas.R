
### identify basal and classical factors
id_bas_clas = function(tops, gene_lists){
  library(clusterProfiler)
  
  # format gene lists
  # gene_lists = unlist(gene_lists,recursive = FALSE)
  
  # run ORA
  labels = ORA(tops, gene_lists, FALSE)
  
  # there may be more than one basal or classical factor
  basal = which(grepl('basal',labels$label,ignore.case = TRUE))
  classical = which(grepl('classical',labels$label,ignore.case = TRUE))
  
  # pick the best one
  if(length(basal>1)){
    basal = basal[which.min(as.numeric(labels$p.adjust[basal]))]
  }
  if(length(classical)>1){
    classical = classical[which.min(as.numeric(labels$p.adjust[classical]))]
  }
  if(length(basal)==0){
    basal = NA
  }
  if(length(classical)==0){
    classical = NA
  }
  bas_clas = data.frame(type=c("basal","classical"),factor=c(basal,classical))
  
  return(bas_clas)
  
}