create_table <- function(tops,gene_lists,which.lists,color.lists,text=NULL,
                            save=FALSE){
  library(dplyr)
  library(gt)
  
  if(which.lists=="all"){
    which.lists=names(gene_lists)
  }
  
  if(!(all(which.lists %in% names(gene_lists)))){
    stop("Some requested lists were not supplied in gene_lists")
  }
  
  idx = idx = which(names(gene_lists) %in% which.lists)
  gene_lists = gene_lists[which.lists]
  
  builder <- function(x, genes){cells_body(columns = !!rlang::sym(x), rows = !!rlang::sym(x) %in% genes)}
  
  tabs=list()
  for(i in 1:length(gene_lists)){
    filename=paste0("gene.overlap_",names(gene_lists)[i],".png")
    # if(!file.exists(paste0(save_dir,"/",filename))){
    colors = color.lists[[idx[i]]]
    text = c("black","white")[1+apply(col2rgb(colors),2,function(x) (sum(x * c(299,587,114))/1000) < 123)]
    
    tab1 <- tops %>% gt() %>% tab_header(title=names(gene_lists)[i],subtitle=paste(paste(colors,names(gene_lists[[i]]),sep="="),collapse = " "))
    for(j in 1:length(colors)){
      tab1 <- tab1 %>%
        tab_style(style=list(cell_fill(color=colors[j]),cell_text(color=text[j])),
                  locations=lapply(colnames(tops), builder, genes = gene_lists[[i]][[j]]))
      
    }
    tabs[[i]] = tab1

  }
  return(tabs)
}
