
determine_clustering_type <- function(tops, gene_lists=NULL, type, facs=NULL) {
  if(!is.null(facs)){
    type=facs
  }
  if (type == "bas/clas") {
    bas_clas <- id_bas_clas(tops, gene_lists)
    facs <- na.omit(bas_clas$factor)
    
    if (length(facs) == 0) stop("No basal or classical factor detected")
    if (length(facs) == 1) {
      kept <- bas_clas$type[facs]
      if (kept == "basal") {
        warning("No clear classical factor detected, clustering on Basal only")
        name <- "basal"
      } else if (kept == "classical") {
        warning("No clear basal factor detected, clustering on Classical only")
        name <- "classical"
      } else {
        stop("Unrecognized factor type")
      }
    } else {
      name <- "basal.and.classical"
    }
  } else if(type=="tumor"){
    # find tumor factors
    labels = results$labels$label
    cf = results$labels[which(grepl('classical',labels, ignore.case = TRUE)),]
    bf = results$labels[which(grepl('basal', labels, ignore.case = TRUE)),]
    
    facs = c(cf$factor[which.min(cf$p.adjust)],bf$factor[which.min(bf$p.adjust)])
    facs=facs[order(facs)]
    if(length(facs)==0){
      warning('no tumor factors detected')
    }
    name = paste0("k=",results$model.params$k,"_top",results$ntop,results$top.type,"_tumor_",paste0(facs,collapse="_"),"_alpha=",results$alpha)
    title = paste0("k=",results$model.params$k," top ",results$ntop,results$top.type," alpha=",results$alpha, " tumor factors ",paste0(facs,collapse = ", "))
  } else if(type=="stroma"){
    # find stromal factors
    labels = results$labels$label
    facs = which(Reduce(`|`, lapply(c('caf','stroma','PurISS','MS'), grepl, x = labels, ignore.case = TRUE)))
    facs=facs[order(facs)]
    if(length(facs)==0){
      warning("no stromal factors detected")
    }
    name = paste0("k=",results$model.params$k,"_stroma")
    name = paste0("k=",results$model.params$k,"_top",results$ntop,results$top.type,"_stroma_",paste0(facs,collapse = "_"),"_alpha=",results$alpha)
    title = paste0("k=",results$model.params$k," top ",results$ntop,results$top.type," alpha=",results$alpha, " stroma factors ",paste0(facs,collapse = ", "))
  } else if(is.numeric(type)){
    facs = type
    name = paste0("k=",results$model.params$k,"_top",results$ntop,results$top.type,"_factor",type,"_alpha=",results$alpha)
    title = paste0("k=",results$model.params$k," top ",results$ntop,results$top.type," alpha=",
                   results$alpha, " factor ",type,", ",results$labels$label[type])
    
  } else if(!is.null(facs)){
    facs = facs
    facs=facs[order(facs)]
    # name = paste0("k=",results$model.params$k,"_top",results$ntop,results$top.type,"_",type,"_",paste0(facs,collapse="_"),"_alpha=",results$alpha)
    # title = paste0("k=",results$model.params$k," top ",results$ntop,results$top.type," alpha=",results$alpha, " ", type," ",paste0(facs,collapse = ", "))
  }else{
    stop("Unsupported clustering type")
  }
  
  
  list(facs = facs)
}
