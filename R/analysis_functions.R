## key functions

#source('code/utils/Plot_heatmap_CC.R')
#source('code/utils/heatmap.3.R')



survival_summary <- function(data,k,ntop,type,top.type="") {
  library(tidyr)
  library(dplyr)
  library(stringr)
  
  if (!"sampInfo" %in% names(data)) stop("data must contain 'sampInfo'")
  
  sampInfo <- data$sampInfo %>%
    dplyr::filter(whitelist, pdac)
  
  # grab all columns of sampInfo with ntop and type
  cols = grep(paste0("k=",k,"_top",ntop,top.type,"_",type,"\\d*_alpha=\\w+_\\d+clusters"), 
              colnames(sampInfo), value = TRUE)
  
  if(length(cols)==0){
    print(paste0("clustering results not available for k=",k))
    return(NULL)
  }
  
  res = data.frame(matrix(ncol=7,nrow=0))
  colnames(res) = c("type","ntop","k","factor","alpha","nclus","bic")
  
  for(c in cols){
    
    # get nclus, k, alpha, and factor number for current col
    pattern <- paste0("k=(\\d+)_top",ntop,top.type,"_",type,"(\\d*)_alpha=(\\w+)_(\\d+)clusters")
    matches <- regexec(pattern, c)
    params = regmatches(c, matches)[[1]][-1]
    
    if (length(params) != 4) next
    
    k = as.numeric(params[1])
    fac = as.numeric(params[2])
    alpha = params[3]
    nclus = as.numeric(params[4])
    
    sampInfo$cluster = sampInfo[[c]]
    
    if(length(table(sampInfo$cluster))>1){
      
      # Fit survival curves
      fit <- coxph(Surv(time, event) ~ as.factor(cluster), data = sampInfo)
      bic <- BIC(fit)
      # hr <- summary(fit2)$coefficients[, 2]
      
    }else{
      bic = as.numeric(NA)
    }
    
    res[nrow(res)+1,] = c(type,ntop,k,fac,alpha,nclus,bic)
    
    sampInfo$cluster = NULL
  }
  
  return(res)
}









full_pipeline = function(data_orig,results_a0,results_best,top_genes,colors,pairs,
                         replace=FALSE,weight=FALSE,pItem=.8,reps=1000){
  
  if(results_a0$ntop != results_best$ntop){
    stop("Number of top genes for the two models must match")
  }
  ntop=basename(results_a0$model_save_dir)
  base_dir = dirname(dirname(dirname(results_a0$model_save_dir)))
  mod1 = basename(dirname(results_a0$model_save_dir))
  mod2 = basename(dirname(results_best$model_save_dir))
  doc_dir = file.path(base_dir,"docs",paste0(mod1,"_VS_",mod2),ntop)
  if(!dir.exists(doc_dir)){
    dir.create(doc_dir,recursive = TRUE)
  }
  doc_name = paste0(data_orig$dataname,".RData")
  
  data_dir1 = file.path(base_dir,"data",mod1,ntop)
  if(!dir.exists(data_dir1)){
    dir.create(data_dir1,recursive = TRUE)
  }
  data_dir2 = file.path(base_dir,"data",mod2,ntop)
  if(!dir.exists(data_dir2)){
    dir.create(data_dir2,recursive = TRUE)
  }
  data_name = paste0(data_orig$dataname,".RData")
  

  
  #if(!file.exists(paste0(doc_dir,doc_name)) | replace){
  # suppressMessages(
  #   suppressWarnings({
  
  ###### Model at alpha=0 #######
  if(file.exists(file.path(data_dir1,data_name))){
    load(file.path(data_dir1,data_name))
  }else{
    data1=data_orig
  }
  
  ## Cluster on these factors
  print("alpha=0, each")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_factor",1:model.params$k,"_alpha=0")
  # if(any(sapply(clus_res,function(x) is.null(data[[x]]))) | replace){
  data1 = run_clustering(results_a0,data1,top_genes,colors,type="each",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  # save(data1,file=file.path(data_dir1,data_name))
  # }
  
  ## Cluster on tumor
  print("alpha=0 tumor")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_tumor","_alpha=0")
  #if(is.null(data[[clus_res]]) | replace){
  data1 = run_clustering(results_a0,data1,top_genes,colors,type="tumor",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  # save(data1,file=file.path(data_dir1,data_name))
  #}
  
  ## Cluster on stroma
  print("alpha=0 stroma")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_stroma","_alpha=0")
  #if(is.null(data[[clus_res]]) | replace){
  data1 = run_clustering(results_a0,data1,top_genes,colors,type="stroma",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  save(data1,file=file.path(data_dir1,data_name))
  #}
  
  ## Cluster on nonzero betas
  # print("alpha=0 survival")
  # # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_survival","_alpha=0")
  # # if(is.null(data[[clus_res]]) | replace){
  #   data = run_clustering(results_a0,data,top_genes,colors,type="surv",save=TRUE,
  #                         pItem=pItem, reps=reps, weight=weight, replace=replace)
  #   save(data,file=data_file)
  # # }
  
  ## Compare clustering results to known labels
  comp_class_known(data1,results_a0)
  
  comp_class_known(data1,results_a0,"tumor")
  
  comp_class_known(data1,results_a0,"stroma")
  
  # comp_class_known(data,results_a0,"surv")
  
  
  
  ####### model at best alpha ##########
  if(file.exists(file.path(data_dir2,data_name))){
    load(file.path(data_dir2,data_name))
  }else{
    data2=data_orig
  }
  
  
  
  ## Cluster on these factors
  print("alpha=best, each")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_factor",1:model.params$k,"_alpha=best")
  # if(any(sapply(clus_res,function(x) is.null(data[[x]]))) | replace){
  data2 = run_clustering(results_best,data2,top_genes,colors,type="each",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  # save(data2,file=file.path(data_dir2,data_name))
  # }
  
  ## Cluster on tumor
  print("alpha=best, tumor")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_tumor","_alpha=best")
  #if(is.null(data[[clus_res]]) | replace){
  data2 = run_clustering(results_best,data2,top_genes,colors,type="tumor",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  # save(data2,file=file.path(data_dir2,data_name))
  #}
  
  ## Cluster on stroma
  print("alpha=best, stroma")
  # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_stroma","_alpha=best")
  #if(is.null(data[[clus_res]]) | replace){
  data2 = run_clustering(results_best,data2,top_genes,colors,type="stroma",save=TRUE,
                        pItem=pItem, reps=reps, weight=weight, replace=replace)
  save(data2,file=file.path(data_dir2,data_name))
  #}
  
  ## Cluster on nonzero betas
  # print("alpha=best survival")
  # # clus_res = paste0("k=",model.params$k,"_top",results_a0$ntop,"_survival","_alpha=best")
  # # if(is.null(data[[clus_res]]) | replace){
  #   data = run_clustering(results_best,data,top_genes,colors,type="surv",save=TRUE,
  #                         pItem=pItem, reps=reps, weight=weight, replace=replace)
  #   save(data,file=data_file)
  # # }
  
  ## Compare clustering results to known labels
  
  comp_class_known(data2,results_best)
  
  comp_class_known(data2,results_best,"tumor")
  
  comp_class_known(data2,results_best,"stroma")
  
  # comp_class_known(data,results_best,"surv")
  
  
  
  ######## comparison between alphas #########
  # comp_class(data=data,pairs=pairs,results_a0=results_a0,results_best=results_best)
  # comp_class(data=data,results_a0=results_a0,results_best=results_best,type="tumor")
  # comp_class(data=data,results_a0=results_a0,results_best=results_best,type="stroma")
  # comp_class(data=data,results_a0=results_a0,results_best=results_best,type="surv")
  
  #   })
  # )
  
  ###### generate R markdown ######
  # model.params = results_a0$model.params
  # ntop = results_a0$ntop
  # 
  # 
  # 
  # tryCatch(
  #   expr = rmarkdown::render('code/generate figures/full_results.Rmd',
  #                            output_file = doc_name,
  #                            output_dir = doc_dir,
  #                            params = list(doc_title = paste0(data$dataname,' results k=',model.params$k),
  #                                          dir=file.path(base_dir,"figures"),
  #                                          wd='/work/users/a/y/ayoung31/DeSurv-paper/',
  #                                          k=model.params$k,
  #                                          pairs=pairs,
  #                                          dataname=data$dataname,
  #                                          results1=results_a0,
  #                                          results2=results_best,
  #                                          mod1=mod1,
  #                                          mod2=mod2,
  #                                          ntop=ntop)),
  #   error = function(cond) {
  #     message(cond)
  #   },
  #   warning = function(cond) {
  #     message(cond)
  #   }
  # )
  #}
  return()
}

