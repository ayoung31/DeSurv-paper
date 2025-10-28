# _targets.R
library(targets)
library(tarchetypes)  
library(crew.cluster)
library(crew)
suppressWarnings(suppressMessages(library(dplyr)))

NINIT=20
NINIT_COLD = 5
NINIT_FULL = 100

# ------ Slurm controllers ------
default_controller = crew_controller_sequential()

low_mem_controller = crew_controller_slurm(
  name = "low_mem",
  workers = 202,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 1,
    time_minutes = 120,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

cv_comp_controller = crew_controller_slurm(
  name = "cv",
  workers = 202,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_per_cpu = 2,
    cpus_per_task = NINIT,
    time_minutes = 600,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

full_run_controller = crew_controller_slurm(
  name = "full",
  workers = 202,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 32,
    cpus_per_task = NINIT_FULL,
    time_minutes = 600,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

med_mem_controller = crew_controller_slurm(
  name = "med_mem",
  workers = 120,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 8,
    cpus_per_task = 1,
    time_minutes = 200,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

cv_cold_controller = crew_controller_slurm(
  name = "cv_cold",
  workers = 120,
  seconds_idle = 120,
  seconds_interval = 0.25,
  options_cluster = crew_options_slurm(
    memory_gigabytes_required = 12,
    cpus_per_task = NINIT_COLD,
    time_minutes = 600,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out",
    script_lines = "module load r/4.4.0"
  )
)

# ---- Global options ----
tar_option_set(
  packages = c("coxNMF","NMF","tidyverse","survival","cvwrapr","rmarkdown","dplyr",
               "parallel","foreach", "doParallel", "doMC", "pec", "glmnet"),
  format = "rds",
  controller = crew_controller_group(default_controller, 
                                     low_mem_controller,
                                     cv_comp_controller,
                                     cv_cold_controller,
                                     full_run_controller,
                                     med_mem_controller),
  error = "continue"
)

# ---- Source helper functions ----
purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

# ---- Reusable constants ----

VAL_DATASETS       = c("Dijk","Moffitt_GEO_array",
                       "PACA_AU_array","PACA_AU_seq","Puleo_array")
METHOD_SELECT_INIT = "surv" ### method for selecting best initialization surv=PL, nmf=recon
ALPHA         = seq(0, 1, by = .05)

# ---- Training parameters ----
TRAIN_DATASETS     = c("TCGA_PAAD","CPTAC")  
TRAIN_PREFIX       = paste0(TRAIN_DATASETS, collapse = ".")
METHOD_TRANS_TRAIN = "rank"
NGENE              = 2000#c(1000,2500)#c(1000,5000)
TOL                = 1e-4
MAXIT              = 6000
K_VALS             = 3:4#2:12#2:16   #= c(2,3,4,5)
LAMBDA_VALS        = c(.1,1)#c(1e-4,1e-3,1e-2,.1,1,10)#,seq(.2,.9,by=.1)   #10^seq(-3,3)#10^seq(-4,4)
ETA_VALS           = .7#seq(0,1,by=.1)#c(0,.01,.1,.5,.9)#c(0,.01)#c(.01,.1,.5,.9)#seq(.1,.9,by=.1)
LAMBDAW_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
LAMBDAH_VALS       = 0#10^seq(-3,3)#10^seq(-4,4)
NTOP               = 50#c(25,50,75)
NFOLD              = 5
PKG_VERSION        = utils::packageDescription("coxNMF", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()

# load top_genes and colors into global environment 
# (known gene lists and corresponding color schemes)
load("data/derv/cmbSubtypes_formatted.RData")

# ---- Targets ----
list(
  
  ################# Training ###################
  
  # set up parameter grid
  tar_target(
    param_list,
    {
      p = create_param_grid_cv(k=K_VALS,lambda=LAMBDA_VALS,eta=ETA_VALS,
                           lambdaW=LAMBDAW_VALS,lambdaH=LAMBDAH_VALS,
                           nfold=NFOLD)
    }
             
  ),
  
  # Raw data files for training datasets
  tar_target(
    raw_data,
    c(
      file.path("data/original", paste0(TRAIN_DATASETS, ".rds")),
      file.path("data/original", paste0(TRAIN_DATASETS, ".survival_data.rds")),
      file.path("data/original", paste0(TRAIN_DATASETS, "_subtype.csv"))
    ),
    format = "file"
  ),
  
  # Load and format raw data
  tar_target(
    data,
    { 
      raw_data
      load_data(TRAIN_DATASETS) 
    }
  ),
  
  # validation datasets
  tar_target(
    val_datasets,
    VAL_DATASETS,
    iteration="vector"
  ),
  
  tar_target(
    raw_data_val,
    c(
      file.path("data/original", paste0(val_datasets, ".rds")),
      file.path("data/original", paste0(val_datasets, ".survival_data.rds")),
      file.path("data/original", paste0(val_datasets, "_subtype.csv"))
    ),
    format = "file",
    pattern = map(val_datasets)
  ),
    
  #note that I map over val datasets here instead of loading them into one merged set
  tar_target(
    data_val,
    setNames(
      lapply(val_datasets, load_data),
      val_datasets
    )
  ),
  
  tar_target(
    data_val_comb,
    {
      raw_data_val
      load_data(val_datasets)
    },
  ),

  
  tar_map(
    values = tibble(ngene=NGENE),
    names = ngene,
    
    # Preprocess and filter data
    tar_target(data_filtered, preprocess_data(data = data,
                                              ngene = ngene,
                                              method_trans_train = METHOD_TRANS_TRAIN)),
    
    
    # Split data into folds, preprocess and filter each fold
    tar_target(
      data_folds,
      {
        set_folds(data = data, nfold = NFOLD, ngene=ngene, seed=123,
                  method_trans_train = METHOD_TRANS_TRAIN)
      }
    ),
    
    
    # Perform cross validation across parameter grid
    tar_target(
      cv_runs,
      {
        print("running cv...")
        # browser()
        param_grid = param_list[[1]]
        path_fits = create_filepath_warmstart_runs_CV(params=param_grid,
                                                      ngene=ngene, tol=TOL, 
                                                      maxit=MAXIT, nfold=NFOLD,
                                                      pkg_version=PKG_VERSION, 
                                                      git_branch=GIT_BRANCH, 
                                                      train_prefix=TRAIN_PREFIX,
                                                      method_trans_train=METHOD_TRANS_TRAIN)
        
        f = param_grid$fold
        X = data_folds$data_train[[f]]$ex
        y = data_folds$data_train[[f]]$sampInfo$time
        d = data_folds$data_train[[f]]$sampInfo$event
        
        run_warmstarts_cv(X=X,y=y,delta=d,
                          params=param_grid,
                          ninit=NINIT, alpha=ALPHA, tol=TOL, maxit=MAXIT,
                          verbose=FALSE,
                          path_fits = path_fits)
        
      },
      pattern = map(param_list),
      iteration = "list",
      format    = "file",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "cv")
      )
    ),
    
    # compile cv metrics into list
    tar_target(
      cv_metrics_list,
      {
        # browser()
        param_grid = param_list[[1]]
        bundle = readRDS(cv_runs)
        compute_metrics_bundle(data_folds,bundle,param_grid$fold)
        
      },
      pattern=map(cv_runs,param_list),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "low_mem")
      )
    ),
    
    # create dataframes of training and testing cv metrics
    tar_target(
      cv_metrics,
      {
        mets_train = dplyr::bind_rows(lapply(cv_metrics_list,`[[`,"mets_train"))
        mets_test = dplyr::bind_rows(lapply(cv_metrics_list,`[[`,"mets_test"))
        list(mets_train=mets_train,mets_test=mets_test)
      }
      
    ),
    
    # select the best parameter combo based on cv metrics
    tar_target(
      best_params,
      {
        mets_test = cv_metrics$mets_test
        avg_init = mets_test %>% group_by(k,fold,alpha,lambda,eta,lambdaW,lambdaH) %>%
          summarize(c_mean_f = mean(c), pl_mean_f=mean(sloss))%>%
          ungroup()
        avg_fold = avg_init %>% group_by(k,alpha,lambda,eta,lambdaW,lambdaH) %>%
          summarize(c_mean = mean(c_mean_f), pl_mean = mean(pl_mean_f),
                    c_sd = sqrt(sum((c_mean_f-c_mean)^2)/(NFOLD*(NFOLD-1))), 
                    pl_sd = sqrt(sum((pl_mean_f-pl_mean)^2)/(NFOLD*(NFOLD-1))))%>%
          ungroup()
        max_c = avg_fold[which.max(avg_fold$c_mean),]
        oneSE =  avg_fold %>% filter(c_mean > max_c$c_mean - max_c$c_sd)
        oneSE_top = oneSE[order(oneSE$k,oneSE$alpha),]
        oneSE_top[1,]
      }
    ),
    
    # get consensus W,H,beta for seeding at best param combo
    tar_target(
      consensus,
      {
        # browser()
        mats = collect_W_H(best_params=best_params,data_folds=data_folds,
                           ngene=ngene,tol=TOL,maxit=MAXIT,nfold=NFOLD,
                           pkg_version=PKG_VERSION,git_branch=GIT_BRANCH,
                           train_prefix=TRAIN_PREFIX,ninit=NINIT,
                           method_trans_train=METHOD_TRANS_TRAIN)
        Ws = mats$Ws
        Hs = mats$Hs
        
        res1 <- align_W_folds(Ws, ref = "auto", rows="union",tau = 0.25, 
                              consensus = "median")
        res2 <- align_W_folds(Hs, ref = "auto", rows="union",tau = 0.25, 
                              consensus = "median")
        
        W_con = res1$consensus_W
        H_con = res2$consensus_W
        
        list(H_con=H_con,W_con=W_con)
        
      }
    ),
    
    tar_target(
      fit_consensus,
      {
        H_con = consensus$H_con
        W_con = consensus$W_con
        
        genes=intersect(rownames(W_con),rownames(data_filtered$ex))
        X_con = data_filtered$ex[genes,]
        W_con = W_con[genes,]
        
        glmnet_fit = glmnet(t(X_con)%*%W_con,
                            Surv(data_filtered$sampInfo$time,
                                 data_filtered$sampInfo$event),
                            family="cox",
                            alpha=best_params$eta,
                            lambda=LAMBDA_VALS)
        beta_con = as.vector(coef(glmnet_fit,s=best_params$lambda))
        
        run_coxNMF(
          X=X_con, y=data_filtered$sampInfo$time, delta=data_filtered$sampInfo$event, 
          k=best_params$k, alpha=best_params$alpha, 
          lambda=best_params$lambda, eta=best_params$eta,
          lambdaW=0, lambdaH=0,
          tol=1e-7, maxit=6000, verbose=FALSE,
          ninit=1, imaxit=6000,
          W0=W_con, H0=t(H_con), beta0=beta_con
        )
        
      }
    ),
    
    
    #################### competing methods ###############
    
    
    ##### standard NMF package in R #####
    tar_target(
      fit_std,
      nmf(data_filtered$ex,K_VALS,nrun=NINIT,method="lee",.options=paste0("p",NINIT)),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "cv")
      )
    ),
    
    tar_target(
      std_nmf_k_selection,
      {
        path <- "data/derv/std_nmf_k_selection.csv"
        df_new <- build_std_nmf_k_selection_table(clusters_all)
        
        if (file.exists(path)) {
          # 1. make a timestamped backup of the old file
          timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
          backup_path <- paste0("std_nmf_k_selection_", timestamp, ".csv")
          file.copy(from = path, to = backup_path)
          
          message("Existing std_nmf_k_selection.csv backed up as: ", backup_path)
          
          # 2. overwrite cluster_review.csv with a clean table
          #    (i.e. you must re-curate for this new clustering)
          write.csv(df_new, path, row.names = FALSE)
          
          message("New std_nmf_k_selection.csv created from current clustering. Please re-curate.")
        } else {
          # first time: just create it
          write.csv(df_new, path, row.names = FALSE)
          message("std_nmf_k_selection.csv created. Please curate and mark keep_k / global_label.")
        }
        
        path
      },
      format = "file"
    ),
    
    
    tar_target(
      fit_std_beta,
      {
        W = fit_std$fit[[as.character(best_params$k)]]@fit@W
        H = fit_std$fit[[as.character(best_params$k)]]@fit@H
        X = data_filtered$ex
        genes = intersect(rownames(W),rownames(X))
        W=W[genes,]
        X=X[genes,]
        Z = t(X)%*%W
        y = data_filtered$sampInfo$time
        d = data_filtered$sampInfo$event
        
        gfit = list()
        mets = list()
        i=1
        for(e in ETA_VALS){
          cv_fit = cv.glmnet(Z,
                                Surv(y,d),family="cox",type.measure="C",
                                alpha=e, lambda=LAMBDA_VALS, foldid=data_folds$folds)
          
          gfit[[as.character(e)]] = cv_fit
          ind = which(cv_fit$lambda==cv_fit$lambda.1se)
          mets[[i]] = data.frame(eta = e, lambda=cv_fit$lambda.1se, cvm = cv_fit$cvm[ind])
          
          i=i+1
        }
        
        mets_all = dplyr::bind_rows(mets)
        mets_best = mets_all[which.max(mets_all$cvm),]
        
        best_fit = gfit[[as.character(mets_best$eta)]]
        beta = coef(best_fit,s=best_fit$lambda.1se)
        
        list(W = W, H = H, beta = beta)
      }
    ),
    
    
    # # lasso cox
    # tar_target(
    #   tune_enet_cox,
    #   {
    #     X = data_filtered$ex
    #     y_surv = Surv(data_filtered$sampInfo$time,
    #                  data_filtered$sampInfo$event)
    #     tune_cox_elasticnet_1se(
    #       x = t(X), y_surv = y_surv, 
    #       alpha_grid   = ETA_VALS, 
    #       lambda_grid  = LAMBDA_VALS,        # if NULL, coxnet picks its own decreasing sequence 
    #       nfolds       = 5, 
    #       type_measure = "C", 
    #       parallel     = FALSE, 
    #       standardize  = TRUE, 
    #       seed         = 1, 
    #       foldid       = data_folds$folds         # optional pre-defined fold ids (1..nfolds) 
    #     )
    #   }
    # ),
    
    # tar_target(
    #   fit_enet_cox,
    #   
    # ),
    # 
    
    
    ########## Validation Data #########
    
    tar_target(
      data_val_filtered,
      {
        genes_train = rownames(data_filtered$ex)
        setNames(
          lapply(data_val,preprocess_data_val,ngene=2500,
                 method_trans_train = METHOD_TRANS_TRAIN),
          val_datasets
        )
      }
    ),
    
    tar_target(
      data_val_comb_filtered,
      {
        genes_train = rownames(data_filtered$ex)
        preprocess_data_val(data = data_val_comb, genes = genes_train,
                            method_trans_train = METHOD_TRANS_TRAIN)
      },
    ),
    
    
    
    tar_map(
      values = tibble(ntop=NTOP),
      names = ntop,
      
      #####  get top genes #####
      
      ### DeSurv
      # get the top genes for each factor
      tar_target(
        tops_best,
        get_top_genes(W=fit_consensus$W,ntop=ntop)
      ),
      
      # visualize overlap of top genes with known lists
      tar_target(
        gene_overlap_best,
        {
          create_table(tops = tops_best$top_genes, gene_lists = top_genes,
                       which.lists = "DECODER", color.lists = colors)
        }
      ),
      
      ### NMF
      # get the top genes for each factor
      tar_target(
        tops_std,
        {
          if(length(K_VALS)>1){
            W=fit_std$fit[[as.character(best_params$k)]]@fit@W
          }else if(length(K_VALS)==1){
            W=fit_std@fit@W
          }else{
            stop("need at least one k value supplied")
          }
          
          
          get_top_genes(W=W,ntop=ntop)
        }
        
      ),
      
      # table of top gene overlap
      tar_target(
        gene_overlap_std,
        {
          create_table(tops = tops_std$top_genes, gene_lists = top_genes,
                       which.lists = "DECODER", color.lists = colors)
        }
      ),
      
      
      ######## clustering #########
      
      ### deSurv model ###
      tar_target(
        clusters_best,
        {
          data=data_val_filtered[[1]]
          val_dataset = data$dataname
          dir = create_filepath_clustering_output(ngene, TOL, MAXIT, PKG_VERSION, 
                        GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN, 
                        ntop, val_dataset, "DeSurv")
          run_clustering(tops_best$top_genes,data,top_genes,colors,
                         facs=c(1,3),plot=FALSE,dir=dir,
                         maxKcol = 5, maxKrow = 5)
        },
        pattern=map(data_val_filtered),
        iteration = "list",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "med_mem")
        )
      ),
      
      
      ### std NMF ###
      tar_target(
        clusters_std,
        {
          data=data_val_filtered[[1]]
          val_dataset = data$dataname
          dir = create_filepath_clustering_output(ngene, TOL, MAXIT, PKG_VERSION, 
                                                  GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN, 
                                                  ntop, val_dataset, "stdNMF")
          run_clustering(tops_std$top_genes,data,top_genes,colors,
                         facs=c(2,3),plot=FALSE,dir=dir,
                         maxKcol = 5, maxKrow = 5)
        },
        pattern=map(data_val_filtered),
        iteration="list",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "med_mem")
        )
      ),
      
      
      ### cox lasso ###
      
      
      ############### generate paper ###################
      tar_target(
        paper,
        {
          clusters_best
          clusters_std
          out_name =  paste0("paper_ntop=",ntop,"_ngene=",ngene,".pdf")
          rmarkdown::render("paper/paper.Rmd",knit_root_dir = "..",
                            params = list(ntop = ntop,
                                          ngene=ngene),
                            output_dir = "paper",
                            output_file = out_name)
         
        },
        format = "file",
        cue = tar_cue(mode = "always")
      )
      
      
      
    )#end map over ntop
    
    

    
  ) #end map over ngene
  

  
  ###################### supplementary material ###########################
  
  
  
  ### cold starts for paper comp
  # tar_target(param_grid_cold_cv,
  #            tidyr::expand_grid(k = c(4,6,8),
  #                               lambda = 10,
  #                               eta = 0,
  #                               lambdaW=0, 
  #                               lambdaH=0,
  #                               fold = 1:NFOLD)
  # ),
  # 
  # tar_target(
  #   cv_runs_cold,
  #   {
  #     print("running cv cold...")
  #     path_fits = create_filepath_coldstart_runs_CV(params=param_grid_cold_cv)
  #     
  #     f = param_grid_cold_cv$fold
  #     X = data_folds$data_train[[f]]$ex
  #     y = data_folds$data_train[[f]]$sampInfo$time
  #     d = data_folds$data_train[[f]]$sampInfo$event
  #     
  #     run_coldstarts_cv(X=X,y=y,delta=d,
  #                       params=param_grid_cold_cv,verbose=FALSE,
  #                       path_fits = path_fits,
  #                       ninit = NINIT_COLD)
  #     
  #   },
  #   pattern = map(param_grid_cold_cv),
  #   iteration = "list",
  #   format    = "file",
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "cv_cold")
  #   )
  # ),
  # 
  # tar_target(
  #   cv_metrics_list_cold,
  #   {
  #     bundle = readRDS(cv_runs_cold)
  #     path_mets = create_filepath_CV_metrics(params = param_grid_cold_cv)
  #     f = param_grid_cold_cv$fold
  #     k = param_grid_cold_cv$k
  #     lambda = param_grid_cold_cv$lambda
  #     eta = param_grid_cold_cv$eta
  #     lambdaW = param_grid_cold_cv$lambdaW
  #     
  #     X = data_folds$data_train[[f]]$ex
  #     y = data_folds$data_train[[f]]$sampInfo$time
  #     d = data_folds$data_train[[f]]$sampInfo$event
  #     
  #     Xtest = data_folds$data_test[[f]]$ex
  #     ytest = data_folds$data_test[[f]]$sampInfo$time
  #     dtest = data_folds$data_test[[f]]$sampInfo$event
  #     
  #     mets_tr=list()
  #     mets_te=list()
  #     j=1
  #     for(i in 1:length(bundle)){
  #       b = bundle[[i]]
  #       for(a in names(b$fits)){
  #         alpha = as.numeric(a)
  #         fit = b$fits[[a]]
  # 
  #         mets_tr[[j]] = tryCatch({
  #           m = compute_metrics(fit,X,y,d,alpha,f,test=FALSE)
  #           m$k = k
  #           m$lambda = lambda
  #           m$eta = eta
  #           m$lambdaW=lambdaW
  #           m$seed = i
  #           m
  #         },error=function(e) NULL)
  # 
  #         mets_te[[j]] = tryCatch({
  #           m = compute_metrics(fit,Xtest,ytest,dtest,alpha,f,test=TRUE)
  #           m$k = k
  #           m$lambda = lambda
  #           m$eta = eta
  #           m$lambdaW=lambdaW
  #           m$seed = i
  #           m
  #         },error=function(e) NULL)
  #         
  #         j=j+1
  #       }
  #     }
  #     mets_train = dplyr::bind_rows(mets_tr)
  #     mets_test = dplyr::bind_rows(mets_te)
  #     setNames(list(list(mets_train=mets_train, mets_test=mets_test)), 
  #              paste0(param_grid_cold_cv$k,param_grid_cold_cv$lambda,param_grid_cold_cv$eta,
  #                     param_grid_cold_cv$lambdaW,param_grid_cold_cv$lambdaH,param_grid_cold_cv$fold))
  #   },
  #   pattern=map(cv_runs_cold,param_grid_cold_cv),
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "low_mem")
  #   )
  # ),
  # 
  # tar_target(
  #   cv_metrics_cold,
  #   {
  #     mets_train = dplyr::bind_rows(lapply(cv_metrics_list_cold,`[[`,"mets_train"))
  #     mets_test = dplyr::bind_rows(lapply(cv_metrics_list_cold,`[[`,"mets_test"))
  #     list(mets_train=mets_train,mets_test=mets_test)
  #   }
  #   
  # ),
  # 
  # # bonus runs for cold and warm start full models
  # tar_target(param_grid_bonus,
  #            tidyr::expand_grid(k = c(4,6,8,10),
  #                               lambda = 10,
  #                               eta = 0,
  #                               lambdaW=0, 
  #                               lambdaH=0)
  # ),
  # 
  # tar_target(
  #   full_model_cold,
  #   {
  # 
  #     path_fits = create_filepath_coldstart_runs(param_grid_bonus)
  # 
  #     X = data_filtered$ex
  #     y = data_filtered$sampInfo$time
  #     d = data_filtered$sampInfo$event
  # 
  #     run_coldstarts(X=X,y=y,delta=d,
  #                    params=param_grid_bonus,verbose=FALSE,
  #                    path_fits = path_fits, ninit=NINIT_COLD)
  # 
  #   },
  #   format="file",
  #   pattern=map(param_grid_bonus),
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "cv_cold")
  #   )
  # ),
  # 
  # tar_target(
  #   full_metrics_cold,
  #   {
  #     bundle = readRDS(full_model_cold)
  #     path_mets = create_filepath_metrics(params = param_grid_bonus)
  #     
  #     X = data_filtered$ex
  #     y = data_filtered$sampInfo$time
  #     d = data_filtered$sampInfo$event
  #     
  #     k = param_grid_bonus$k
  #     lambda = param_grid_bonus$lambda
  #     eta = param_grid_bonus$eta
  #     lambdaW = param_grid_bonus$lambdaW
  #     
  #     mets=list()
  #     j=1
  #     for(i in 1:length(bundle)){
  #       print(i)
  #       b = bundle[[i]]
  #       for(a in names(b$fits)){
  #         alpha = as.numeric(a)
  #         fit = b$fits[[a]]
  #         
  #         mets[[j]] = tryCatch({
  #           m = compute_metrics(fit,X,y,d,alpha,NA,test=FALSE)
  #           m$k = k
  #           m$lambda = lambda
  #           m$eta = eta
  #           m$lambdaW=lambdaW
  #           m$seed = i
  #           m
  #         },error=function(e) NULL)
  #         
  #         j=j+1
  #       }
  #     }
  #     
  #     dplyr::bind_rows(mets)
  #   },
  #   pattern=map(full_model_cold,param_grid_bonus),
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "low_mem")
  #   )
  # ),
  # 
  # ## bonus warmstart runs
  # tar_target(
  #   full_model_bonus,
  #   {
  #     
  #     path_fits = create_filepath_warmstart_runs(param_grid_bonus)
  #     
  #     X = data_filtered$ex
  #     y = data_filtered$sampInfo$time
  #     d = data_filtered$sampInfo$event
  #     
  #     run_warmstarts(X=X,y=y,delta=d,
  #                    params=param_grid_bonus,verbose=FALSE,
  #                    path_fits = path_fits)
  #     
  #   },
  #   format="file",
  #   pattern = map(param_grid_bonus),
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "cv")
  #   )
  # ),
  # 
  # tar_target(
  #   full_metrics_bonus,
  #   {
  #     bundle = readRDS(full_model_bonus)
  #     path_mets = create_filepath_metrics(params = param_grid_bonus)
  #     
  #     X = data_filtered$ex
  #     y = data_filtered$sampInfo$time
  #     d = data_filtered$sampInfo$event
  #     
  #     k = param_grid_bonus$k
  #     lambda = param_grid_bonus$lambda
  #     eta = param_grid_bonus$eta
  #     lambdaW = param_grid_bonus$lambdaW
  #     mets=list()
  #     j=1
  #     for(i in 1:length(bundle)){
  #       print(i)
  #       b = bundle[[i]]
  #       for(a in names(b$fits)){
  #         alpha = as.numeric(a)
  #         fit = b$fits[[a]]
  # 
  #         mets[[j]] = tryCatch({
  #           m = compute_metrics(fit,X,y,d,alpha,NA,test=FALSE)
  #           m$k = k
  #           m$lambda = lambda
  #           m$eta = eta
  #           m$lambdaW=lambdaW
  #           m$seed = i
  #           m
  #         },error=function(e) NULL)
  #         
  #         j=j+1
  #       }
  #     }
  #     
  #     dplyr::bind_rows(mets)
  #   },
  #   pattern = map(param_grid_bonus, full_model_bonus),
  #   iteration="list",
  #   resources = tar_resources(
  #     crew = tar_resources_crew(controller = "low_mem")
  #   )
  # )
  # 
  # tar_render(
  #   paper,
  #   "paper/paper.Rmd"
  # )
  
  
  # 
  #   tar_target(
  #     best_init_per_param_combo,
  #     {
  #       df = readRDS(inits)
  #       select_best_init(df = df, method_select = METHOD_SELECT_INIT)
  #     },
  #     pattern   = map(inits),
  #     iteration = "list"
  #   ),
  #   
  #   tar_target(
  #     model_runs,
  #     {
  #       print("running coldstarts...")
  # 
  #       path = create_filepath_coldstart_runs(params = best_init_per_param_combo)
  # 
  #       run_coldstarts(
  #         X = data_filtered$ex, y = data_filtered$sampInfo$time, delta = data_filtered$sampInfo$event,
  #         params = best_init_per_param_combo,
  #         verbose = FALSE,
  #         path = path
  #       )
  # 
  #     },
  #     pattern   = map(best_init_per_param_combo),
  #     iteration = "list",
  #     format = "file",
  #     resources = tar_resources(
  #       crew = tar_resources_crew(controller = "model_runs")
  #     )
  #   ),
  #   # #
  #   # Compute model metrics
  #   tar_target(
  #     metrics,
  #     compute_metrics(model_runs,
  #                     data_filtered$ex,
  #                     data_filtered$sampInfo$time,
  #                     data_filtered$sampInfo$event),
  #     pattern   = map(model_runs),
  #     iteration = "list"
  #   ),
  #   # 
  #   tar_target(
  #     metrics_table,
  #     dplyr::bind_rows(metrics)
  #   )
  # 
  
  
  # run and save initializations for each parameter combo
  
  # 
  # # read in all initializations
  # tar_target(
  #   alpha0_scores,
  #   {
  #     df = dplyr::bind_rows(lapply(alpha0_score_file, function(path) {
  #       readRDS(path)
  #     }))
  #     df
  #   }
  # ),
  # 
  # tarchetypes::tar_group_by(
  #   alpha0_groups,                 # = grouped target name
  #   alpha0_scores,                 # = your tibble from reading *.rds
  #   groups = c(k, lambda, eta, lambdaW, lambdaH)
  # ),
  # 
  
  
  
  #
  # # initialize full model run for selected params
  
  
  ##### Summaries #####
)

