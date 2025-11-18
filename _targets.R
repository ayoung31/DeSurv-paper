# _targets.R
library(targets)
library(tarchetypes)  
library(crew.cluster)
library(crew)
suppressWarnings(suppressMessages(library(dplyr)))

if (dir.exists("../DeSurv") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressMessages(
    pkgload::load_all(
      path = "../DeSurv",
      export_all = FALSE,
      helpers = FALSE,
      attach_testthat = FALSE
    )
  )
}

NINIT=20
NINIT_COLD = 5
NINIT_FULL = 100
BO_N_INIT = 12
BO_N_ITER = 30
BO_CANDIDATE_POOL = 2000

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
  packages = c("DeSurv","pheatmap","NMF","tidyverse","tidyselect","survival","cvwrapr","rmarkdown","dplyr",
               "parallel","foreach", "doParallel", "doMC", "pec", "glmnet","webshot2"),
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
NGENE_CONFIG       = c(2000,5000) # set length 1 to fix, or range to tune via BO
NGENE_DEFAULT      = NGENE_CONFIG[[1]]
TUNE_NGENE         = length(unique(NGENE_CONFIG)) > 1
NTOP_CONFIG        = c(50) # set length 1 to fix or range to tune via BO
NTOP_DEFAULT       = NTOP_CONFIG[[1]]
TUNE_NTOP          = length(unique(NTOP_CONFIG)) > 1
LAMBDAW_CONFIG     = c(0)
LAMBDAW_DEFAULT    = LAMBDAW_CONFIG[[1]]
TUNE_LAMBDAW       = length(unique(LAMBDAW_CONFIG)) > 1
LAMBDAH_CONFIG     = c(0)
LAMBDAH_DEFAULT    = LAMBDAH_CONFIG[[1]]
TUNE_LAMBDAH       = length(unique(LAMBDAH_CONFIG)) > 1
DESURV_BO_BOUNDS   = list(
  k_grid = list(lower = 2L, upper = 12L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
  lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)
TOL                = 1e-4
MAXIT              = 6000
STD_NMF_K_GRID     = 2:12#2:16   #= c(2,3,4,5)
COXNET_LAMBDA_GRID = c(1e-4,1e-3,1e-2,.1,1,10)#,seq(.2,.9,by=.1)   #10^seq(-3,3)#10^seq(-4,4)
COXNET_ALPHA_GRID  = seq(0,1,by=.1)#c(0,.01,.1,.5,.9)#c(0,.01)#c(.01,.1,.5,.9)#seq(.1,.9,by=.1)
NFOLD              = 5
PKG_VERSION        = utils::packageDescription("DeSurv", fields = "RemoteRef")
GIT_BRANCH         = gert::git_branch()

# load top_genes and colors into global environment 
# (known gene lists and corresponding color schemes)
load("data/derv/cmbSubtypes_formatted.RData")

# ---- Targets ----
list(
  
  ################# Training ###################
  
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
  # Preprocess and filter data
  # Split data into folds, preprocess and filter each fold
  # Bayesian optimisation of hyperparameters
  tar_target(
    desurv_bo_results,
    {
      bounds <- DESURV_BO_BOUNDS
      bounds <- maybe_add_numeric_bound(bounds, NGENE_CONFIG, "ngene", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, NTOP_CONFIG, "ntop", type = "integer")
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAW_CONFIG, "lambdaW_grid", log_scale = TRUE)
      bounds <- maybe_add_numeric_bound(bounds, LAMBDAH_CONFIG, "lambdaH_grid", log_scale = TRUE)

      extra_args <- list(n_starts = NINIT)
      if (!TUNE_NGENE) {
        extra_args$ngene <- NGENE_DEFAULT
      }
      if (!TUNE_NTOP) {
        extra_args$ntop <- NTOP_DEFAULT
      }
      if (!TUNE_LAMBDAW) {
        extra_args$lambdaW_grid <- LAMBDAW_DEFAULT
      }
      if (!TUNE_LAMBDAH) {
        extra_args$lambdaH_grid <- LAMBDAH_DEFAULT
      }

      do.call(
        DeSurv::desurv_cv_bayesopt,
        c(
          list(
            X = data$ex,
            y = data$sampInfo$time,
            d = data$sampInfo$event,
            dataset = data$sampInfo$dataset,
            samp_keeps = data$samp_keeps,
            preprocess = TRUE,
            method_trans_train = METHOD_TRANS_TRAIN,
            engine = "warmstart",
            nfolds = NFOLD,
            tol = TOL,
            maxit = MAXIT,
            bo_bounds = bounds,
            n_init = BO_N_INIT,
            n_iter = BO_N_ITER,
            candidate_pool = BO_CANDIDATE_POOL,
            exploration_weight = 0.01,
            seed = 123,
            cv_verbose = FALSE,
            verbose = TRUE
          ),
          extra_args
        )
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  
  tar_target(
    params_best,
    standardize_bo_params(desurv_bo_results$best$params)
  ),
  
  tar_target(
    ngene_value,
    {
      value <- params_best$ngene
      if (is.null(value) || is.na(value)) {
        as.integer(NGENE_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    ntop_value,
    {
      value <- params_best$ntop
      if (is.null(value) || is.na(value)) {
        as.integer(NTOP_DEFAULT)
      } else {
        as.integer(round(value))
      }
    }
  ),

  tar_target(
    lambdaW_value,
    {
      value <- params_best$lambdaW
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAW_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),

  tar_target(
    lambdaH_value,
    {
      value <- params_best$lambdaH
      if (is.null(value) || is.na(value)) {
        as.numeric(LAMBDAH_DEFAULT)
      } else {
        as.numeric(value)
      }
    }
  ),
  
  tar_target(
    training_results_dir,
    results_root_dir(
      ngene = ngene_value,
      tol = TOL,
      maxit = MAXIT,
      pkg_version = PKG_VERSION,
      git_branch = GIT_BRANCH,
      train_prefix = TRAIN_PREFIX,
      method_trans_train = METHOD_TRANS_TRAIN
    )
  ),
  
  tar_target(
    desurv_bo_history,
    {
      path <- file.path(training_results_dir, "desurv_bo_history.csv")
      utils::write.csv(desurv_bo_results$history, path, row.names = FALSE)
      path
    },
    format = "file"
  ),
  
  tar_target(
    data_filtered,
    {
      prep <- DeSurv::preprocess_data(
        X = data$ex,
        y = data$sampInfo$time,
        d = data$sampInfo$event,
        dataset = data$sampInfo$dataset,
        samp_keeps = data$samp_keeps,
        ngene = ngene_value,
        method_trans_train = METHOD_TRANS_TRAIN,
        verbose = FALSE
      )
      prep$dataname <- data$dataname
      prep
    }
  ),
  
  tar_target(
    fit_desurv,
    {
      seeds <- seq_len(NINIT_FULL)
      fits <- vector("list", length(seeds))
      scores <- rep(NA_real_, length(seeds))
      for (i in seq_along(seeds)) {
        fits[[i]] <- try(
          desurv_fit(
            X = data_filtered$ex,
            y = data_filtered$sampInfo$time,
            d = data_filtered$sampInfo$event,
            k = params_best$k,
            alpha = params_best$alpha,
            lambda = params_best$lambda,
            nu = params_best$nu,
            lambdaW = lambdaW_value,
            lambdaH = lambdaH_value,
            seed = seeds[i],
            tol = TOL / 100,
            tol_init = TOL,
            maxit = MAXIT,
            imaxit = MAXIT,
            ninit = 1,
            parallel_init = FALSE,
            verbose = FALSE
          ),
          silent = TRUE
        )
        if (!inherits(fits[[i]], "try-error") && !is.null(fits[[i]]$cindex)) {
          scores[i] <- fits[[i]]$cindex
        }
      }
      best_idx <- which.max(scores)
      if (!length(best_idx) || is.na(scores[best_idx])) {
        stop("No successful full-model fits were obtained.")
      }
      fits[[best_idx]]
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = "full")
    )
  ),
  
  
  #################### competing methods ###############
  
  
  ##### standard NMF package in R #####
  tar_target(
    fit_std,
    nmf(data_filtered$ex,STD_NMF_K_GRID,nrun=NINIT,method="lee",.options=paste0("p",NINIT)),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  
  tar_target(
    std_nmf_k_selection_plots,
    {
      save_dir=file.path(training_results_dir,"std_nmf_k_selection")
      dir.create(save_dir,showWarnings = FALSE)
      path = file.path(save_dir,paste0("std_nmf_k_selection_plots.png"))
      if (file.exists(path)) {
        # 1. make a timestamped backup of the old file
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        backup_path <- file.path(dirname(path),
                                 paste0("std_nmf_k_selection_plots_", timestamp, ".png")
        )
        file.copy(from = path, to = backup_path)
      }
      png(filename=path)
      print(NMF::plot(fit_std))
      dev.off()
      path
    },
    format="file"
  ),
  
  tar_target(
    std_nmf_k_selection_table,
    {
      save_dir=dirname(std_nmf_k_selection_plots)
      path <- file.path(save_dir,
                        paste0("std_nmf_k_selection_table.csv")
                        )
      if(!file.exists(std_nmf_k_selection_plots)){
        stop("Standard NMF k selection plots not generated")
      }
      build_std_nmf_k_selection_table(STD_NMF_K_GRID,path,fit_std)
    },
  ),
  
  tar_target(
    fit_std_selected_k,
    {
      std_nmf_k_table = read.csv(std_nmf_k_selection_table,
                                 stringsAsFactors = FALSE)
      if(all(std_nmf_k_table$selected==FALSE)){
        stop("Please edit std_nmf_k_selection.csv to select a k")
      }
      k_selected = std_nmf_k_table$rank[which(std_nmf_k_table$selected)]
      print(paste0("You have selected ",k_selected," factors for standard NMF"))
      fit_std$fit[[as.character(k_selected)]]
    }
  ),
  
  
  tar_target(
    fit_std_beta,
    {
      W = fit_std_selected_k@fit@W
      H = fit_std_selected_k@fit@H
      X = data_filtered$ex
      genes = intersect(rownames(W),rownames(X))
      W=W[genes,]
      X=X[genes,]
      Z = t(X)%*%W
      y = data_filtered$sampInfo$time
      d = data_filtered$sampInfo$event
      strata <- interaction(d, data_filtered$sampInfo$dataset, drop = FALSE)
      foldid <- caret::createFolds(strata, NFOLD, list = FALSE)
      
      gfit = list()
      mets = list()
      i=1
      for(alpha in COXNET_ALPHA_GRID){
        cv_fit = cv.glmnet(Z,
                           Surv(y,d),family="cox",type.measure="C",
                           alpha=alpha, lambda=COXNET_LAMBDA_GRID, foldid=foldid)
        
        gfit[[as.character(alpha)]] = cv_fit
        ind = which(cv_fit$lambda==cv_fit$lambda.1se)
        mets[[i]] = data.frame(eta = alpha, lambda=cv_fit$lambda.1se, cvm = cv_fit$cvm[ind])
        
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
  #       alpha_grid   = COXNET_ALPHA_GRID, 
  #       lambda_grid  = COXNET_LAMBDA_GRID,        # if NULL, coxnet picks its own decreasing sequence 
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
      setNames(
        lapply(
          data_val,
          function(dataset) {
            prep <- DeSurv::preprocess_data(
              X = dataset$ex,
              y = dataset$sampInfo$time,
              d = dataset$sampInfo$event,
              dataset = dataset$sampInfo$dataset,
              samp_keeps = dataset$samp_keeps,
              ngene = ngene_value,
              method_trans_train = METHOD_TRANS_TRAIN,
              verbose = FALSE
            )
            prep$dataname <- dataset$dataname
            prep
          }
        ),
        val_datasets
      )
    }
  ),
  
  tar_target(
    data_val_comb_filtered,
    {
      genes_train = rownames(data_filtered$ex)
      prep <- DeSurv::preprocess_data(
        X = data_val_comb$ex,
        y = data_val_comb$sampInfo$time,
        d = data_val_comb$sampInfo$event,
        dataset = data_val_comb$sampInfo$dataset,
        samp_keeps = data_val_comb$samp_keeps,
        genes = genes_train,
        method_trans_train = METHOD_TRANS_TRAIN,
        verbose = FALSE
      )
      prep$dataname <- data_val_comb$dataname
      prep
    },
  ),
  
  
  #####  get top genes #####
  
  ### DeSurv
  # get the top genes for each factor
  tar_target(
    tops_desurv,
    get_top_genes(W = fit_desurv$W, ntop = ntop_value)
  ),
  
  # visualize overlap of top genes with known lists
  tar_target(
    gene_overlap_desurv,
    {
      create_table(tops = tops_desurv$top_genes, gene_lists = top_genes,
                   which.lists = "DECODER", color.lists = colors)
    }
  ),
  
  ### NMF
  # get the top genes for each factor
  tar_target(
    tops_std,
    {
      W = fit_std_selected_k@fit@W
      get_top_genes(W = W, ntop = ntop_value)
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
  
  #### select factors ####
  tar_target(
    selected_factors_desurv,
    {
      save_dir = file.path(training_results_dir,"factor_selection","desurv")
      dir.create(save_dir,showWarnings = FALSE,recursive = TRUE)
      path = file.path(save_dir,"factor_selection.csv")
      build_factor_selection_table(tops_desurv$top_genes,path)
    },
    format="file"
  ),

  tar_target(
    selected_factors_std,
    {
      save_dir = file.path(training_results_dir,"factor_selection","std")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"factor_selection.csv")
      build_factor_selection_table(tops_std$top_genes,path)
    },
    format="file"
  ),
  # 
  ######## clustering #########
  
  ### deSurv model ###
    tar_target(
      clusters_desurv,
      {
        if(!file.exists(selected_factors_desurv)){
          stop("first generate factor selection table for desurv")
        }
        tbl = read.csv(selected_factors_desurv)
        sel = tbl$factor[tbl$selected]
        if(all(!sel)){
          stop("No factor selected for desurv clustering, please select at least 1")
        }
        data <- data_val_filtered
        val_dataset = data$dataname
        dir = create_filepath_clustering_output(ngene_value, TOL, MAXIT, PKG_VERSION, 
                                                GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN, 
                                                ntop_value, val_dataset, "DeSurv")
        run_clustering(tops_desurv$top_genes,data,top_genes,colors,
                       facs=sel,plot=FALSE,dir=dir,
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
        if(!file.exists(selected_factors_std)){
          stop("first generate factor selection table for std nmf")
        }
        tbl = read.csv(selected_factors_std)
        sel = tbl$factor[tbl$selected]
        if(all(!sel)){
          stop("No factor selected for std nmf clustering, please select at least 1")
        }
        data <- data_val_filtered
        val_dataset = data$dataname
        dir = create_filepath_clustering_output(ngene_value, TOL, MAXIT, PKG_VERSION, 
                                                GIT_BRANCH, TRAIN_PREFIX, METHOD_TRANS_TRAIN, 
                                                ntop_value, val_dataset, "stdNMF")
        run_clustering(tops_std$top_genes,data,top_genes,colors,
                       facs=sel,plot=FALSE,dir=dir,
                     maxKcol = 5, maxKrow = 5)
    },
    pattern=map(data_val_filtered),
    iteration="list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "med_mem")
    )
  ),
  
  tar_target(
    selected_nclusters_desurv,
    {
      save_dir = file.path(training_results_dir,"ncluster_selection","desurv")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"ncluster_selection.csv")
      build_ncluster_selection_table(clusters_desurv,path)
    },
    format="file"
  ),
  
  tar_target(
    selected_nclusters_std,
    {
      save_dir = file.path(training_results_dir,"ncluster_selection","std")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"ncluster_selection.csv")
      build_ncluster_selection_table(clusters_std,path)
    },
    format="file"
  ),
  
  tar_target(
    cluster_alignment_plot_desurv,
    {
      nclus_tbl = read.csv(selected_nclusters_desurv)
      nclus_vec <- setNames(nclus_tbl$nclus, nclus_tbl$dataset)
      save_dir = file.path(training_results_dir,"cluster_alignment","desurv")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"cluster_alignment.pdf")
      pdf(path)
      plot_cluster_cor(clusters_desurv,tops_desurv$top_genes,nclus_vec)
      dev.off()
    }
  ),
  
  tar_target(
    cluster_alignment_plot_std,
    {
      nclus_tbl = read.csv(selected_nclusters_std)
      nclus_vec <- setNames(nclus_tbl$nclus, nclus_tbl$dataset)
      save_dir = file.path(training_results_dir,"cluster_alignment","std")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"cluster_alignment.pdf")
      pdf(path)
      plot_cluster_cor(clusters_std,tops_std$top_genes,nclus_vec)
      dev.off()
    }
  ), 
  
  tar_target(
    cluster_alignment_table_desurv,
    {
      nclus_tbl = read.csv(selected_nclusters_desurv)
      save_dir = file.path(training_results_dir,"cluster_alignment","desurv")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"cluster_alignment.csv")
      build_cluster_alignment_table(nclus_tbl,clusters_desurv,path)
    },
    format="file"
  ),
  
  tar_target(
    cluster_alignment_table_std,
    {
      nclus_tbl = read.csv(selected_nclusters_std)
      save_dir = file.path(training_results_dir,"cluster_alignment","std")
      dir.create(save_dir,showWarnings = FALSE, recursive=TRUE)
      path = file.path(save_dir,"cluster_alignment.csv")
      build_cluster_alignment_table(nclus_tbl,clusters_std,path)
    },
    format="file"
  ),
  
  tar_target(
    aligned_clusters_desurv,
    {
      nclus_tbl = read.csv(selected_nclusters_desurv)
      samp_clus = read.csv(cluster_alignment_table_desurv)
      clus=clusters_desurv
      datasets <- vapply(clus, function(x) x$data$dataname, character(1))
      for(i in seq_along(clus)){
        dataname=datasets[[i]]
        nclus_i <- nclus_tbl$nclus[nclus_tbl$dataset==dataname]
        if(length(nclus_i) != 1 || is.na(nclus_i)){
          stop("Selected ncluster table missing entry for dataset: ", dataname)
        }
        nms = names(clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass)
        y = samp_clus[clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass,dataname]
        names(y) = nms
        clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass=y

        clus[[i]]$data$sampInfo$samp_cluster = clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass
      }
      clus
    }
  ),
  
  tar_target(
    aligned_clusters_std,
    {
      nclus_tbl = read.csv(selected_nclusters_std)
      samp_clus = read.csv(cluster_alignment_table_std)
       clus=clusters_std
      datasets <- vapply(clus, function(x) x$data$dataname, character(1))
      for(i in seq_along(clus)){
        dataname=datasets[[i]]
        nclus_i <- nclus_tbl$nclus[nclus_tbl$dataset==dataname]
        if(length(nclus_i) != 1 || is.na(nclus_i)){
          stop("Selected ncluster table missing entry for dataset: ", dataname)
        }
        nms = names(clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass)
        y = samp_clus[clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass,dataname]
        names(y) = nms
        clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass=y
        
        clus[[i]]$data$sampInfo$samp_cluster = clus[[i]]$clus_res$clusCol[[nclus_i]]$consensusClass
      }
      clus
    }
  )
  
    
    ### cox lasso ###
    
    
    ############### generate paper ###################
    # tar_target(
    #   paper,
    #   {
    #     clusters_desurv
    #     clusters_std
    #     out_name =  paste0("paper_ntop=",ntop_value,"_ngene=",ngene_value,".pdf")
    #     rmarkdown::render("paper/paper.Rmd",knit_root_dir = "..",
    #                       params = list(ntop = ntop_value,
    #                                     ngene=ngene_value),
    #                       output_dir = "paper",
    #                       output_file = out_name)
    #     
    #   },
    #   format = "file",
    #   cue = tar_cue(mode = "always")
    # )
    # 
  
  
  
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
