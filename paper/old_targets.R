# 
# # run multiple seeds of warmstart for selected param combo on full training data
# tar_target(
#   full_model,
#   {
#     
#     path_fits = create_filepath_warmstart_runs(best_params, ngene=ngene, 
#                                                tol=TOL, maxit=MAXIT, 
#                                                pkg_version=PKG_VERSION, 
#                                                git_branch=GIT_BRANCH, 
#                                                train_prefix=TRAIN_PREFIX,
#                                                method_trans_train=METHOD_TRANS_TRAIN)
#     
#     X = data_filtered$ex
#     y = data_filtered$sampInfo$time
#     d = data_filtered$sampInfo$event
#     
#     run_warmstarts(X=X,y=y,delta=d,
#                    params=best_params,
#                    ninit=NINIT, alpha=ALPHA, tol=TOL, maxit=MAXIT,
#                    verbose=FALSE,
#                    path_fits = path_fits)
#     
#   },
#   format="file",
#   resources = tar_resources(
#     crew = tar_resources_crew(controller = "full")
#   )
# ),
# 
# # compile metrics for full data run 
# tar_target(
#   full_metrics,
#   {
#     bundle = readRDS(full_model)
#     path_mets = create_filepath_metrics(params = best_params, ngene=ngene, 
#                                         tol=TOL, maxit=MAXIT, nfold=NFOLD, 
#                                         pkg_version=PKG_VERSION, 
#                                         git_branch=GIT_BRANCH, 
#                                         train_prefix=TRAIN_PREFIX, 
#                                         method_trans_train=METHOD_TRANS_TRAIN)
#     
#     k = best_params$k
#     lambda = best_params$lambda
#     eta = best_params$eta
#     lambdaW = best_params$lambdaW
#     lambdaH = best_params$lambdaH
#     mets=list()
#     j=1
#     for(i in 1:length(bundle)){
#       print(i)
#       b = bundle[[i]]
#       for(a in names(b$fits)){
#         alpha = as.numeric(a)
#         fit = b$fits[[a]]
#         
#         X = data_filtered$ex
#         y = data_filtered$sampInfo$time
#         d = data_filtered$sampInfo$event
#         
#         mets[[j]] = tryCatch({
#           m = compute_metrics(fit,X,y,d,test=FALSE)
#           m$alpha = alpha
#           m$k = k
#           m$lambda = lambda
#           m$eta = eta
#           m$lambdaW=lambdaW
#           m$lambdaH=lambdaH
#           m$seed = i
#           m
#         },error=function(e) NULL)
#         
#         j=j+1
#       }
#     }
#     
#     dplyr::bind_rows(mets)
#   }
# ), 
# 
# 
# # reduce to selected alpha and best seed
# tar_target(
#   fit_train,
#   {
#     alpha_best = full_metrics %>% filter(alpha==best_params$alpha)
#     init_best = select_best_init(alpha_best,method_select=METHOD_SELECT_INIT)
#     fit_all = readRDS(full_model)
#     fit_init = fit_all[[init_best$seed]]
#     meta = fit_init$meta
#     meta$alpha = init_best$alpha
#     fit = fit_init$fits[[as.character(init_best$alpha)]]
#     fit$meta = meta
#     
#     fit
#   }
# ),
#


##### alpha = 0 #####
# tar_target(
#   fit_a0,
#   {
#     alpha0 = full_metrics %>% filter(alpha==0)
#     init_best = select_best_init(alpha0,method_select="loss")
#     # same k as selected model
#     fit_all = readRDS(full_model)
#     fit_init = fit_all[[init_best$seed]]
#     meta = fit_init$meta
#     meta$alpha = 0
#     fit = fit_init$fits[[as.character(0)]]
#     fit$meta = meta
#     
#     fit
#   }
# ),
# 
# # get the top genes for each factor
# tar_target(
#   tops_a0,
#   get_top_genes(W=fit_a0$W,ntop=NTOP)
# ),
# 
# # table of top gene overlap
# tar_target(
#   gene_overlap_a0,
#   {
#     known_env = new.env()
#     load("data/derv/cmbSubtypes_formatted.RData",envir = known_env)
#     create_table(tops = tops_a0$top_genes, gene_lists = known_env$top_genes,
#                  which.lists = "DECODER", color.lists = known_env$colors)
#   }
# ),



# Load and format raw data
# tar_target(
#   data_val,
#   {
#     raw_data_val
#     load_data(VAL_DATASETS)
#   }
# ),
# 
# # Preprocess and filter data
# tar_target(data_val_filtered,
#   {
#     train_genes = rownames(fit_consensus$W)
#     preprocess_data(data = data_val,
#                     genes = train_genes,
#                     method_trans_train = METHOD_TRANS_TRAIN)
#   }
# ),