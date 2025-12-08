run_full_simulation_replicate <- function(rep_cfg) {
  simulator <- get(rep_cfg$simulator, mode = "function")
  set.seed(rep_cfg$sim_seed)
  sim <- do.call(simulator, rep_cfg$args)
  sim$surv$time <- pmax(sim$surv$time, 1e-3)
  expression_matrix <- transform_simulation_expression(
    sim$counts,
    rep_cfg$expression_transform
  )
  set.seed(rep_cfg$split_seed)
  split <- split_train_test(ncol(expression_matrix), train_frac = rep_cfg$train_fraction %||% 0.7)
  train_dataset <- build_simulation_dataset(
    expression_matrix = expression_matrix,
    survival_df = sim$surv,
    column_indices = split$train,
    dataname = paste0(rep_cfg$run_id, "_train")
  )
  test_dataset <- build_simulation_dataset(
    expression_matrix = expression_matrix,
    survival_df = sim$surv,
    column_indices = split$test,
    dataname = paste0(rep_cfg$run_id, "_test")
  )

  bo_summary <- run_single_bo(
    data = train_dataset,
    run_prefix = rep_cfg$run_id,
    bo_bounds = rep_cfg$bo_bounds,
    bo_extra_args = rep_cfg$bo_extra_args,
    hyper_defaults = rep_cfg$hyper_defaults,
    coarse_control = rep_cfg$bo_coarse_control,
    refine_control = rep_cfg$bo_refine_control,
    max_refinements = rep_cfg$bo_max_refinements,
    tol_gain = rep_cfg$bo_tol_gain,
    plateau = rep_cfg$bo_plateau,
    top_k = rep_cfg$bo_top_k,
    shrink_base = rep_cfg$bo_shrink_base,
    importance_gain = rep_cfg$bo_importance_gain,
    parallel_grid = rep_cfg$parallel_grid
  )
  alpha0_bounds <- rep_cfg$bo_bounds
  alpha0_bounds$alpha_grid <- NULL
  alpha0_extra_args <- rep_cfg$bo_extra_args %||% list()
  alpha0_extra_args$alpha_grid <- 0
  bo_summary_alpha0 <- run_single_bo(
    data = train_dataset,
    run_prefix = paste0(rep_cfg$run_id, "_alpha0"),
    bo_bounds = alpha0_bounds,
    bo_extra_args = alpha0_extra_args,
    hyper_defaults = rep_cfg$hyper_defaults,
    coarse_control = rep_cfg$bo_coarse_control,
    refine_control = rep_cfg$bo_refine_control,
    max_refinements = rep_cfg$bo_max_refinements,
    tol_gain = rep_cfg$bo_tol_gain,
    plateau = rep_cfg$bo_plateau,
    top_k = rep_cfg$bo_top_k,
    shrink_base = rep_cfg$bo_shrink_base,
    importance_gain = rep_cfg$bo_importance_gain,
    parallel_grid = rep_cfg$parallel_grid
  )
  params_best <- bo_summary$params_best
  ngene_value <- extract_param_value(params_best$ngene, NGENE_DEFAULT, integer = TRUE)
  ntop_value <- extract_param_value(params_best$ntop, NTOP_DEFAULT, integer = TRUE)
  lambdaW_value <- extract_param_value(params_best$lambdaW, LAMBDAW_DEFAULT)
  lambdaH_value <- extract_param_value(params_best$lambdaH, LAMBDAH_DEFAULT)
  params_best_alpha0 <- bo_summary_alpha0$params_best
  params_best_alpha0$alpha <- 0
  ngene_value_alpha0 <- extract_param_value(params_best_alpha0$ngene, NGENE_DEFAULT, integer = TRUE)
  ntop_value_alpha0 <- extract_param_value(params_best_alpha0$ntop, NTOP_DEFAULT, integer = TRUE)
  lambdaW_value_alpha0 <- extract_param_value(params_best_alpha0$lambdaW, LAMBDAW_DEFAULT)
  lambdaH_value_alpha0 <- extract_param_value(params_best_alpha0$lambdaH, LAMBDAH_DEFAULT)

  data_filtered <- DeSurv::preprocess_data(
    X = train_dataset$ex,
    y = train_dataset$sampInfo$time,
    d = train_dataset$sampInfo$event,
    dataset = train_dataset$sampInfo$dataset,
    samp_keeps = train_dataset$samp_keeps,
    ngene = ngene_value,
    method_trans_train = METHOD_TRANS_TRAIN,
    verbose = FALSE
  )
  data_filtered_alpha0 <- DeSurv::preprocess_data(
    X = train_dataset$ex,
    y = train_dataset$sampInfo$time,
    d = train_dataset$sampInfo$event,
    dataset = train_dataset$sampInfo$dataset,
    samp_keeps = train_dataset$samp_keeps,
    ngene = ngene_value_alpha0,
    method_trans_train = METHOD_TRANS_TRAIN,
    verbose = FALSE
  )

  seed_fits <- run_seed_fits(data_filtered, params_best, lambdaW_value, lambdaH_value)
  consensus <- DeSurv::desurv_consensus_seed(
    fits = seed_fits$fits,
    X = data_filtered$ex,
    ntop = ntop_value,
    k = params_best$k,
    min_frequency = max(1, ceiling(0.3 * length(seed_fits$fits)))
  )
  seed_fits_alpha0 <- run_seed_fits(data_filtered_alpha0, params_best_alpha0, lambdaW_value_alpha0, lambdaH_value_alpha0)
  consensus_alpha0 <- DeSurv::desurv_consensus_seed(
    fits = seed_fits_alpha0$fits,
    X = data_filtered_alpha0$ex,
    ntop = ntop_value_alpha0,
    k = params_best_alpha0$k,
    min_frequency = max(1, ceiling(0.3 * length(seed_fits_alpha0$fits)))
  )

  final_fit <- desurv_fit(
    X = data_filtered$ex,
    y = data_filtered$sampInfo$time,
    d = data_filtered$sampInfo$event,
    k = params_best$k,
    alpha = params_best$alpha,
    lambda = params_best$lambda,
    nu = params_best$nu,
    lambdaW = lambdaW_value,
    lambdaH = lambdaH_value,
    W0 = consensus$W0,
    H0 = consensus$H0,
    beta0 = consensus$beta0,
    tol = TOL / 100,
    maxit = MAXIT,
    verbose = FALSE
  )
  final_fit_alpha0 <- desurv_fit(
    X = data_filtered_alpha0$ex,
    y = data_filtered_alpha0$sampInfo$time,
    d = data_filtered_alpha0$sampInfo$event,
    k = params_best_alpha0$k,
    alpha = params_best_alpha0$alpha,
    lambda = params_best_alpha0$lambda,
    nu = params_best_alpha0$nu,
    lambdaW = lambdaW_value_alpha0,
    lambdaH = lambdaH_value_alpha0,
    W0 = consensus_alpha0$W0,
    H0 = consensus_alpha0$H0,
    beta0 = consensus_alpha0$beta0,
    tol = TOL / 100,
    maxit = MAXIT,
    verbose = FALSE
  )

  test_processed <- DeSurv::preprocess_data(
    X = test_dataset$ex,
    y = test_dataset$sampInfo$time,
    d = test_dataset$sampInfo$event,
    dataset = test_dataset$sampInfo$dataset,
    samp_keeps = test_dataset$samp_keeps,
    genes = data_filtered$featInfo,
    method_trans_train = METHOD_TRANS_TRAIN,
    verbose = FALSE
  )
  test_processed$dataname <- test_dataset$dataname
  if (is.null(test_processed$sampInfo$ID)) {
    test_processed$sampInfo$ID <- colnames(test_processed$ex)
  }
  rownames(test_processed$sampInfo) <- test_processed$sampInfo$ID
  if (is.null(test_processed$sampInfo$dataset)) {
    test_processed$sampInfo$dataset <- test_dataset$sampInfo$dataset
  }

  truth <- list(
    counts = sim$counts,
    expression = expression_matrix,
    surv = sim$surv,
    W_true = sim$W_true,
    H_true = sim$H_true,
    beta_true = sim$beta_true,
    marker_info = sim$marker_info,
    config_name = rep_cfg$config_name
  )
  metrics <- tryCatch(
    compute_test_metrics(
      truth,
      test_processed,
      final_fit
    ),
    error = function(e) default_metrics()
  )
  metrics_alpha0 <- tryCatch(
    compute_test_metrics(
      truth,
      test_processed,
      final_fit_alpha0
    ),
    error = function(e) default_metrics()
  )
  model_results <- list(
    list(
      model = "desurv",
      metrics = metrics,
      params = params_best,
      bo_best_score = bo_summary$best_score,
      bo_history = bo_summary$history_path
    ),
    list(
      model = "desurv_alpha0",
      metrics = metrics_alpha0,
      params = params_best_alpha0,
      bo_best_score = bo_summary_alpha0$best_score,
      bo_history = bo_summary_alpha0$history_path
    )
  )

  list(
    config_name = rep_cfg$config_name,
    config_signature = rep_cfg$config_signature,
    design_name = rep_cfg$design_name,
    replicate = rep_cfg$replicate,
    split = split,
    params = params_best,
    bo_history = bo_summary$history_path,
    bo_best_score = bo_summary$best_score,
    metrics = metrics,
    test_dataset = test_dataset$dataname,
    alpha0 = list(
      params = params_best_alpha0,
      bo_history = bo_summary_alpha0$history_path,
      bo_best_score = bo_summary_alpha0$best_score,
      metrics = metrics_alpha0
    ),
    model_results = model_results
  )
}
