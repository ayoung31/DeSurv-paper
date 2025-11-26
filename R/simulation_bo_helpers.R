`%||%` <- function(x, y) if (is.null(x)) y else x

format_simulation_args <- function(args_list) {
  if (!length(args_list)) {
    return("none")
  }
  parts <- vapply(
    names(args_list),
    function(name) {
      value <- args_list[[name]]
      paste0(name, "=", paste(value, collapse = ","))
    },
    character(1)
  )
  paste(parts, collapse = ";")
}

merge_simulation_args <- function(base_args, override_args) {
  base_args <- base_args %||% list()
  override_args <- override_args %||% list()
  if (!length(base_args)) {
    return(override_args)
  }
  if (!length(override_args)) {
    return(base_args)
  }
  modifyList(base_args, override_args, keep.null = TRUE)
}

standardize_design_spec <- function(design_spec,
                                    provided_name = NULL,
                                    fallback_label = NULL,
                                    index = 1L) {
  if (is.null(design_spec)) {
    design_spec <- list()
  }
  if (!is.list(design_spec)) {
    stop("Simulation design specifications must be lists.")
  }
  reserved <- c("label", "design_label", "design_name", "name", "args", "params")
  label <- design_spec$label %||% design_spec$design_label %||%
    design_spec$design_name %||% provided_name %||%
    fallback_label %||% sprintf("design_%02d", index)
  params <- design_spec$params %||% design_spec$args
  if (is.null(params)) {
    keep <- setdiff(names(design_spec), reserved)
    if (length(keep)) {
      params <- design_spec[keep]
    } else {
      params <- list()
    }
  }
  if (!is.list(params)) {
    stop("Simulation design parameter sets must be lists.")
  }
  list(label = label, params = params)
}

extract_simulation_designs <- function(cfg, scenario_name) {
  base_args <- cfg$args %||% list()
  design_list <- cfg$designs
  explicit <- !is.null(design_list)
  if (is.null(design_list) && !is.null(cfg$design)) {
    design_list <- list(cfg$design)
    explicit <- TRUE
  }
  if (is.null(design_list)) {
    params <- cfg$design_params %||% list()
    args <- merge_simulation_args(base_args, params)
    return(list(list(
      label = cfg$design_name,
      params = params,
      args = args,
      explicit = FALSE
    )))
  }
  if (!is.list(design_list)) {
    stop(
      "`designs` for scenario `", scenario_name,
      "` must be a list of design specifications."
    )
  }
  if (!length(design_list)) {
    design_list <- list(list())
  }
  idx <- 0L
  purrr::imap(
    design_list,
    function(spec, name) {
      idx <<- idx + 1L
      normalized <- standardize_design_spec(
        spec,
        provided_name = if (!is.null(name) && nzchar(name)) name else NULL,
        fallback_label = cfg$design_name,
        index = idx
      )
      list(
        label = normalized$label,
        params = normalized$params,
        args = merge_simulation_args(base_args, normalized$params),
        explicit = explicit
      )
    }
  )
}

transform_simulation_expression <- function(counts, transform_spec = "log2") {
  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  }
  if (is.character(transform_spec)) {
    transform_spec <- match.arg(transform_spec, c("log2", "rank", "none"))
    expr <- switch(
      transform_spec,
      log2 = log2(counts + 1),
      rank = apply(counts, 2, rank, ties.method = "average"),
      none = counts
    )
  } else if (is.function(transform_spec)) {
    expr <- transform_spec(counts)
  } else {
    stop("Unsupported transform specification for simulation expression.")
  }
  expr <- as.matrix(expr)
  dimnames(expr) <- dimnames(counts)
  expr
}

build_simulation_dataset <- function(expression_matrix,
                                     survival_df,
                                     column_indices,
                                     dataname) {
  if (!length(column_indices)) {
    stop("Simulation split produced an empty index set.")
  }
  expression_matrix <- as.matrix(expression_matrix)
  sample_ids <- colnames(expression_matrix)[column_indices]
  surv_idx <- match(sample_ids, survival_df$patient)
  if (anyNA(surv_idx)) {
    missing <- sample_ids[is.na(surv_idx)]
    stop("Missing survival entries for samples: ", paste(missing, collapse = ", "))
  }
  surv_subset <- survival_df[surv_idx, , drop = FALSE]
  sampInfo <- surv_subset
  sampInfo$ID <- sample_ids
  sampInfo$dataset <- dataname
  sampInfo$keep <- 1L
  sampInfo$event <- as.numeric(sampInfo$status)
  rownames(sampInfo) <- sampInfo$ID
  list(
    ex = expression_matrix[, column_indices, drop = FALSE],
    sampInfo = sampInfo,
    featInfo = rownames(expression_matrix),
    samp_keeps = seq_len(length(column_indices)),
    dataname = dataname
  )
}

compute_simulation_bounds <- function(base_bounds = DESURV_BO_BOUNDS) {
  base_bounds %||% list()
}

coalesce_config_values <- function(values, fallback) {
  if (is.null(values) || !length(values)) {
    return(fallback)
  }
  values
}

apply_simulation_hyperparam <- function(cfg,
                                        values,
                                        fallback,
                                        bound_name,
                                        arg_name = bound_name,
                                        config_key = NULL,
                                        default_key = arg_name,
                                        type = c("continuous", "integer"),
                                        log_scale = FALSE) {
  type <- match.arg(type)
  values <- coalesce_config_values(values, fallback)
  if (is.null(values) || !length(values)) {
    stop("No values supplied for hyperparameter `", bound_name, "`.")
  }
  numeric_values <- unique(as.numeric(values))
  if (type == "integer") {
    numeric_values <- unique(as.integer(round(numeric_values)))
  }
  if (!length(numeric_values)) {
    stop("Unable to coerce values for `", bound_name, "` to numeric.")
  }
  default_value <- numeric_values[[1]]
  if (type == "integer") {
    default_value <- as.integer(default_value)
  }
  cfg$hyper_defaults[[default_key]] <- default_value
  tune <- length(numeric_values) > 1
  if (tune) {
    lower <- min(numeric_values)
    upper <- max(numeric_values)
    spec <- list(lower = lower, upper = upper, type = type)
    if (isTRUE(log_scale) && lower > 0) {
      spec$scale <- "log10"
    }
    cfg$bo_bounds[[bound_name]] <- spec
    cfg$bo_extra_args[[arg_name]] <- NULL
  } else {
    cfg$bo_bounds[[bound_name]] <- NULL
    cfg$bo_extra_args[[arg_name]] <- default_value
  }
  if (!is.null(config_key)) {
    cfg[[config_key]] <- numeric_values
  }
  cfg
}

finalize_simulation_hyperparams <- function(cfg) {
  cfg$bo_bounds <- cfg$bo_bounds %||% DESURV_BO_BOUNDS
  cfg$bo_extra_args <- cfg$bo_extra_args %||% list()
  cfg$hyper_defaults <- cfg$hyper_defaults %||% list()

  cfg <- apply_simulation_hyperparam(
    cfg,
    values = cfg$ngene_config,
    fallback = NGENE_CONFIG,
    bound_name = "ngene",
    arg_name = "ngene",
    config_key = "ngene_config",
    default_key = "ngene",
    type = "integer"
  )
  cfg <- apply_simulation_hyperparam(
    cfg,
    values = cfg$ntop_config,
    fallback = NTOP_CONFIG,
    bound_name = "ntop",
    arg_name = "ntop",
    config_key = "ntop_config",
    default_key = "ntop",
    type = "integer"
  )
  cfg <- apply_simulation_hyperparam(
    cfg,
    values = cfg$lambdaW_config,
    fallback = LAMBDAW_CONFIG,
    bound_name = "lambdaW_grid",
    arg_name = "lambdaW_grid",
    config_key = "lambdaW_config",
    default_key = "lambdaW",
    log_scale = TRUE
  )
  cfg <- apply_simulation_hyperparam(
    cfg,
    values = cfg$lambdaH_config,
    fallback = LAMBDAH_CONFIG,
    bound_name = "lambdaH_grid",
    arg_name = "lambdaH_grid",
    config_key = "lambdaH_config",
    default_key = "lambdaH",
    log_scale = TRUE
  )
  cfg
}

extract_param_value <- function(value, default_value, integer = FALSE) {
  if (is.null(value) || is.na(value)) {
    return(default_value)
  }
  if (integer) {
    return(as.integer(round(value)))
  }
  as.numeric(value)
}

run_single_bo <- function(data,
                          run_prefix,
                          bo_bounds = DESURV_BO_BOUNDS,
                          bo_extra_args = list(),
                          hyper_defaults = list(),
                          coarse_control = BO_COARSE_CONTROL,
                          refine_control = BO_REFINE_CONTROL,
                          max_refinements = BO_MAX_REFINEMENTS,
                          tol_gain = BO_TOL_GAIN,
                          plateau = BO_PLATEAU,
                          top_k = BO_TOP_K,
                          shrink_base = BO_SHRINK_BASE,
                          importance_gain = BO_IMPORTANCE_GAIN,
                          parallel_grid = TRUE) {
  bounds <- compute_simulation_bounds(bo_bounds)
  defaults <- list(
    ngene = hyper_defaults$ngene %||% NGENE_DEFAULT,
    ntop = hyper_defaults$ntop %||% NTOP_DEFAULT,
    lambdaW = hyper_defaults$lambdaW %||% LAMBDAW_DEFAULT,
    lambdaH = hyper_defaults$lambdaH %||% LAMBDAH_DEFAULT
  )
  bo_fixed <- c(list(n_starts = NINIT), bo_extra_args %||% list())
  if (length(bo_fixed)) {
    keep <- !vapply(bo_fixed, is.null, logical(1))
    bo_fixed <- bo_fixed[keep]
  }
  parallel_grid <- if (is.null(parallel_grid)) TRUE else isTRUE(parallel_grid)
  bo_results <- DeSurv::desurv_cv_bayesopt_refine(
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
    coarse_bounds = bounds,
    bo_fixed = bo_fixed,
    max_refinements = max_refinements,
    tol_gain = tol_gain,
    plateau = plateau,
    top_k = top_k,
    shrink_base = shrink_base,
    importance_gain = importance_gain,
    coarse_control = coarse_control,
    refine_control = refine_control,
    verbose = TRUE,
    parallel_grid = parallel_grid,
    ncores_grid = NINIT
  )
  params_best <- standardize_bo_params(bo_results$overall_best$params)
  training_dir <- results_root_dir(
    ngene = extract_param_value(params_best$ngene, defaults$ngene, integer = TRUE),
    tol = TOL,
    maxit = MAXIT,
    pkg_version = PKG_VERSION,
    git_branch = GIT_BRANCH,
    train_prefix = run_prefix,
    method_trans_train = METHOD_TRANS_TRAIN
  )
  dir.create(training_dir, recursive = TRUE, showWarnings = FALSE)
  history_path <- file.path(training_dir, "desurv_bo_history.csv")
  utils::write.csv(bo_results$history, history_path, row.names = FALSE)
  list(
    params_best = params_best,
    history_path = history_path,
    training_dir = training_dir,
    best_score = bo_results$overall_best$mean_cindex
  )
}

run_simulation_config <- function(cfg) {
  parallel_grid <- cfg$parallel_grid
  if (is.null(parallel_grid)) {
    parallel_grid <- TRUE
  } else {
    parallel_grid <- isTRUE(parallel_grid)
  }
  replicate_specs <- expand_simulation_replicates(cfg)
  # Preserve legacy parallel flags on each replicate spec so other callers
  # (e.g., targets pipelines) can inspect them, but still route the actual
  # toggle through the explicit `parallel_grid` argument below.
  replicate_specs <- lapply(
    replicate_specs,
    function(spec) {
      spec$parallel_grid <- parallel_grid
      spec
    }
  )
  rep_results <- purrr::map(
    replicate_specs,
    run_simulation_replicate,
    parallel_grid = parallel_grid
  )
  list(
    config_name = cfg$config_name,
    config_signature = cfg$config_signature,
    results = rep_results
  )
}

expand_simulation_replicates <- function(cfg) {
  n_reps <- cfg$n_reps %||% 1L
  base_seed <- cfg$base_seed %||% 1L
  split_seed <- cfg$split_seed %||% (base_seed + 1L)
  parallel_flag <- cfg$parallel_grid
  if (is.null(parallel_flag)) {
    parallel_flag <- TRUE
  } else {
    parallel_flag <- isTRUE(parallel_flag)
  }
  dir.create(SIMULATION_RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)
  purrr::map(
    seq_len(n_reps),
    function(rep_id) {
      list(
        config_name = cfg$config_name,
        config_signature = cfg$config_signature,
        simulator = cfg$simulator,
        args = cfg$args %||% list(),
        design_name = cfg$design_name,
        design_params = cfg$design_params %||% list(),
        expression_transform = cfg$expression_transform %||% "log2",
        train_fraction = cfg$train_fraction %||% 0.7,
        run_prefix = cfg$run_prefix,
        replicate = rep_id,
        run_id = paste0(cfg$run_prefix, "_rep", sprintf("%03d", rep_id)),
        sim_seed = base_seed + rep_id - 1L,
        split_seed = split_seed + rep_id - 1L,
        bo_bounds = cfg$bo_bounds %||% DESURV_BO_BOUNDS,
        bo_extra_args = cfg$bo_extra_args %||% list(),
        hyper_defaults = cfg$hyper_defaults %||% list(),
        bo_coarse_control = cfg$bo_coarse_control %||% BO_COARSE_CONTROL,
        bo_refine_control = cfg$bo_refine_control %||% BO_REFINE_CONTROL,
        bo_max_refinements = cfg$bo_max_refinements %||% BO_MAX_REFINEMENTS,
        bo_tol_gain = cfg$bo_tol_gain %||% BO_TOL_GAIN,
        bo_plateau = cfg$bo_plateau %||% BO_PLATEAU,
        bo_top_k = cfg$bo_top_k %||% BO_TOP_K,
        bo_shrink_base = cfg$bo_shrink_base %||% BO_SHRINK_BASE,
        bo_importance_gain = cfg$bo_importance_gain %||% BO_IMPORTANCE_GAIN,
        parallel_grid = parallel_flag
      )
    }
  )
}

expand_all_simulation_replicates <- function(config_list) {
  config_list <- config_list %||% list()
  replicate_lists <- purrr::map(config_list, expand_simulation_replicates)
  purrr::flatten(replicate_lists)
}

run_simulation_replicate <- function(rep_cfg, parallel_grid = NULL) {
  if (is.null(parallel_grid)) {
    parallel_grid <- rep_cfg$parallel_grid
  }
  if (is.null(parallel_grid)) {
    parallel_grid <- TRUE
  } else {
    parallel_grid <- isTRUE(parallel_grid)
  }
  simulator <- get(rep_cfg$simulator, mode = "function")
  set.seed(rep_cfg$sim_seed)
  sim <- do.call(simulator, rep_cfg$args)
  expression_matrix <- transform_simulation_expression(
    sim$counts,
    rep_cfg$expression_transform
  )
  survival_df <- sim$surv
  training_dataset <- build_simulation_dataset(
    expression_matrix = expression_matrix,
    survival_df = survival_df,
    column_indices = seq_len(ncol(expression_matrix)),
    dataname = rep_cfg$run_prefix
  )
  bo_summary <- run_single_bo(
    training_dataset,
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
    parallel_grid = parallel_grid
  )
  list(
    config_name = rep_cfg$config_name,
    config_signature = rep_cfg$config_signature,
    replicate = rep_cfg$replicate,
    run_id = rep_cfg$run_id,
    simulator = rep_cfg$simulator,
    args = rep_cfg$args,
    design_name = rep_cfg$design_name,
    design_params = rep_cfg$design_params,
    expression_transform = rep_cfg$expression_transform,
    train_fraction = rep_cfg$train_fraction,
    bo_summary = bo_summary
  )
}

assemble_simulation_config_results <- function(config_spec, replicate_results) {
  replicate_results <- replicate_results %||% list()
  config_reps <- purrr::keep(
    replicate_results,
    ~ identical(.x$config_signature, config_spec$config_signature)
  )
  expected_n <- config_spec$n_reps %||% 1L
  if (length(config_reps) != expected_n) {
    warning(
      "Config `", config_spec$config_name, "` expected ", expected_n,
      " replicates but found ", length(config_reps), "."
    )
  }
  list(
    config_name = config_spec$config_name,
    config_signature = config_spec$config_signature,
    design_name = config_spec$design_name,
    design_params = config_spec$design_params,
    results = config_reps
  )
}

build_simulation_config_list <- function(scenarios) {
  scenarios <- scenarios %||% list()
  config_lists <- purrr::imap(
    scenarios,
    function(cfg, name) {
      design_specs <- extract_simulation_designs(cfg, name)
      purrr::map(
        design_specs,
        function(design_spec) {
          cfg_copy <- rlang::duplicate(cfg, shallow = FALSE)
          cfg_copy$config_name <- name
          cfg_copy$design_name <- design_spec$label
          cfg_copy$design_params <- design_spec$params
          cfg_copy$args <- design_spec$args
          args_label <- format_simulation_args(cfg_copy$args %||% list())
          signature_parts <- c(
            name,
            cfg_copy$simulator,
            args_label,
            cfg_copy$expression_transform %||% "log2",
            cfg_copy$train_fraction %||% 0.7
          )
          if (isTRUE(design_spec$explicit) &&
            !is.null(cfg_copy$design_name) &&
            nzchar(cfg_copy$design_name)) {
            signature_parts <- c(signature_parts, cfg_copy$design_name)
          }
          signature <- paste(signature_parts, collapse = "__")
          cfg_copy$config_signature <- signature
          cfg_copy$run_prefix <- paste0("SIM_", gsub("[^A-Za-z0-9]+", "_", signature))
          cfg_copy <- finalize_simulation_hyperparams(cfg_copy)
          cfg_copy
        }
      )
    }
  )
  purrr::flatten(config_lists)
}

summarize_simulation_config <- function(config_result) {
  config_name <- config_result$config_name
  config_signature <- config_result$config_signature
  config_design_name <- config_result$design_name
  config_design_params <- config_result$design_params %||% list()
  replicate_results <- config_result$results %||% list()
  if (!length(replicate_results)) {
    return(
      tibble::tibble(
        config_name = config_name,
        config_signature = config_signature,
        design_name = config_design_name %||% NA_character_,
        design_args = format_simulation_args(config_design_params),
        design_params = list(config_design_params)
      )[0, ]
    )
  }
  purrr::map_dfr(
    replicate_results,
    function(rep_result) {
      summary <- rep_result$bo_summary %||% list()
      params_best <- summary$params_best %||% list()
      param_cols <- if (length(params_best)) {
        tibble::as_tibble_row(params_best)
      } else {
        tibble::tibble()
      }
      rep_design_params <- rep_result$design_params %||% config_design_params
      rep_design_name <- rep_result$design_name %||% config_design_name
      result_row <- tibble::tibble(
        config_name = config_name,
        config_signature = config_signature,
        design_name = rep_design_name %||% NA_character_,
        design_args = format_simulation_args(rep_design_params %||% list()),
        design_params = list(rep_design_params %||% list()),
        simulator = rep_result$simulator %||% NA_character_,
        expression_transform = rep_result$expression_transform %||% NA_character_,
        train_fraction = rep_result$train_fraction %||% NA_real_,
        replicate = rep_result$replicate %||% NA_integer_,
        run_id = rep_result$run_id %||% NA_character_,
        best_score = summary$best_score %||% NA_real_,
        training_dir = summary$training_dir %||% NA_character_,
        history_path = summary$history_path %||% NA_character_
      )
      if (ncol(param_cols)) {
        result_row[ names(param_cols) ] <- as.list(param_cols)
      }
      result_row
    }
  )
}
