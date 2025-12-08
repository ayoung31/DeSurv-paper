library(targets)
library(tarchetypes)
library(crew)
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(DeSurv)
})

sim_files <- list.files(
  "R/simulation_functions",
  pattern = "[.]R$",
  full.names = TRUE
)
purrr::walk(sim_files, source)

tar_option_set(
  packages = c("dplyr", "purrr", "tibble", "DeSurv", "survival"),
  format = "rds",
  controller = crew::crew_controller_sequential(),
  error = "continue"
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

standardize_bayes_params <- function(params) {
  if (is.null(params)) {
    return(list())
  }
  out <- as.list(params)
  if (!length(out)) {
    return(out)
  }
  names(out) <- sub("_grid$", "", names(out))
  out
}

SIM_DATASETS_PER_SCENARIO <- 100L
SIM_GLOBAL_SEED <- 101L
SIM_METHOD_TRANSFORM <- "none"
SIM_DESURV_TOL <- 1e-5
SIM_DESURV_MAXIT <- 2000L
SIM_FINAL_NINIT <- 1L
SIM_DEFAULT_NGENE <- NULL
SIM_DEFAULT_NTOP <- NULL
SIM_DEFAULT_LAMBDAW <- 0
SIM_DEFAULT_LAMBDAH <- 0
SIM_DEFAULT_K <- 3L
SIM_DEFAULT_ALPHA <- 0.6
SIM_DEFAULT_LAMBDA <- 0.1
SIM_DEFAULT_NU <- 0.3
SIM_CV_NFOLDS <- 5L
SIM_CV_NSTARTS <- 1L
SIM_TRAIN_FRACTION <- 0.7
SIM_SPLIT_SEED_OFFSET <- 10000L

SIM_DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 6L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e2, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

SIM_BO_N_INIT <- 8L
SIM_BO_N_ITER <- 10L
SIM_BO_CANDIDATE_POOL <- 1000L
SIM_BO_EXPLORATION_WEIGHT <- 0.01

SIMULATION_SCENARIOS <- list(
  list(
    scenario_id = "R0_easy",
    scenario = "R0",
    description = "Default easy/sanity scenario",
    replicates = SIM_DATASETS_PER_SCENARIO,
    seed_offset = SIM_GLOBAL_SEED,
    overrides = list()
  )#,
  # list(
  #   scenario_id = "R2_correlated",
  #   scenario = "R2",
  #   description = "Correlated programs with moderate noise",
  #   replicates = SIM_DATASETS_PER_SCENARIO,
  #   seed_offset = SIM_GLOBAL_SEED + 1000L,
  #   overrides = list()
  # ),
  # list(
  #   scenario_id = "R3_overlap",
  #   scenario = "R3",
  #   description = "Overlapping markers and strong background",
  #   replicates = SIM_DATASETS_PER_SCENARIO,
  #   seed_offset = SIM_GLOBAL_SEED + 2000L,
  #   overrides = list()
  # )
)

SIM_FIXED_PARAMS <- list(
  k = NULL,
  alpha = SIM_DEFAULT_ALPHA,
  lambda = SIM_DEFAULT_LAMBDA,
  nu = SIM_DEFAULT_NU,
  lambdaW = SIM_DEFAULT_LAMBDAW,
  lambdaH = SIM_DEFAULT_LAMBDAH,
  ngene = SIM_DEFAULT_NGENE
)

SIM_ANALYSIS_SPECS <- list(
  list(
    analysis_id = "fixed",
    label = "Fixed parameters",
    mode = "fixed",
    params = SIM_FIXED_PARAMS
  ),
  list(
    analysis_id = "fixed_alpha0",
    label = "Fixed with alpha = 0",
    mode = "fixed",
    params = modifyList(SIM_FIXED_PARAMS, list(alpha = 0))
  ),
  list(
    analysis_id = "bo",
    label = "Bayesian optimization",
    mode = "bayesopt",
    bounds = SIM_DESURV_BO_BOUNDS,
    bo_fixed = list(
      ngene = SIM_DEFAULT_NGENE,
      ntop = SIM_DEFAULT_NTOP,
      lambdaW_grid = SIM_DEFAULT_LAMBDAW,
      lambdaH_grid = SIM_DEFAULT_LAMBDAH
    )
  )#,
  # list(
  #   analysis_id = "bo_alpha0",
  #   label = "Bayesian optimization with alpha = 0",
  #   mode = "bayesopt",
  #   bounds = modifyList(
  #     SIM_DESURV_BO_BOUNDS,
  #     list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"))
  #   ),
  #   bo_fixed = list(
  #     ngene = SIM_DEFAULT_NGENE,
  #     ntop = SIM_DEFAULT_NTOP,
  #     alpha_grid = 0,
  #     lambdaW_grid = SIM_DEFAULT_LAMBDAW,
  #     lambdaH_grid = SIM_DEFAULT_LAMBDAH
  #   ),
  #   final_overrides = list(alpha = 0)
  # )
)

build_simulation_dataset_specs <- function(
    scenarios,
    replicates_per_scenario = SIM_DATASETS_PER_SCENARIO,
    base_seed = SIM_GLOBAL_SEED) {
  purrr::imap(scenarios, function(scenario, idx) {
    scenario_id <- scenario$scenario_id %||%
      names(scenarios)[[idx]] %||%
      sprintf("scenario_%02d", idx)
    scenario_name <- scenario$scenario %||% scenario_id
    replicates <- scenario$replicates %||% replicates_per_scenario
    seed_offset <- scenario$seed_offset %||% (base_seed + (idx - 1L) * 1000L)
    overrides <- scenario$overrides %||% list()
    description <- scenario$description %||% ""
    purrr::map(seq_len(replicates), function(rep_id) {
      list(
        scenario_id = scenario_id,
        scenario_name = scenario_name,
        replicate = rep_id,
        description = description,
        seed = as.integer(seed_offset + rep_id - 1L),
        overrides = overrides
      )
    })
  }) |> purrr::flatten()
}

generate_simulation_dataset <- function(spec) {
  args <- c(
    list(
      scenario = spec$scenario_name,
      seed = spec$seed
    ),
    spec$overrides %||% list()
  )
  sim <- do.call(simulate_desurv_scenario, args)
  list(
    spec = spec,
    simulation = sim,
    data = format_simulation_data(sim, spec),
    metadata = list(
      scenario_id = spec$scenario_id,
      scenario_name = spec$scenario_name,
      replicate = spec$replicate,
      seed = spec$seed,
      description = spec$description
    )
  )
}

format_simulation_data <- function(sim, spec) {
  expr <- sim$X
  if (is.null(rownames(expr))) {
    rownames(expr) <- sprintf("gene_%05d", seq_len(nrow(expr)))
  }
  if (is.null(colnames(expr))) {
    colnames(expr) <- sprintf(
      "%s_rep%03d_s%03d",
      spec$scenario_id,
      spec$replicate,
      seq_len(ncol(expr))
    )
  }
  n_samples <- ncol(expr)
  ids <- colnames(expr)
  time_vec <- sim$time
  status_vec <- sim$status
  if (length(time_vec) != n_samples) {
    stop(
      sprintf(
        "Simulated survival time length (%d) does not match sample count (%d).",
        length(time_vec),
        n_samples
      ),
      call. = FALSE
    )
  }
  if (length(status_vec) != n_samples) {
    stop(
      sprintf(
        "Simulated survival status length (%d) does not match sample count (%d).",
        length(status_vec),
        n_samples
      ),
      call. = FALSE
    )
  }
  samp_info <- tibble::tibble(
    ID = ids,
    dataset = spec$scenario_id,
    scenario = spec$scenario_name,
    replicate = spec$replicate,
    time = time_vec,
    event = status_vec
  )
  samp_info <- as.data.frame(samp_info, stringsAsFactors = FALSE)
  rownames(samp_info) <- samp_info$ID
  list(
    ex = expr,
    sampInfo = samp_info,
    samp_keeps = rep(TRUE, ncol(expr)),
    dataname = paste0(spec$scenario_id, "_rep", sprintf("%03d", spec$replicate))
  )
}

split_simulation_samples <- function(dataset_entry,
                                     train_fraction = SIM_TRAIN_FRACTION,
                                     seed_offset = SIM_SPLIT_SEED_OFFSET) {
  ids <- colnames(dataset_entry$data$ex)
  if (length(ids) < 2) {
    stop("Need at least two samples to perform a train/test split.", call. = FALSE)
  }
  train_fraction <- min(max(train_fraction, 0), 1)
  n_train <- round(length(ids) * train_fraction)
  n_train <- min(max(n_train, 1L), length(ids) - 1L)
  base_seed <- dataset_entry$metadata$seed %||% 1L
  split_seed <- as.integer(base_seed + seed_offset)
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- NULL
  if (seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(split_seed)
  train_ids <- sample(ids, size = n_train, replace = FALSE)
  test_ids <- ids[!(ids %in% train_ids)]
  if (!length(test_ids)) {
    stop("Failed to allocate samples to the test split.", call. = FALSE)
  }
  list(train_ids = train_ids, test_ids = test_ids)
}

prepare_simulation_data <- function(dataset_entry,
                                    sample_ids = NULL,
                                    ngene = SIM_DEFAULT_NGENE,
                                    transform_method = SIM_METHOD_TRANSFORM,
                                    genes = NULL) {
  X <- dataset_entry$data$ex
  samp_info <- dataset_entry$data$sampInfo
  transform_method <- match.arg(transform_method, c("rank", "none"))
  if (identical(transform_method, "rank")) {
    ranks <- apply(X, 2, rank, ties.method = "average")
    dim(ranks) <- dim(X)
    dimnames(ranks) <- dimnames(X)
    X <- ranks
  }
  if (!is.null(sample_ids)) {
    missing_ids <- setdiff(sample_ids, colnames(X))
    if (length(missing_ids)) {
      stop("Requested sample IDs are not present in the simulation data.", call. = FALSE)
    }
    keep <- match(sample_ids, colnames(X))
    X <- X[, keep, drop = FALSE]
    samp_rows <- match(sample_ids, samp_info$ID)
    if (any(is.na(samp_rows))) {
      stop("Sample metadata missing for requested IDs.", call. = FALSE)
    }
    samp_info <- samp_info[samp_rows, , drop = FALSE]
  } else {
    ordered_ids <- colnames(X)
    samp_rows <- match(ordered_ids, samp_info$ID)
    if (any(is.na(samp_rows))) {
      stop("Sample metadata missing for some simulation IDs.", call. = FALSE)
    }
    samp_info <- samp_info[samp_rows, , drop = FALSE]
  }
  if (!is.null(genes)) {
    keep_genes <- genes[genes %in% rownames(X)]
    if (!length(keep_genes)) {
      stop("Requested genes are not present in the simulation data.", call. = FALSE)
    }
    X <- X[keep_genes, , drop = FALSE]
  } else if (!is.null(ngene) && length(ngene) > 0 && !is.na(ngene)) {
    ng <- as.integer(min(ngene, nrow(X)))
    vars <- apply(X, 1, stats::var)
    vars[is.na(vars)] <- 0
    order_idx <- order(vars, decreasing = TRUE)
    keep_idx <- head(order_idx, ng)
    X <- X[keep_idx, , drop = FALSE]
  }
  list(
    ex = X,
    sampInfo = samp_info,
    samp_keeps = rep(TRUE, ncol(X)),
    dataname = dataset_entry$data$dataname
  )
}

predict_dataset_risk <- function(fit, dataset) {
  stopifnot(!is.null(dataset$ex), !is.null(fit$W), !is.null(fit$beta))
  genes <- intersect(rownames(fit$W), rownames(dataset$ex))
  if (!length(genes)) {
    warning("No overlapping genes between fit and dataset: ", dataset$dataname %||% "split")
    return(tibble::tibble())
  }
  W_sub <- fit$W[genes, , drop = FALSE]
  X_sub <- dataset$ex[genes, , drop = FALSE]
  if (nrow(W_sub) == 0L || ncol(W_sub) == 0L) {
    warning("Model does not contain usable factors for prediction.")
    return(tibble::tibble())
  }
  Z <- t(X_sub) %*% W_sub
  if (ncol(Z) == 0L) {
    warning("No latent representation available for prediction.")
    return(tibble::tibble())
  }
  beta <- as.numeric(fit$beta)
  if (length(beta) != ncol(Z)) {
    stop("Mismatch between beta length and latent space dimension.", call. = FALSE)
  }
  risk <- drop(Z %*% beta)
  tibble::tibble(
    sample_id = colnames(dataset$ex),
    risk_score = -risk
  )
}

compute_dataset_cindex <- function(fit, dataset) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required to compute concordance.", call. = FALSE)
  }
  preds <- predict_dataset_risk(fit, dataset)
  if (nrow(preds) == 0L) {
    return(NA_real_)
  }
  samp <- dataset$sampInfo
  rownames(samp) <- samp$ID
  idx <- match(preds$sample_id, rownames(samp))
  if (any(is.na(idx))) {
    stop("Missing survival metadata for predicted samples.", call. = FALSE)
  }
  surv_obj <- survival::Surv(samp$time[idx], samp$event[idx])
  cc <- tryCatch(
    survival::concordance(surv_obj ~ preds$risk_score),
    error = function(e) survival::survConcordance(surv_obj ~ preds$risk_score)
  )
  as.numeric(cc$concordance)
}

coerce_int <- function(value, default) {
  has_value <- !is.null(value) && length(value) > 0 && !is.na(value)
  if (!has_value) {
    if (is.null(default) || length(default) == 0 || is.na(default)) {
      return(NULL)
    }
    return(as.integer(round(default)))
  }
  as.integer(round(value))
}

coerce_num <- function(value, default) {
  has_value <- !is.null(value) && length(value) > 0 && !is.na(value)
  if (!has_value) {
    if (is.null(default) || length(default) == 0 || is.na(default)) {
      return(NULL)
    }
    return(as.numeric(default))
  }
  as.numeric(value)
}

get_simulation_k <- function(dataset_entry) {
  sim <- dataset_entry$simulation %||% list()
  params <- sim$params %||% list()
  k_value <- params$K %||% params$k
  if (is.null(k_value) && !is.null(sim$W)) {
    k_value <- ncol(sim$W)
  }
  if (is.null(k_value) && !is.null(sim$H)) {
    k_value <- ncol(sim$H)
  }
  if (is.null(k_value)) {
    k_value <- SIM_DEFAULT_K
  }
  as.integer(k_value)
}

build_result_row <- function(dataset_entry,
                             analysis_spec,
                             processed,
                             fit,
                             params,
                             bo_details = NULL,
                             split = NULL,
                             train_cindex = NULL,
                             test_cindex = NULL) {
  split_info <- split %||% list()
  train_val <- if (is.null(train_cindex)) fit$cindex %||% NA_real_ else train_cindex
  test_val <- if (is.null(test_cindex)) fit$cindex %||% NA_real_ else test_cindex
  n_train <- if (!is.null(split_info$train_ids)) length(split_info$train_ids) else NA_integer_
  n_test <- if (!is.null(split_info$test_ids)) length(split_info$test_ids) else NA_integer_
  tibble::tibble(
    scenario_id = dataset_entry$metadata$scenario_id,
    scenario = dataset_entry$metadata$scenario_name,
    replicate = dataset_entry$metadata$replicate,
    seed = dataset_entry$metadata$seed,
    dataset_name = dataset_entry$data$dataname,
    analysis_id = analysis_spec$analysis_id,
    analysis_mode = analysis_spec$mode,
    n_train = n_train,
    n_test = n_test,
    train_cindex = train_val,
    cindex = test_val,
    params = list(params),
    bo_history = list(bo_details$history %||% NULL),
    bo_summary = list(bo_details$summary %||% NULL),
    fit = list(fit),
    processed = list(processed),
    truth = list(dataset_entry$simulation),
    split = list(split_info)
  )
}

format_param_value <- function(value) {
  if (is.null(value) || !length(value)) {
    return(NA_character_)
  }
  if (is.numeric(value)) {
    formatted <- format(value, digits = 6, trim = TRUE, scientific = FALSE)
  } else {
    formatted <- as.character(value)
  }
  formatted <- formatted[!is.na(formatted)]
  if (!length(formatted)) {
    return(NA_character_)
  }
  paste(formatted, collapse = ";")
}

format_result_params <- function(params) {
  if (is.null(params) || !length(params)) {
    return(NA_character_)
  }
  keep <- !vapply(
    params,
    function(x) is.null(x) || !length(x) || all(is.na(x)),
    logical(1)
  )
  if (!any(keep)) {
    return(NA_character_)
  }
  params <- params[keep]
  pieces <- purrr::imap_chr(params, function(value, name) {
    value_str <- format_param_value(value)
    if (is.na(value_str)) {
      sprintf("%s=NA", name)
    } else {
      sprintf("%s=%s", name, value_str)
    }
  })
  paste(pieces, collapse = ", ")
}

summarize_simulation_results <- function(result_list) {
  empty_tbl <- tibble::tibble(
    scenario = character(),
    model_run = character(),
    parameters = character(),
    cindex = numeric()
  )
  result_list <- purrr::compact(result_list)
  if (!length(result_list)) {
    return(empty_tbl)
  }
  required_cols <- c(
    "scenario_id",
    "scenario",
    "replicate",
    "dataset_name",
    "analysis_id",
    "params",
    "cindex"
  )
  is_valid <- purrr::map_lgl(
    result_list,
    ~ inherits(.x, "data.frame") && all(required_cols %in% names(.x))
  )
  if (!any(is_valid)) {
    stop(
      "No valid simulation analysis results were produced. ",
      "Inspect sim_analysis_result for errors.",
      call. = FALSE
    )
  }
  invalid_count <- sum(!is_valid)
  if (invalid_count > 0) {
    warning(
      sprintf(
        "Dropping %d simulation result(s) missing required metadata.",
        invalid_count
      ),
      call. = FALSE
    )
  }
  result_list <- result_list[is_valid]
  if (!length(result_list)) {
    return(empty_tbl)
  }
  result_tbl <- dplyr::bind_rows(result_list)
  if (!nrow(result_tbl)) {
    return(empty_tbl)
  }
  dataset_name <- result_tbl$dataset_name
  dataset_valid <- !is.na(dataset_name) & nzchar(dataset_name)
  scenario_id <- result_tbl$scenario_id
  scenario_valid <- !is.na(scenario_id) & nzchar(scenario_id)
  replicate_vals <- result_tbl$replicate
  replicate_label <- ifelse(
    !is.na(replicate_vals),
    sprintf("rep%03d", as.integer(replicate_vals)),
    NA_character_
  )
  scenario_fallback <- scenario_id
  scenario_fallback[scenario_valid & !is.na(replicate_label)] <- paste0(
    scenario_fallback[scenario_valid & !is.na(replicate_label)],
    "_",
    replicate_label[scenario_valid & !is.na(replicate_label)]
  )
  scenario_fallback[!scenario_valid] <- NA_character_
  default_ids <- sprintf("run_%03d", seq_len(nrow(result_tbl)))
  base_run <- ifelse(dataset_valid, dataset_name, scenario_fallback)
  base_missing <- is.na(base_run) | !nzchar(base_run)
  base_run[base_missing] <- default_ids[base_missing]
  analysis_id <- result_tbl$analysis_id
  analysis_valid <- !is.na(analysis_id) & nzchar(analysis_id)
  model_run <- base_run
  model_run[analysis_valid] <- paste(model_run[analysis_valid], analysis_id[analysis_valid], sep = "::")
  tibble::tibble(
    scenario = result_tbl$scenario,
    model_run = model_run,
    parameters = purrr::map_chr(result_tbl$params, format_result_params),
    cindex = result_tbl$cindex
  )
}

run_fixed_analysis <- function(dataset_entry, analysis_spec) {
  params <- analysis_spec$params %||% SIM_FIXED_PARAMS
  k_default <- get_simulation_k(dataset_entry)
  k_value <- coerce_int(params$k, k_default)
  params$k <- k_value
  ngene_value <- coerce_int(params$ngene, SIM_DEFAULT_NGENE)
  split_ids <- split_simulation_samples(dataset_entry)
  processed_train <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = ngene_value,
    transform_method = SIM_METHOD_TRANSFORM
  )
  processed_test <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$test_ids,
    transform_method = SIM_METHOD_TRANSFORM,
    genes = rownames(processed_train$ex)
  )
  fit <- DeSurv::desurv_fit(
    X = processed_train$ex,
    y = processed_train$sampInfo$time,
    d = processed_train$sampInfo$event,
    k = k_value,
    alpha = coerce_num(params$alpha, SIM_DEFAULT_ALPHA),
    lambda = coerce_num(params$lambda, SIM_DEFAULT_LAMBDA),
    nu = coerce_num(params$nu, SIM_DEFAULT_NU),
    lambdaW = coerce_num(params$lambdaW, SIM_DEFAULT_LAMBDAW),
    lambdaH = coerce_num(params$lambdaH, SIM_DEFAULT_LAMBDAH),
    seed = dataset_entry$metadata$seed,
    tol = SIM_DESURV_TOL,
    tol_init = SIM_DESURV_TOL,
    maxit = SIM_DESURV_MAXIT,
    imaxit = SIM_DESURV_MAXIT,
    ninit = SIM_FINAL_NINIT,
    parallel_init = FALSE,
    verbose = TRUE
  )
  fit$data <- list(X = processed_train$ex, sampInfo = processed_train$sampInfo)
  train_cindex <- compute_dataset_cindex(fit, processed_train)
  test_cindex <- compute_dataset_cindex(fit, processed_test)
  build_result_row(
    dataset_entry = dataset_entry,
    analysis_spec = analysis_spec,
    processed = list(train = processed_train, test = processed_test),
    fit = fit,
    params = params,
    split = split_ids,
    train_cindex = train_cindex,
    test_cindex = test_cindex
  )
}

run_bayesopt_analysis <- function(dataset_entry, analysis_spec) {
  bounds <- analysis_spec$bounds %||% SIM_DESURV_BO_BOUNDS
  bo_fixed <- modifyList(
    list(
      n_starts = SIM_CV_NSTARTS,
      nfolds = SIM_CV_NFOLDS,
      tol = SIM_DESURV_TOL,
      maxit = SIM_DESURV_MAXIT
    ),
    analysis_spec$bo_fixed %||% list()
  )
  n_init <- analysis_spec$n_init %||% SIM_BO_N_INIT
  n_iter <- analysis_spec$n_iter %||% SIM_BO_N_ITER
  candidate_pool <- analysis_spec$candidate_pool %||% SIM_BO_CANDIDATE_POOL
  exploration_weight <- analysis_spec$exploration_weight %||% SIM_BO_EXPLORATION_WEIGHT
  split_ids <- split_simulation_samples(dataset_entry)
  bo_data <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = bo_fixed$ngene %||% SIM_DEFAULT_NGENE,
    transform_method = SIM_METHOD_TRANSFORM
  )
  bo_results <- DeSurv::desurv_bayesopt(
    X = bo_data$ex,
    y = bo_data$sampInfo$time,
    d = bo_data$sampInfo$event,
    dataset = bo_data$sampInfo$dataset,
    samp_keeps = bo_data$samp_keeps,
    preprocess = FALSE,
    method_trans_train = SIM_METHOD_TRANSFORM,
    engine = "warmstart",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    n_init = n_init,
    n_iter = n_iter,
    candidate_pool = candidate_pool,
    exploration_weight = exploration_weight,
    seed = dataset_entry$metadata$seed,
    cv_verbose = analysis_spec$cv_verbose %||% FALSE,
    verbose = TRUE
  )
  params_best <- standardize_bayes_params(bo_results$best$params)
  final_params <- modifyList(params_best, analysis_spec$final_overrides %||% list())
  dataset_k <- get_simulation_k(dataset_entry)
  final_params$ngene <- coerce_int(final_params$ngene, bo_fixed$ngene %||% SIM_DEFAULT_NGENE)
  final_params$k <- coerce_int(final_params$k, dataset_k)
  final_params$alpha <- coerce_num(final_params$alpha, SIM_DEFAULT_ALPHA)
  final_params$lambda <- coerce_num(final_params$lambda, SIM_DEFAULT_LAMBDA)
  final_params$nu <- coerce_num(final_params$nu, SIM_DEFAULT_NU)
  final_params$lambdaW <- coerce_num(final_params$lambdaW, SIM_DEFAULT_LAMBDAW)
  final_params$lambdaH <- coerce_num(final_params$lambdaH, SIM_DEFAULT_LAMBDAH)
  processed_train <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = final_params$ngene,
    transform_method = SIM_METHOD_TRANSFORM
  )
  processed_test <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$test_ids,
    transform_method = SIM_METHOD_TRANSFORM,
    genes = rownames(processed_train$ex)
  )
  fit <- DeSurv::desurv_fit(
    X = processed_train$ex,
    y = processed_train$sampInfo$time,
    d = processed_train$sampInfo$event,
    k = final_params$k,
    alpha = final_params$alpha,
    lambda = final_params$lambda,
    nu = final_params$nu,
    lambdaW = final_params$lambdaW,
    lambdaH = final_params$lambdaH,
    seed = dataset_entry$metadata$seed,
    tol = SIM_DESURV_TOL,
    tol_init = SIM_DESURV_TOL,
    maxit = SIM_DESURV_MAXIT,
    imaxit = SIM_DESURV_MAXIT,
    ninit = SIM_FINAL_NINIT,
    parallel_init = FALSE,
    verbose = FALSE
  )
  fit$data <- list(X = processed_train$ex, sampInfo = processed_train$sampInfo)
  train_cindex <- compute_dataset_cindex(fit, processed_train)
  test_cindex <- compute_dataset_cindex(fit, processed_test)
  bo_details <- list(
    history = bo_results$history,
    summary = list(best_score = bo_results$best$mean_cindex)
  )
  build_result_row(
    dataset_entry = dataset_entry,
    analysis_spec = analysis_spec,
    processed = list(train = processed_train, test = processed_test),
    fit = fit,
    params = final_params,
    bo_details = bo_details,
    split = split_ids,
    train_cindex = train_cindex,
    test_cindex = test_cindex
  )
}

run_simulation_analysis <- function(dataset_entry, analysis_spec) {
  mode <- analysis_spec$mode %||% "fixed"
  if (identical(mode, "fixed")) {
    run_fixed_analysis(dataset_entry, analysis_spec)
  } else if (identical(mode, "bayesopt")) {
    run_bayesopt_analysis(dataset_entry, analysis_spec)
  } else {
    stop(sprintf("Unsupported analysis mode '%s'.", mode), call. = FALSE)
  }
}

targets_list <- list(
  tar_target(
    sim_dataset_specs,
    build_simulation_dataset_specs(
      scenarios = SIMULATION_SCENARIOS,
      replicates_per_scenario = SIM_DATASETS_PER_SCENARIO,
      base_seed = SIM_GLOBAL_SEED
    )
  ),
  tar_target(
    sim_dataset,
    {
      sim_dataset_specs=sim_dataset_specs[[1]]
      generate_simulation_dataset(sim_dataset_specs)
    },
    pattern = map(sim_dataset_specs),
    iteration = "list"
  ),
  tar_target(
    sim_analysis_spec,
    SIM_ANALYSIS_SPECS,
    iteration = "list"
  ),
  tar_target(
    sim_analysis_result,
    {
      run_simulation_analysis(sim_dataset, sim_analysis_spec)
    },
    pattern = cross(sim_dataset, sim_analysis_spec),
    iteration = "list"
  ),
  tar_target(
    sim_results_table,
    {
      summarize_simulation_results(sim_analysis_result)
    }
    
  )
)

targets_list
