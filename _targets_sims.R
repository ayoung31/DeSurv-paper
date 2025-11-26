NINIT <- 50
NINIT_FULL <- 100
BO_N_INIT <- 20
BO_N_ITER <- 100
BO_CANDIDATE_POOL <- 2000
BO_MAX_REFINEMENTS <- 2
BO_TOL_GAIN <- 0.002
BO_PLATEAU <- 1
BO_TOP_K <- 10
BO_SHRINK_BASE <- 0.5
BO_IMPORTANCE_GAIN <- 0.3
BO_COARSE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 123,
  cv_verbose = FALSE
)
BO_REFINE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 456,
  cv_verbose = FALSE
)

SIMULATION_FUNCTION <- "simulate_desurv_data"
SIMULATION_ARGS <- list(G = 5000, N = 200, K = 4)
SIMULATION_SEED <- 123
SIMULATION_SPLIT_SEED <- 456
SIMULATION_TRAIN_FRACTION <- 0.7
SIMULATION_TRANSFORM <- "log2"
SIMULATION_DATASET_LABEL <- "SIM"
USE_TRAIN_GENES_FOR_VAL <- TRUE

TRAIN_DATASETS <- SIMULATION_DATASET_LABEL
TRAIN_PREFIX <- paste0(TRAIN_DATASETS, collapse = ".")

DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 8L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

NGENE_CONFIG <- c(1000, 5000)
NTOP_CONFIG <- c(25, 200)
LAMBDAW_CONFIG <- c(0)
LAMBDAH_CONFIG <- c(1e-3, 1e3)

source("targets_setup.R")
purrr::walk(
  list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)
source("targets_common_pipeline.R")

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
      rank = {
        ranked <- apply(counts, 2, function(col) rank(col, ties.method = "average"))
        ranked
      },
      none = counts
    )
  } else if (is.function(transform_spec)) {
    expr <- transform_spec(counts)
    expr <- as.matrix(expr)
  } else {
    stop("`transform_spec` must be \"log2\", \"rank\", \"none\", or a function.")
  }

  if (is.null(dim(expr))) {
    expr <- matrix(expr, nrow = nrow(counts), ncol = ncol(counts))
  }
  dimnames(expr) <- dimnames(counts)
  expr
}

build_simulation_dataset <- function(expression_matrix,
                                     survival_df,
                                     column_indices,
                                     dataname) {
  if (length(column_indices) == 0) {
    stop("No samples available for dataset `", dataname, "`.")
  }
  expression_matrix <- as.matrix(expression_matrix)
  survival_df <- as.data.frame(survival_df)

  sample_ids <- colnames(expression_matrix)[column_indices]

  surv_idx <- match(sample_ids, survival_df$patient)
  if (anyNA(surv_idx)) {
    stop("Missing survival rows for samples: ",
         paste(sample_ids[is.na(surv_idx)], collapse = ", "))
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
    samp_keeps = seq_len(ncol(expression_matrix[, column_indices, drop = FALSE])),
    dataname = dataname
  )
}

targets_list <- list(
  tar_target(
    simulation_raw,
    {
      set.seed(SIMULATION_SEED)
      sim_fun <- get(SIMULATION_FUNCTION, mode = "function")
      do.call(sim_fun, SIMULATION_ARGS)
    }
  ),

  tar_target(
    simulation_truth,
    list(
      W = simulation_raw$W_true,
      H = simulation_raw$H_true,
      beta = simulation_raw$beta_true
    )
  ),

  tar_target(
    simulation_expression,
    transform_simulation_expression(simulation_raw$counts, SIMULATION_TRANSFORM)
  ),

  tar_target(
    simulation_survival,
    simulation_raw$surv
  ),

  tar_target(
    simulation_split,
    {
      set.seed(SIMULATION_SPLIT_SEED)
      split_train_test(ncol(simulation_expression), train_frac = SIMULATION_TRAIN_FRACTION)
    }
  ),

  tar_target(
    simulation_data_train,
    build_simulation_dataset(
      expression_matrix = simulation_expression,
      survival_df = simulation_survival,
      column_indices = simulation_split$train,
      dataname = paste0(SIMULATION_DATASET_LABEL, "_train")
    )
  ),

  tar_target(
    simulation_data_val,
    build_simulation_dataset(
      expression_matrix = simulation_expression,
      survival_df = simulation_survival,
      column_indices = simulation_split$test,
      dataname = paste0(SIMULATION_DATASET_LABEL, "_val")
    )
  ),

  tar_target(
    data,
    simulation_data_train
  ),

  tar_target(
    data_val,
    {
      val_dataset <- simulation_data_val
      setNames(list(val_dataset), val_dataset$dataname)
    }
  )
)

c(
  targets_list,
  COMMON_DESURV_TARGETS
)
