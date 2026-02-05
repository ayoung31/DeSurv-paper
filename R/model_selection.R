#' Generate metrics summary file path
build_metrics_path <- function(version, ngene, k, ninit, imaxit, tol, maxit) {
  file.path('results', 
            version, 
            'full', 
            paste0('ngene', ngene), 
            'summary',
            paste0('metrics_ngene', ngene, '_k=', k,
                   '_ninit=', ninit, '_imaxit=', imaxit,
                   '_tol=', tol, '_maxit=', maxit, '.RData'))
}

#' Generate model filename (raw result file)
build_model_filename <- function(k, alpha, lambda, eta, lambdaW, lambdaH,
                                 ninit, imaxit, tol, maxit) {
  paste0("k=", k, 
         "_alpha", alpha,
         "_lambda", lambda, 
         "_eta", eta,
         "_lambdaW", lambdaW, 
         "_lambdaH", lambdaH,
         "_full_ninit", ninit, 
         "_imaxit", imaxit,
         "_tol", tol, 
         "_maxit", maxit, ".RData")
}


#' Generate full path to the raw model
build_model_path <- function(version, ngene, filename) {
  file.path('results', 
            version, 
            'full', 
            paste0('ngene', ngene), 
            'raw', 
            filename)
}


#' Build (and create) the save directory for a selected model
build_model_save_dir <- function(version, ngene, ninit, imaxit, tol, maxit,
                                 k, alpha, lambda, eta, lambdaW, lambdaH, ntop) {
  path <- file.path(
    "results",
    version,
    "full",
    paste0("ngene", ngene),
    "selected models",
    paste0(
      "ninit=", ninit,
      "_imaxit=", imaxit,
      "_tol=", tol,
      "_maxit=", maxit
    ),
    "figures",
    paste0(
      "k=", k,
      "_alpha",  alpha,
      "_lambda", lambda,
      "_eta",    eta,
      "_lambdaW",lambdaW,
      "_lambdaH",lambdaH
    ),
    paste0("top", ntop)
  )
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  return(path)
}


#' Select the best model (or a model at a specified alpha) and prepare results
#'
#' @description
#' Loads a metrics table, selects a model according to `which.alpha`, loads the
#' corresponding fitted object, and returns a results list with paths and metadata.
#'
#' @param model.params list
#'   A list containing (at minimum) the fields:
#'   `dataset`, `version`, `ngene`, `ninit`, `imaxit`, `tol`, `maxit`, `k`.
#' @param which.alpha character or numeric
#'   Either `"best"` to select the best row by the score `c / nloss`,
#'   or a numeric value (including `0`) to select rows with `alpha == which.alpha`.
#' @param ntop integer
#'   The number of top genes (or features) to carry forward in downstream steps.
#' @details
#' - The metrics file is expected to contain a data frame named `metrics` with
#'   columns `alpha`, `lambda`, `eta`, `lambdaW`, `lambdaH`, `c`, and `nloss`.
#' - The following helper functions are assumed to exist and return valid paths:
#'   `build_metrics_path(version, ngene, k, ninit, imaxit, tol, maxit)`,
#'   `build_model_filename(k, alpha, lambda, eta, lambdaW, lambdaH, ninit, imaxit, tol, maxit)`,
#'   `build_model_path(version, ngene, model_filename)`,
#'   `build_model_save_dir(version, ngene, ninit, imaxit, tol, maxit,
#'                         k, alpha, lambda, eta, lambdaW, lambdaH, ntop)`.
#'
#' @return list
#'   A named list with elements:
#'   - `ntop`, `fit_cox`, `model.params`, `alpha`
#'   - `model_save_dir`: directory where outputs for this selection should be saved
#'   - `selection`: the single-row data frame of chosen hyperparameters/metrics
#'
#' @examples
#' \dontrun{
#' res <- select_model(model.params, which.alpha = "best", ntop = 100)
#' }
#'
#' @export
select_model <- function(model.params, which.alpha, ntop) {
  # ---- Validate inputs -------------------------------------------------------
  required_fields <- c("dataset","version","ngene","ninit","imaxit","tol","maxit","k")
  missing_fields <- setdiff(required_fields, names(model.params))
  if (length(missing_fields) > 0) {
    stop("`model.params` is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }

  if (!is.numeric(ntop) || length(ntop) != 1L || ntop < 1) {
    stop("`ntop` must be a single positive integer.")
  }
  
  version <- model.params$version
  ngene   <- model.params$ngene
  ninit   <- model.params$ninit
  imaxit  <- model.params$imaxit
  tol     <- model.params$tol
  maxit   <- model.params$maxit
  k       <- model.params$k
  
  # ---- Load metrics safely ---------------------------------------------------
  metrics_path <- build_metrics_path(version, ngene, k, ninit, imaxit, tol, maxit)
  if (!file.exists(metrics_path)) stop("Metrics file not found: ", metrics_path)
  
  # prefer readRDS if possible; otherwise load into a private env
  metrics <- tryCatch({
    # If your pipeline writes an .RDS, you can swap to readRDS(metrics_path)
    env <- new.env(parent = emptyenv())
    loaded <- load(metrics_path, envir = env)
    if (!"metrics" %in% loaded && !"metrics" %in% ls(env)) {
      stop("The metrics file does not contain an object named `metrics`.")
    }
    get("metrics", envir = env)
  }, error = function(e) {
    stop("Failed to load metrics from: ", metrics_path, "\n", conditionMessage(e))
  })
  
  # Basic schema checks
  needed_cols <- c("alpha","lambda","eta","lambdaW","lambdaH","c","nloss")
  missing_cols <- setdiff(needed_cols, colnames(metrics))
  if (length(missing_cols) > 0) {
    stop("`metrics` is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # ---- Filter metrics by alpha ----------------------------------------------
  if (is.character(which.alpha)) {
    if (!identical(which.alpha, "best")) {
      stop("Invalid `which.alpha`: use \"best\" or a numeric alpha (e.g., 0).")
    }
    metrics_filtered <- metrics
  } else if (is.numeric(which.alpha) && length(which.alpha) == 1L) {
    metrics_filtered <- metrics[metrics$alpha == which.alpha, , drop = FALSE]
  } else {
    stop("`which.alpha` must be \"best\" or a single numeric value.")
  }
  
  if (nrow(metrics_filtered) == 0) {
    stop("No matching models found for the requested alpha (", which.alpha, ").")
  }
  
  # ---- Choose the row to use -------------------------------------------------
  # Score = c / nloss (guard against zero/NA)
  eps <- 1e-12
  denom <- ifelse(is.na(metrics_filtered$nloss) | metrics_filtered$nloss == 0,
                  eps, metrics_filtered$nloss)
  score <- metrics_filtered$c / denom
  
  # Tie-breaker: prefer lower lambda, then lower eta, then lower lambdaW/H
  ord <- order(score, -metrics_filtered$c,  # favor larger c if same ratio
               -1/denom,                    # smaller nloss is better
               metrics_filtered$lambda,
               metrics_filtered$eta,
               metrics_filtered$lambdaW,
               metrics_filtered$lambdaH,
               decreasing = TRUE)
  best_row <- metrics_filtered[ord[1L], , drop = FALSE]
  
  # ---- Build model filename & path ------------------------------------------
  model_filename <- build_model_filename(
    k            = k,
    alpha        = best_row$alpha,
    lambda       = best_row$lambda,
    eta          = best_row$eta,
    lambdaW      = best_row$lambdaW,
    lambdaH      = best_row$lambdaH,
    ninit        = ninit,
    imaxit       = imaxit,
    tol          = tol,
    maxit        = maxit
  )
  
  model_path <- build_model_path(version, ngene, model_filename)
  if (!file.exists(model_path)) stop("Model file not found: ", model_path)
  
  # ---- Load fitted model safely ---------------------------------------------
  fit_cox <- tryCatch({
    env <- new.env(parent = emptyenv())
    loaded <- load(model_path, envir = env)
    if (!"fit_cox" %in% loaded && !"fit_cox" %in% ls(env)) {
      stop("Model file did not contain an object named `fit_cox`.")
    }
    get("fit_cox", envir = env)
  }, error = function(e) {
    stop("Failed to load fitted model from: ", model_path, "\n", conditionMessage(e))
  })
  
  # ---- Determine the save directory -----------------------------------------
  model_save_dir <- build_model_save_dir(
    version, ngene, ninit, imaxit, tol, maxit,
    k,
    best_row$alpha, best_row$lambda, best_row$eta,
    best_row$lambdaW, best_row$lambdaH, ntop
  )
  if (!dir.exists(model_save_dir)) dir.create(model_save_dir, recursive = TRUE)
  
  # ---- Assemble results ------------------------------------------------------
  results <- list(
    ntop          = ntop,
    fit_cox       = fit_cox,
    model.params  = model.params,
    alpha         = which.alpha,     # what the caller requested ("best" or numeric)
    model_save_dir = model_save_dir,
    selection     = best_row         # the actual chosen hyperparameters/metrics
  )
  
  return(results)
}


model_description = function(model.params,which.alpha,ntop,
                             top_genes,p.adj=.05,colors,plot=TRUE,do_ora=TRUE){
  ####### model at alpha=0 ########
  
  results = select_model(model.params,which.alpha=which.alpha,ntop=ntop)
  
  #if selected model exists, load it
  model_name=file.path(results$model_save_dir,"model_summary.RData")

  ## get top genes
  results = get_top_genes(results)
  results=ORA(results,top_genes,plot=plot,p.adj=p.adj,msigdbr=FALSE)
  
  save(results,file=file.path(results$model_save_dir,"model_summary.RData"))
  
  return(results)
}
