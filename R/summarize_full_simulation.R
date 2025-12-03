`%||%` <- function(x, y) if (is.null(x)) y else x

summarize_full_simulation <- function(config_results) {
  reps <- config_results$replicates %||% list()
  empty_summary <- data.frame(
    scenario = config_results$config_name,
    design = config_results$design_name,
    model = "desurv",
    n_reps = 0L,
    mean_test_cindex = NA_real_,
    sd_test_cindex = NA_real_,
    mean_reconstruction_error = NA_real_,
    sd_reconstruction_error = NA_real_,
    mean_effect_correlation = NA_real_,
    sd_effect_correlation = NA_real_,
    column_cor = I(list(data.frame())),
    factor_detection = I(list(data.frame())),
    gene_set_recovery = I(list(data.frame())),
    global_gene_set_recovery = I(list(data.frame())),
    stringsAsFactors = FALSE
  )
  if (!length(reps)) {
    return(list(summary = empty_summary, replicates = data.frame()))
  }
  normalize_table_metric <- function(x) {
    if (is.null(x)) {
      return(data.frame())
    }
    x
  }
  rep_metrics <- lapply(
    reps,
    function(rep) {
      models <- rep$model_results
      if (is.null(models) || !length(models)) {
        models <- list(list(
          model = "desurv",
          metrics = rep$metrics,
          bo_best_score = rep$bo_best_score
        ))
      }
      rows <- lapply(
        models,
        function(model_entry) {
          if (is.null(model_entry)) {
            return(NULL)
          }
          metrics <- model_entry$metrics %||% list()
          scalar_or_na <- function(value) {
            if (is.null(value) || !length(value)) {
              return(NA_real_)
            }
            as.numeric(value)[1]
          }
          data.frame(
            scenario = rep$config_name,
            design = rep$design_name,
            replicate = rep$replicate,
            model = model_entry$model %||% "desurv",
            test_cindex = scalar_or_na(metrics$c_index),
            reconstruction_error = scalar_or_na(metrics$reconstruction),
            effect_correlation = scalar_or_na(metrics$effect_correlation),
            bo_best_score = model_entry$bo_best_score %||% rep$bo_best_score %||% NA_real_,
            column_cor = I(list(normalize_table_metric(metrics$column_cor))),
            factor_detection = I(list(normalize_table_metric(metrics$factor_detection))),
            gene_set_recovery = I(list(normalize_table_metric(metrics$gene_set_recovery))),
            global_gene_set_recovery = I(list(normalize_table_metric(metrics$global_gene_set_recovery))),
            stringsAsFactors = FALSE
          )
        }
      )
      rows <- rows[!vapply(rows, is.null, logical(1))]
      if (!length(rows)) {
        return(data.frame(
          scenario = character(),
          design = character(),
          replicate = integer(),
          model = character(),
          test_cindex = numeric(),
          reconstruction_error = numeric(),
          effect_correlation = numeric(),
          bo_best_score = numeric(),
          column_cor = I(list()),
          factor_detection = I(list()),
          gene_set_recovery = I(list()),
          global_gene_set_recovery = I(list()),
          stringsAsFactors = FALSE
        ))
      }
      do.call(rbind, rows)
    }
  )
  rep_table <- do.call(rbind, rep_metrics)
  if (is.null(rep_table)) {
    rep_table <- data.frame(
      scenario = character(),
      design = character(),
      replicate = integer(),
      model = character(),
      test_cindex = numeric(),
      reconstruction_error = numeric(),
      effect_correlation = numeric(),
      bo_best_score = numeric(),
      column_cor = I(list()),
      factor_detection = I(list()),
      gene_set_recovery = I(list()),
      global_gene_set_recovery = I(list()),
      stringsAsFactors = FALSE
    )
  }
  if (!nrow(rep_table)) {
    return(list(summary = empty_summary, replicates = rep_table))
  }
  safe_mean <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  }
  safe_sd <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) <= 1) {
      return(NA_real_)
    }
    stats::sd(x)
  }
  bind_metric_tables <- function(tbl, colname) {
    if (!nrow(tbl)) {
      return(data.frame())
    }
    entries <- tbl[[colname]]
    if (is.null(entries) || !length(entries)) {
      return(data.frame())
    }
    combined <- lapply(seq_along(entries), function(i) {
      entry <- entries[[i]]
      if (is.null(entry) || !nrow(entry)) {
        return(NULL)
      }
      entry$scenario <- tbl$scenario[[i]]
      entry$design <- tbl$design[[i]]
      entry$model <- tbl$model[[i]]
      entry$replicate <- tbl$replicate[[i]]
      entry
    })
    combined <- combined[!vapply(combined, is.null, logical(1))]
    if (!length(combined)) {
      return(data.frame())
    }
    do.call(rbind, combined)
  }
  split_models <- split(rep_table, rep_table$model)
  summary_rows <- lapply(
    split_models,
    function(tbl) {
      data.frame(
        scenario = tbl$scenario[[1]],
        design = tbl$design[[1]],
        model = tbl$model[[1]],
        n_reps = nrow(tbl),
        mean_test_cindex = safe_mean(tbl$test_cindex),
        sd_test_cindex = safe_sd(tbl$test_cindex),
        mean_reconstruction_error = safe_mean(tbl$reconstruction_error),
        sd_reconstruction_error = safe_sd(tbl$reconstruction_error),
        mean_effect_correlation = safe_mean(tbl$effect_correlation),
        sd_effect_correlation = safe_sd(tbl$effect_correlation),
        column_cor = I(list(bind_metric_tables(tbl, "column_cor"))),
        factor_detection = I(list(bind_metric_tables(tbl, "factor_detection"))),
        gene_set_recovery = I(list(bind_metric_tables(tbl, "gene_set_recovery"))),
        global_gene_set_recovery = I(list(bind_metric_tables(tbl, "global_gene_set_recovery"))),
        stringsAsFactors = FALSE
      )
    }
  )
  summary <- do.call(rbind, summary_rows)
  list(summary = summary, replicates = rep_table)
}

