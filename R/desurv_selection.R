arrange_desurv_preference <- function(df) {
  if (!all(c("k", "alpha", "lambda", "nu", "lambdaW", "lambdaH") %in% names(df))) {
    stop("Summary data frame is missing required hyperparameter columns.")
  }
  ord <- order(df$k, df$alpha, -df$lambda, df$nu, df$lambdaW, df$lambdaH)
  df[ord, , drop = FALSE]
}

select_one_se <- function(summary_df) {
  if (!all(c("mean_cindex", "se_cindex") %in% names(summary_df))) {
    stop("Summary data frame must contain `mean_cindex` and `se_cindex`.")
  }
  best <- summary_df[which.max(summary_df$mean_cindex), , drop = FALSE]
  if (nrow(best) == 0L) {
    stop("No rows available to select from.")
  }

  if (is.na(best$se_cindex)) {
    return(best)
  }

  threshold <- best$mean_cindex - best$se_cindex
  eligible <- summary_df[summary_df$mean_cindex >= threshold, , drop = FALSE]
  arrange_desurv_preference(eligible)[1, , drop = FALSE]
}
