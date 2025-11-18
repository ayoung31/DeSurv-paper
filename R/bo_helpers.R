standardize_bo_params <- function(params) {
  if (is.null(params)) {
    return(list())
  }
  params_list <- as.list(params)
  if (length(params_list) == 0) {
    return(params_list)
  }
  names(params_list) <- sub("_grid$", "", names(params_list))
  params_list
}

maybe_add_numeric_bound <- function(bounds,
                                    values,
                                    name,
                                    type = c("continuous", "integer"),
                                    log_scale = FALSE) {
  type <- match.arg(type)
  if (length(values) == 0) {
    stop(sprintf("No values provided for optional bound '%s'.", name), call. = FALSE)
  }
  vals <- unique(as.numeric(values))
  if (length(vals) <= 1) {
    return(bounds)
  }
  lower <- min(vals)
  upper <- max(vals)
  spec <- list(lower = lower, upper = upper, type = type)
  if (isTRUE(log_scale) && lower > 0) {
    spec$scale <- "log10"
  }
  bounds[[name]] <- spec
  bounds
}
