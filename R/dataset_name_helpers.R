infer_validation_dataset_name <- function(data) {
  name <- NULL
  if (!is.null(data$dataname) && length(data$dataname)) {
    name <- data$dataname[[1]]
  }
  if (is.null(name) || is.na(name) || !nzchar(trimws(as.character(name)))) {
    if (!is.null(data$sampInfo) && !is.null(data$sampInfo$dataset)) {
      candidates <- unique(data$sampInfo$dataset)
      candidates <- candidates[!is.na(candidates)]
      if (length(candidates)) {
        name <- candidates[[1]]
      }
    }
  }
  if (is.null(name) || is.na(name) || !nzchar(trimws(as.character(name)))) {
    name <- "validation"
  }
  trimws(as.character(name))
}

unwrap_validation_dataset <- function(data) {
  if (is.list(data) && length(data) == 1L && is.list(data[[1]]) && !is.null(data[[1]]$ex)) {
    return(data[[1]])
  }
  data
}
