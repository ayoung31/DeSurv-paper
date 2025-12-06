sanitize_path_component <- function(x, label = "path component") {
  if (length(x) != 1) {
    stop(
      sprintf(
        "%s must be length 1, but received a vector of length %d.",
        label,
        length(x)
      )
    )
  }
  if (is.null(x) || is.na(x)) {
    stop(sprintf("%s cannot be NULL or NA.", label))
  }
  value <- trimws(as.character(x))
  if (!nzchar(value)) {
    stop(sprintf("%s cannot be empty.", label))
  }
  value <- gsub("[[:cntrl:]]", "", value)
  value <- gsub("[/\\\\:*?\"<>|]+", "_", value)
  value
}
