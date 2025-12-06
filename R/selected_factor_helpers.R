sanitize_selection_flags <- function(flags) {
  if (is.null(flags)) {
    return(logical())
  }
  if (is.logical(flags)) {
    clean <- flags
  } else if (is.numeric(flags)) {
    clean <- as.logical(flags)
  } else {
    clean <- tolower(trimws(as.character(flags))) %in% c("true", "t", "1", "yes", "y")
  }
  clean[is.na(clean)] <- FALSE
  clean
}

clean_factor_ids <- function(facs) {
  if (is.null(facs)) {
    return(integer())
  }
  facs <- facs[!is.na(facs)]
  if (!length(facs)) {
    return(integer())
  }
  if (!is.numeric(facs)) {
    suppressWarnings(facs <- as.integer(as.character(facs)))
  } else {
    facs <- as.integer(facs)
  }
  facs <- facs[!is.na(facs) & facs > 0]
  unique(facs)
}

read_selected_factor_indices <- function(path, label = "factor selection table") {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop(sprintf("Path to %s is missing; run the selection target first.", label))
  }
  if (!file.exists(path)) {
    stop(sprintf("Could not find %s at %s. Run the selection target first.", label, path))
  }
  
  tbl <- read.csv(path, stringsAsFactors = FALSE)
  if (!("factor" %in% names(tbl))) {
    stop(sprintf("%s at %s is missing a `factor` column.", label, path))
  }
  if (!("selected" %in% names(tbl))) {
    stop(sprintf("%s at %s is missing a `selected` column.", label, path))
  }
  
  flags <- sanitize_selection_flags(tbl$selected)
  if (length(flags) != nrow(tbl)) {
    stop(sprintf("`selected` column in %s (%s) is inconsistent with number of rows.", label, path))
  }
  idx <- which(flags)
  if (!length(idx)) {
    stop(sprintf("No factors marked as selected in %s (%s). Please mark at least one entry as selected.", label, path))
  }
  
  facs <- clean_factor_ids(tbl$factor[idx])
  if (!length(facs)) {
    stop(sprintf("Selected factors in %s (%s) could not be parsed into positive integer indices.", label, path))
  }
  facs
}
