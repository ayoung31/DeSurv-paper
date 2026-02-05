# Combine PDAC datasets from raw file inputs and return a formatted data list.
load_data_pdac_split <- function(raw_files) {
  if (is.null(raw_files) || !length(raw_files)) {
    stop("raw_files must be a non-empty character vector.")
  }
  if (!is.character(raw_files)) {
    stop("raw_files must be a character vector of file paths.")
  }

  base_names <- basename(raw_files)
  datasets <- sub("\\.survival_data\\.rds$", "", base_names)
  datasets <- sub("_subtype\\.csv$", "", datasets)
  datasets <- sub("\\.rds$", "", datasets)
  datasets <- unique(datasets)
  if (!length(datasets)) {
    stop("No dataset names could be inferred from raw_files.")
  }

  load_data(datasets)
}
