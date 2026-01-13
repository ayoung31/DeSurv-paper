targets_val_configs <- function() {
  bo_labels <- names(targets_bo_configs())
  base_val <- list(
    mode = "external",
    val_datasets = c(
      "CPTAC",
      "Dijk",
      "Moffitt_GEO_array",
      "PACA_AU_array",
      "PACA_AU_seq",
      "Puleo_array"
    )
  )
  configs <- lapply(bo_labels, function(label) base_val)
  names(configs) <- bo_labels
  configs
}
