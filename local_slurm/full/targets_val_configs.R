# targets_val_configs.R - Local desktop FULL mode
# Same validation datasets as HPC version (unchanged from quick mode)

targets_val_config <- function(label) {
  list(
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
}

targets_val_configs <- function(bo_labels = names(targets_bo_configs())) {
  if (is.null(bo_labels)) {
    bo_labels <- character(0)
  }
  configs <- lapply(bo_labels, targets_val_config)
  names(configs) <- bo_labels
  configs
}
