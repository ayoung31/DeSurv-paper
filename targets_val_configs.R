targets_val_configs <- function() {
  list(
    # easy = list(
    #   run_key = "easy",
    #   mode = "external",
    #   val_datasets = c(
    #     "CPTAC",
    #     "Dijk",
    #     "Moffitt_GEO_array",
    #     "PACA_AU_array",
    #     "PACA_AU_seq",
    #     "Puleo_array"
    #   )
    # )
    
    full = list(
      run_key = "full",
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
    # default = list(
    #   run_key = "default",
    #   mode = "external",
    #   val_datasets = c(
    #     "CPTAC",
    #     "Dijk",
    #     "Moffitt_GEO_array",
    #     "PACA_AU_array",
    #     "PACA_AU_seq",
    #     "Puleo_array"
    #   )
    # )
    # bladder_holdout = list(
    #   run_key = "bladder_default",
    #   mode = "train_split",
    #   val_datasets = character(0)
    # )
  )
}
