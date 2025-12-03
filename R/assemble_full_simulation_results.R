assemble_full_simulation_results <- function(config_spec, replicate_results) {
  replicate_results <- replicate_results %||% list()
  reps <- purrr::keep(
    replicate_results,
    ~ identical(.x$config_signature, config_spec$config_signature)
  )
  list(
    config_name = config_spec$config_name,
    design_name = config_spec$design_name,
    config_signature = config_spec$config_signature,
    replicates = reps
  )
}
