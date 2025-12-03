NINIT <- 5
NINIT_FULL <- 10
BO_N_INIT <- 10
BO_N_ITER <- 20
BO_CANDIDATE_POOL <- 200
BO_MAX_REFINEMENTS <- 1
BO_TOL_GAIN <- 0.01
BO_PLATEAU <- 1
BO_TOP_K <- 10
BO_SHRINK_BASE <- 0.5
BO_IMPORTANCE_GAIN <- 0.3
BO_COARSE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 123,
  cv_verbose = FALSE
)
BO_REFINE_CONTROL <- BO_COARSE_CONTROL

SIMULATION_SCENARIOS <- list(
  easy = list(
    simulator = "simulate_desurv_easy",
    design_name="simulate_desurv_easy",
    args = list(
      G = 5000,
      N = 500,
      K = 5,
      big_prog = 2,
      correlated_pairs = list(c(1, 2), c(3, 4))
    ),
    n_reps = 5L,
    expression_transform = "log2",
    train_fraction = 0.7,
    base_seed = 1001L,
    split_seed = 2001L
  )#,
  # nuisance = list(
  #   simulator = "simulate_desurv_nuisance",
  #   design_name = "simulate_desurv_nuisance",
  #   args = list(
  #     G = 5000,
  #     N = 500,
  #     K = 5,
  #     subtle_prog = 2,
  #     correlated_pairs = list(c(1, 2), c(3, 4))
  #   ),
  #   n_reps = 5L,
  #   expression_transform = "log2",
  #   train_fraction = 0.7,
  #   base_seed = 3001L,
  #   split_seed = 4001L
  # )
)

SIMULATION_RESULTS_ROOT <- file.path("results", "simulations_full")
SIMULATION_DEFAULT_DATASET <- "SIM_FULL"
TRAIN_DATASETS <- SIMULATION_DEFAULT_DATASET
TRAIN_PREFIX <- SIMULATION_DEFAULT_DATASET
USE_TRAIN_GENES_FOR_VAL <- TRUE

DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 10L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-4, upper = 1e2, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)
NGENE_CONFIG <- c(1000,5000)
NTOP_CONFIG <- c(50,200)
LAMBDAW_CONFIG <- c(0)
LAMBDAH_CONFIG <- c(1e-5)

source("targets_setup.R")
purrr::walk(
  list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)
METHOD_TRANS_TRAIN="none"
source("targets_common_pipeline.R")
tar_option_set(packages = setdiff(tar_option_get("packages"), "dplyr"))
Sys.setenv(CREW_HOST = Sys.getenv("CREW_HOST", unset = "127.0.0.1"))

SIM_FULL_TARGETS <- list(
  tar_target(
    simulation_config_specs,
    build_simulation_config_list(SIMULATION_SCENARIOS),
    iteration = "list"
  ),
  tar_target(
    simulation_replicate_specs,
    expand_all_simulation_replicates(simulation_config_specs),
    iteration = "list"
  ),
  tar_target(
    simulation_full_replicates,
    run_full_simulation_replicate(simulation_replicate_specs),
    pattern = map(simulation_replicate_specs),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  tar_target(
    simulation_full_config_results,
    assemble_full_simulation_results(simulation_config_specs, simulation_full_replicates),
    pattern = map(simulation_config_specs),
    iteration = "list"
  ),
  tar_target(
    simulation_full_config_summary,
    {
      browser()
      summarize_full_simulation(simulation_full_config_results)
    },
    pattern = map(simulation_full_config_results),
    iteration = "list"
  ),
  tar_target(
    simulation_full_summary_table,
    aggregate_summary_tables(simulation_full_config_summary)
  )
)

SIM_FULL_TARGETS
