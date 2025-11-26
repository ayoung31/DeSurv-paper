NINIT <- 20
NINIT_FULL <- 100
BO_N_INIT <- 10
BO_N_ITER <- 10
BO_CANDIDATE_POOL <- 2000
BO_MAX_REFINEMENTS <- 2
BO_TOL_GAIN <- 0.002
BO_PLATEAU <- 1
BO_TOP_K <- 5
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
BO_REFINE_CONTROL <- list(
  n_init = BO_N_INIT,
  n_iter = BO_N_ITER,
  candidate_pool = BO_CANDIDATE_POOL,
  exploration_weight = 0.01,
  seed = 456,
  cv_verbose = FALSE
)

# Each scenario may optionally provide `designs = list(...)` where each entry is a
# named list of simulator arguments (or a list with fields `label`/`args`). When
# omitted, the `args` list below is used for all runs.
SIMULATION_SCENARIOS <- list(
  shared_baseline = list(
    simulator = "simulate_desurv_data_shared_baseline",
    args = list(G = 5000, N = 200, K = 4),
    designs = list(
      default = list(),
      large_cohort = list(
        N = 400,
        lethal_prog = 2,
        lethal_effect = 2.0
      ),
      noisy_programs = list(
        noise_sd_baseline = 0.35,
        bump_sd = 1.0,
        bump_mean = 4
      )
    ),
    n_reps = 5L,
    expression_transform = "rank",
    base_seed = 2024L,
    split_seed = 4048L
  )
)


SIMULATION_RESULTS_ROOT <- file.path("results", "simulations_bo")
SIMULATION_DEFAULT_DATASET <- "SIMULATION"
USE_TRAIN_GENES_FOR_VAL <- TRUE

TRAIN_DATASETS <- SIMULATION_DEFAULT_DATASET
TRAIN_PREFIX <- SIMULATION_DEFAULT_DATASET

DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 8L, type = "integer"),
  alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)
NGENE_CONFIG <- c(1000, 5000)
NTOP_CONFIG <- c(25, 200)
LAMBDAW_CONFIG <- c(0)
LAMBDAH_CONFIG <- c(1e-3, 1e3)

SIMULATION_SCENARIOS_ALPHA0 <- purrr::imap(
  SIMULATION_SCENARIOS,
  function(cfg, name) {
    cfg_alpha0 <- rlang::duplicate(cfg, shallow = FALSE)
    base_bounds <- cfg_alpha0$bo_bounds
    if (is.null(base_bounds)) {
      base_bounds <- DESURV_BO_BOUNDS
    }
    base_bounds$alpha_grid <- NULL
    cfg_alpha0$bo_bounds <- base_bounds
    extra_args <- cfg_alpha0$bo_extra_args
    if (is.null(extra_args)) {
      extra_args <- list()
    }
    extra_args$alpha_grid <- 0
    cfg_alpha0$bo_extra_args <- extra_args
    cfg_alpha0
  }
)

SIMULATION_SCENARIOS <- c(
  SIMULATION_SCENARIOS,
  purrr::set_names(
    SIMULATION_SCENARIOS_ALPHA0,
    paste0(names(SIMULATION_SCENARIOS_ALPHA0), "_alpha0")
  )
)

source("targets_setup.R")
purrr::walk(
  list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$"),
  source
)

SIMULATION_BO_TARGETS <- list(
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
    simulation_bo_replicates,
    run_simulation_replicate(simulation_replicate_specs),
    pattern = map(simulation_replicate_specs),
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "cv")
    )
  ),
  tar_target(
    simulation_bo_runs,
    assemble_simulation_config_results(
      config_spec = simulation_config_specs,
      replicate_results = simulation_bo_replicates
    ),
    pattern = map(simulation_config_specs),
    iteration = "list"
  ),
  tar_target(
    simulation_bo_summaries,
    summarize_simulation_config(simulation_bo_runs),
    pattern = map(simulation_bo_runs),
    iteration = "list"
  )
)

SIMULATION_BO_TARGETS
