# _targets_sims.R - Local Slurm version for desktop (20 CPUs, 31GB RAM)
# Key difference: cpus_per_task is decoupled from SIM_CV_NSTARTS to fit local resources
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(DeSurv)
})

sim_files <- list.files(
  "R/simulation_functions",
  pattern = "[.]R$",
  full.names = TRUE
)
purrr::walk(sim_files, source)
source("R/get_top_genes.R")



`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

standardize_bayes_params <- function(params) {
  if (is.null(params)) {
    return(list())
  }
  out <- as.list(params)
  if (!length(out)) {
    return(out)
  }
  names(out) <- sub("_grid$", "", names(out))
  out
}

SIM_DATASETS_PER_SCENARIO <- 100L
SIM_GLOBAL_SEED <- 101L
SIM_METHOD_TRANSFORM <- "none"
SIM_DESURV_TOL <- 1e-5
SIM_DESURV_MAXIT <- 3000L
# SIM_FINAL_NINIT <- 1L
SIM_DEFAULT_NGENE <- NULL
SIM_DEFAULT_NTOP <- 150L
SIM_DEFAULT_LAMBDAW <- 0
SIM_DEFAULT_LAMBDAH <- 0
SIM_DEFAULT_K <- 3L
SIM_DEFAULT_ALPHA <- 0.6
SIM_DEFAULT_LAMBDA <- 0.1
SIM_DEFAULT_NU <- 0.3
SIM_BETA_NONZERO_TOL <- 1e-8
SIM_CV_NFOLDS <- 5L
SIM_CV_NSTARTS <- 30  # Number of initializations (unchanged from original)
SIM_TRAIN_FRACTION <- 0.7
SIM_SPLIT_SEED_OFFSET <- 10000L

# Local desktop: cpus_per_task limited to 19 (decoupled from SIM_CV_NSTARTS)
# The parallel code will schedule 30 inits across 19 CPUs
SIM_LOCAL_CPUS <- 19L

SIM_ANALYSIS_CONTROLLER <- crew_controller_slurm(

  name = "sim_analysis",
  workers = 20,
  options_cluster = crew_options_slurm(
    cpus_per_task = SIM_LOCAL_CPUS,  # Local limit, not SIM_CV_NSTARTS
    memory_gigabytes_required = 8,
    time_minutes = 240,
    log_error = "logs/crew_log_%A.err",
    log_output = "logs/crew_log_%A.out"
  )
)
default_controller = crew_controller_sequential()

active_controller <- crew_controller_group(default_controller,
                                           SIM_ANALYSIS_CONTROLLER)

tar_option_set(
  packages = c("dplyr", "purrr", "tibble", "DeSurv", "survival"),
  format = "rds",
  controller = active_controller,
  error = "continue"
)


SIM_DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 5L, type = "integer"),
  alpha_grid = list(lower = 0, upper = .95, type = "continuous"),
  lambda_grid = list(lower = 1e-2, upper = 1e2, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

SIM_BO_N_INIT <- 10L
SIM_BO_N_ITER <- 20L
SIM_BO_CANDIDATE_POOL <- 1500L
SIM_BO_EXPLORATION_WEIGHT <- 0.01

SIMULATION_SCENARIOS <- list(
  list(
    scenario_id = "R0_easy",
    scenario = "R0",
    description = "Default easy/sanity scenario",
    replicates = SIM_DATASETS_PER_SCENARIO,
    seed_offset = SIM_GLOBAL_SEED,
    overrides = list()
  ),
  list(
    scenario_id = "R00_null",
    scenario = "R00",
    description = "No survival associated programs",
    replicates = SIM_DATASETS_PER_SCENARIO,
    seed_offset = SIM_GLOBAL_SEED,
    overrides = list()
  )#,
  # list(
  #   scenario_id = "R2_correlated",
  #   scenario = "R2",
  #   description = "Correlated programs with moderate noise",
  #   replicates = SIM_DATASETS_PER_SCENARIO,
  #   seed_offset = SIM_GLOBAL_SEED + 1000L,
  #   overrides = list()
  # ),
  # list(
  #   scenario_id = "R3_overlap",
  #   scenario = "R3",
  #   description = "Overlapping markers and strong background",
  #   replicates = SIM_DATASETS_PER_SCENARIO,
  #   seed_offset = SIM_GLOBAL_SEED + 2000L,
  #   overrides = list()
  # )
)

SIM_FIXED_PARAMS <- list(
  k = NULL,
  alpha = SIM_DEFAULT_ALPHA,
  lambda = SIM_DEFAULT_LAMBDA,
  nu = SIM_DEFAULT_NU,
  lambdaW = SIM_DEFAULT_LAMBDAW,
  lambdaH = SIM_DEFAULT_LAMBDAH,
  ngene = SIM_DEFAULT_NGENE
)

SIM_ANALYSIS_SPECS <- list(
  list(
    analysis_id = "fixed",
    label = "Fixed parameters",
    mode = "fixed",
    params = SIM_FIXED_PARAMS
  ),
  list(
    analysis_id = "fixed_alpha0",
    label = "Fixed with alpha = 0",
    mode = "fixed",
    params = modifyList(SIM_FIXED_PARAMS, list(alpha = 0))
  ),
  list(
    analysis_id = "bo",
    label = "Bayesian optimization",
    mode = "bayesopt",
    bounds = SIM_DESURV_BO_BOUNDS,
    bo_fixed = list(
      ngene = SIM_DEFAULT_NGENE,
      ntop = SIM_DEFAULT_NTOP,
      lambdaW_grid = SIM_DEFAULT_LAMBDAW,
      lambdaH_grid = SIM_DEFAULT_LAMBDAH
    )
  ),
  list(
    analysis_id = "bo_alpha0",
    label = "Bayesian optimization with alpha = 0",
    mode = "bayesopt",
    bounds = modifyList(
      SIM_DESURV_BO_BOUNDS,
      list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"))
    ),
    bo_fixed = list(
      ngene = SIM_DEFAULT_NGENE,
      ntop = SIM_DEFAULT_NTOP,
      alpha_grid = 0,
      lambdaW_grid = SIM_DEFAULT_LAMBDAW,
      lambdaH_grid = SIM_DEFAULT_LAMBDAH
    ),
    final_overrides = list(alpha = 0)
  ),
  list(
    analysis_id = "bo_tune_ntop",
    label = "Bayesian optimization adding ntop to tuning",
    mode = "bayesopt",
    bounds = modifyList(
      SIM_DESURV_BO_BOUNDS,
      list(ntop = list(lower = 50, upper = 250, type = "integer"))),
    bo_fixed = list(
      ngene = SIM_DEFAULT_NGENE,
      lambdaW_grid = SIM_DEFAULT_LAMBDAW,
      lambdaH_grid = SIM_DEFAULT_LAMBDAH
    )
  ),
  list(
    analysis_id = "bo_tune_ntop_alpha0",
    label = "Bayesian optimization with alpha = 0 adding ntop for tuning",
    mode = "bayesopt",
    bounds = modifyList(
      SIM_DESURV_BO_BOUNDS,
      list(alpha_grid = list(lower = 0, upper = 0, type = "continuous"),
           ntop = list(lower = 50, upper = 250, type = "integer"))
    ),
    bo_fixed = list(
      ngene = SIM_DEFAULT_NGENE,
      alpha_grid = 0,
      lambdaW_grid = SIM_DEFAULT_LAMBDAW,
      lambdaH_grid = SIM_DEFAULT_LAMBDAH
    ),
    final_overrides = list(alpha = 0)
  )
)

# NOTE: The rest of this file should be sourced from the main _targets_sims.R
# This local version only overrides the controller configuration at the top.
# To use: copy the remaining functions and targets_list from the original file,
# or source this file first, then source the original (excluding the controller setup).

message("Local Slurm configuration loaded:")
message("  - SIM_CV_NSTARTS (initializations): ", SIM_CV_NSTARTS)
message("  - SIM_LOCAL_CPUS (Slurm cpus_per_task): ", SIM_LOCAL_CPUS)
message("  - workers: 20")
message("  - No module load required")
