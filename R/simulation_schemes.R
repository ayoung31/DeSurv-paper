#' Simulation schemes highlighting DeSurv advantages
#'
#' Each entry describes a scenario we use when benchmarking methods.
#' The comment explains the biology and why DeSurv should outperform a
#' standard NMF pipeline that ignores survival (alpha = 0).
simulation_schemes <- list(
  shared_baseline = list(
    description = "Housekeeping baseline + lethal program. All tumors share a dominant housekeeping program; only a subtle lethal program distinguishes outcomes." ,
    expectation = "NMF tends to focus on variance (housekeeping). DeSurv prioritizes the lethal program, so it recovers the survival driver faster."
  ),
  large_cohort = list(
    description = "Large sample size with aggressive lethal program. Persists across clones but only survival-guided methods exploit the signal fully.",
    expectation = "Both methods improve with N=400, but DeSurv still converges with fewer BO refinements while NMF wastes capacity fitting neutral variation."
  ),
  noisy_programs = list(
    description = "Shared baseline plus program-specific noise and bumps. Survival depends on a program partially obscured by noise.",
    expectation = "DeSurv leverages survival labels to down-weight noisy programs, whereas alpha=0 NMF often chases noise and misses the hazard driver."
  ),
  easy_single_program = list(
    description = "One dominant malignant program drives both expression and mortality. Useful sanity-check scenario.",
    expectation = "Both recover the lethal program, but DeSurv learns it with fewer samples/iterations and yields better rank-ordering of risk."
  ),
  subtle_program = list(
    description = "Lethal program has tiny amplitude amid dominant benign programs (immune/stromal).",
    expectation = "Variance-based NMF ignores the subtle program; DeSurv, tied to survival, isolates it despite small expression effect sizes."
  ),
  nuisance_dominated = list(
    description = "Subtle lethal program plus large batch/purity shifts.",
    expectation = "DeSurv learns to ignore batch/purity because they are survival-neutral; alpha=0 NMF overfits nuisances and loses predictive power."
  ),
  interaction = list(
    description = "Two programs only become lethal when co-occurring (synergy).",
    expectation = "DeSurv spots the interaction via survival feedback. Standard NMF lacks supervision and fails to connect co-activation with hazard."
  )
)
