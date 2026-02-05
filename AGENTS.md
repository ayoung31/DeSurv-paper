# Repository Guidelines

## Project Structure & Module Organization
- `R/`: core analysis helpers and pipeline modules; simulation helpers live under `R/simulation_functions/`.
- `_targets.R`, `_targets_bladder.R`, `_targets_sims.R`: main `targets` pipelines for TCGA/CPTAC, bladder, and simulation runs.
- `targets_setup.R`, `targets_common_pipeline.R`, `targets_configs.R`: shared pipeline setup, common targets, and run/BO/validation configs.
- `tests/testthat/`, `tests/testthat.R`: unit tests for helpers and pipeline behavior.
- `paper/`: R Markdown sources and rendered outputs; `figures/` holds figures used in the manuscript.
- `data/`: input datasets; `_targets/` is the pipeline cache/store; `inst/` contains helper scripts.

## Build, Test, and Development Commands
- `Rscript -e 'devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)'`: install the companion `DeSurv` package required by the pipelines.
- `sbatch _targets.sh` (or `_targets_bladder.sh`, `_targets_sims.sh`): submit pipelines to Slurm (required for distributed execution).
- `Rscript -e 'targets::tar_make(script = "_targets_bladder.R")'`: run a specific pipeline file (requires Slurm environment).
- `Rscript -e 'testthat::test_dir("tests/testthat")'`: run unit tests.

## Coding Style & Naming Conventions
- R code uses 2-space indentation and snake_case for functions and files (match existing `R/*.R`).
- Keep pipeline configuration edits in `targets_configs.R` rather than embedding constants in scripts.
- No formatter/linter is enforced; follow surrounding style and keep helper functions small and composable.

## Testing Guidelines
- Tests use `testthat` and live in `tests/testthat/` with filenames like `test_*.R`.
- Add tests for new helpers or pipeline logic (e.g., configuration validation, bounds helpers).

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and lowercase (e.g., "add flexibility to run pipeline").
- PRs should describe pipeline/config changes, list affected targets, and note expected outputs or figures.
- Call out any long-running steps or data dependencies so reviewers can plan reruns.

## Configuration & Dependency Notes
- The pipelines assume a local checkout of `../DeSurv` and will load it automatically when present.
- ORA/KEGG helpers require `clusterProfiler` and `org.Hs.eg.db` when enrichment targets are enabled.
