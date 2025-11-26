This package contains code to reproduce findings in the DeSurv paper: 
A survival-driven deconvolution tool to prognostic subtype discovery

To run this code you will need to install the R targets package. 
Targets is a pipeline management tool that will run DeSurv from original data to final figures.
```{r}
install.packages("targets")
```

The pipeline now depends on the development version of the `DeSurv` package
located one directory up in `../DeSurv`. Install or update it before running
any targets to pick up the latest C++ and R functionality:

```r
# install.packages("devtools")
devtools::install_local("../DeSurv", upgrade = "never", force = TRUE)
```

Each `_targets*.R` script will automatically `pkgload::load_all("../DeSurv")`
when that directory is present, so the local checkout is always used.

Hyperparameter tuning is handled inside the `DeSurv` package via
`desurv_cv_bayesopt_refine()`, which launches an initial Bayesian optimisation
run, shrinks the promising region, and iterates until improvements plateau
instead of relying on custom grid helpers in this repository. All preprocessing,
validation filtering, and model-selection logic now flow through the package
interface.

The `_targets*.R` scripts expose configuration constants near the top. Two that
you will likely edit are:

- `NGENE_CONFIG`, which controls the gene-filter size passed into preprocessing.
  Provide a single value (e.g., `c(2000)`) to fix the number of genes or a
  range (e.g., `c(2000, 5000)`) to let Bayesian optimisation tune `ngene`
  alongside the other hyperparameters.
- `NTOP_CONFIG`, which controls how many genes per factor are used when
  generating downstream clustering/paper products. As with `NGENE_CONFIG`, set
  one value to keep it fixed or a range to allow BO to optimise it.
- `LAMBDAW_CONFIG`/`LAMBDAH_CONFIG`, which control whether the factor penalty
  terms are tuned. Keep them as singletons (e.g., `c(0)`) to fix the penalties
  at that value or provide positive ranges (for example
  `c(1e-5, 1e-2)`) to let BO search on a log-scale.


### Optional tuning parameters

The configuration constants above feed directly into the Bayesian optimiser.
If you set a vector with more than one unique value, the pipeline will add the
appropriate bound to `desurv_cv_bayesopt_refine()`; otherwise the single value is
passed as a fixed `_grid` argument. The helper `maybe_add_numeric_bound()`
(see `R/bo_helpers.R`) handles this logic and is covered by the unit tests under
`tests/testthat`.

The shared `DESURV_BO_BOUNDS` list defines the non-dataset-specific search space
for the optimiser (k, alpha, lambda, nu). Adjust these entries if the default
range does not match your experiment; the pipeline will merge them with any
`NGENE_CONFIG`/`NTOP_CONFIG`/`LAMBDA*_CONFIG` ranges automatically before
calling `desurv_cv_bayesopt_refine()`. Additional refinement controls live next
to the other global constants in each `_targets*.R`: `BO_MAX_REFINEMENTS`
determines how many shrinks are attempted, while `BO_COARSE_CONTROL` and
`BO_REFINE_CONTROL` describe the per-stage BO settings (e.g., `n_init`, `n_iter`,
`candidate_pool`, `seed`, `cv_verbose`). Update them if you need a longer coarse
search, different seeds, or a quieter optimiser.

### Optional Enrichment Dependencies

The ORA/KEGG helpers in `R/ora_analysis.R` depend on both `clusterProfiler`
and `org.Hs.eg.db`. Make sure those packages are installed and available when
running any targets that trigger enrichment; otherwise the pipeline will stop
early with a clear error message.

## Running tests

Basic automated checks live under `tests/testthat`. To execute them:

```r
Rscript -e 'testthat::test_dir("tests/testthat")'
```
