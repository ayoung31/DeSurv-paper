test_that("resolve_desurv_run_configs_by_bo aligns labels", {
  bo_configs <- resolve_desurv_bo_configs(list(a = list(), b = list()))
  run_raw <- list(
    a = list(ninit_full = 5),
    b = list(ninit_full = 6)
  )
  run_configs <- resolve_desurv_run_configs_by_bo(bo_configs, run_raw)
  expect_named(run_configs, names(bo_configs))
  expect_equal(
    vapply(run_configs, function(cfg) cfg$label, character(1)),
    names(bo_configs)
  )
})

test_that("resolve_desurv_run_configs_by_bo errors on missing labels", {
  bo_configs <- resolve_desurv_bo_configs(list(a = list(), b = list()))
  run_raw <- list(a = list())
  expect_error(
    resolve_desurv_run_configs_by_bo(bo_configs, run_raw),
    "Missing run configs"
  )
})

test_that("resolve_desurv_val_configs_by_bo aligns labels", {
  bo_configs <- resolve_desurv_bo_configs(list(a = list(), b = list()))
  val_raw <- list(
    a = list(mode = "external", val_datasets = "CPTAC"),
    b = list(mode = "external", val_datasets = "CPTAC")
  )
  val_configs <- resolve_desurv_val_configs_by_bo(bo_configs, val_raw)
  expect_named(val_configs, names(bo_configs))
  expect_equal(
    vapply(val_configs, function(cfg) cfg$label, character(1)),
    names(bo_configs)
  )
})

test_that("validate_desurv_configs enforces train_split mode", {
  bo_configs <- resolve_desurv_bo_configs(list(a = list()))
  run_configs <- resolve_desurv_run_configs_by_bo(bo_configs, list(a = list()))
  val_configs <- resolve_desurv_val_configs_by_bo(
    bo_configs,
    list(a = list(mode = "train_split", val_datasets = character(0)))
  )
  expect_error(
    validate_desurv_configs(bo_configs, run_configs, val_configs),
    "train_split"
  )
})
