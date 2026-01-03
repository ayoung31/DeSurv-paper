
can_load_package <- function(pkg) {
  tryCatch({
    requireNamespace(pkg, quietly = TRUE)
  }, error = function(e) {
    FALSE
  })
}

skip_if_packages_unavailable <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, can_load_package, logical(1))]
  if (length(missing)) {
    skip(sprintf("Required packages not available: %s", paste(missing, collapse = ", ")))
  }
}

skip_if_pipeline_data_missing <- function() {
  path <- file.path("data", "derv", "cmbSubtypes_formatted.RData")
  if (!file.exists(path)) {
    skip(sprintf("Required pipeline data missing at %s", path))
  }

}

test_that("targets manifest loads and includes core pipeline targets", {
  skip_if_packages_unavailable(c("targets", "tarchetypes", "crew", "gert"))
  skip_if_pipeline_data_missing()

  manifest <- local({
    old <- Sys.getenv("DESURV_LOCAL_RENDER", unset = NA_character_)
    on.exit({
      if (is.na(old)) {
        Sys.unsetenv("DESURV_LOCAL_RENDER")
      } else {
        Sys.setenv(DESURV_LOCAL_RENDER = old)
      }
    }, add = TRUE)
    Sys.setenv(DESURV_LOCAL_RENDER = "1")
    targets::tar_manifest(script = "_targets.R")
  })

  target_names <- manifest$name
  expect_equal(length(target_names), length(unique(target_names)))
  expect_true(all(c(
    "bo_config",
    "run_config",
    "val_config",
    "bo_bundle",
    "bo_bundles",
    "run_bundle",
    "run_bundles",
    "val_run_bundle",
    "data",
    "data_val",
    "data_val_filtered",
    "desurv_bo_results"
  ) %in% target_names))
})

test_that("targets manifest captures split/external validation behavior", {
  skip_if_packages_unavailable(c("targets", "tarchetypes", "crew", "gert"))
  skip_if_pipeline_data_missing()

  manifest <- local({
    old <- Sys.getenv("DESURV_LOCAL_RENDER", unset = NA_character_)
    on.exit({
      if (is.na(old)) {
        Sys.unsetenv("DESURV_LOCAL_RENDER")
      } else {
        Sys.setenv(DESURV_LOCAL_RENDER = old)
      }
    }, add = TRUE)
    Sys.setenv(DESURV_LOCAL_RENDER = "1")
    targets::tar_manifest(script = "_targets.R")
  })

  data_val_cmd <- manifest$command[manifest$name == "data_val"]
  expect_true(any(grepl("data_split\\$test", data_val_cmd)))
  expect_true(any(grepl("val_config\\$mode", data_val_cmd)))
})

test_that("targets manifest includes tagged results paths", {
  skip_if_packages_unavailable(c("targets", "tarchetypes", "crew", "gert"))
  skip_if_pipeline_data_missing()

  manifest <- local({
    old <- Sys.getenv("DESURV_LOCAL_RENDER", unset = NA_character_)
    on.exit({
      if (is.na(old)) {
        Sys.unsetenv("DESURV_LOCAL_RENDER")
      } else {
        Sys.setenv(DESURV_LOCAL_RENDER = old)
      }
    }, add = TRUE)
    Sys.setenv(DESURV_LOCAL_RENDER = "1")
    targets::tar_manifest(script = "_targets.R")
  })

  training_cmd <- manifest$command[manifest$name == "training_results_dir"]
  bo_cmd <- manifest$command[manifest$name == "bo_results_dir"]
  expect_true(any(grepl("run_config\\$path_tag", training_cmd)))
  expect_true(any(grepl("bo_config\\$path_tag", bo_cmd)))
})
