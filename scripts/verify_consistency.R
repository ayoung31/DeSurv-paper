#!/usr/bin/env Rscript
# DeSurv Consistency Verification Script
# Run: Rscript scripts/verify_consistency.R
#
# Checks:
# 1. Store path consistency across all config files
# 2. Config file validation
# 3. Pipeline validation
# 4. Debug statement detection
# 5. Figure target status

suppressPackageStartupMessages({
  library(targets)
  library(yaml)
})

cat("
╔══════════════════════════════════════════════════════════════╗
║           DeSurv Consistency Verification                     ║
╚══════════════════════════════════════════════════════════════╝
")

errors <- 0
warnings <- 0

# Helper function
check_result <- function(passed, message, is_warning = FALSE) {
  if (passed) {
    cat(sprintf("   ✓ %s\n", message))
  } else if (is_warning) {
    cat(sprintf("   ⚠ %s\n", message))
    warnings <<- warnings + 1
  } else {
    cat(sprintf("   ✗ %s\n", message))
    errors <<- errors + 1
  }
  invisible(passed)
}

# =============================================================================
# 1. Store Consistency
# =============================================================================
cat("\n1. Store Path Consistency\n")
cat("   ─────────────────────\n")

store_refs <- list()

# Check paper/_targets.yaml
tryCatch({
  yaml_content <- yaml::read_yaml("paper/_targets.yaml")
  store_refs$yaml <- yaml_content$main$store
  check_result(TRUE, sprintf("paper/_targets.yaml: %s", store_refs$yaml))
}, error = function(e) {
  check_result(FALSE, sprintf("paper/_targets.yaml: %s", e$message))
})

# Check paper/paper.Rmd params
tryCatch({
  rmd_yaml <- rmarkdown::yaml_front_matter("paper/paper.Rmd")
  store_refs$rmd <- rmd_yaml$params$tar_store
  check_result(TRUE, sprintf("paper/paper.Rmd params: %s", store_refs$rmd))
}, error = function(e) {
  check_result(FALSE, sprintf("paper/paper.Rmd: %s", e$message))
})

# Check _targets_sims_local.sh
tryCatch({
  sims_content <- readLines("_targets_sims_local.sh", warn = FALSE)
  store_line <- grep("tar_config_set.*store", sims_content, value = TRUE)
  if (length(store_line) > 0) {
    store_match <- regmatches(store_line, regexpr("store_[^\"']+", store_line))
    store_refs$sims <- store_match
    check_result(TRUE, sprintf("_targets_sims_local.sh: %s", store_refs$sims))
  }
}, error = function(e) {
  check_result(FALSE, sprintf("_targets_sims_local.sh: %s", e$message), is_warning = TRUE)
})

# Verify all match
unique_stores <- unique(unlist(store_refs))
if (length(unique_stores) == 1) {
  check_result(TRUE, sprintf("All stores match: %s", unique_stores))
} else {
  check_result(FALSE, sprintf("Store MISMATCH: %s", paste(unique_stores, collapse = " vs ")))
}

# =============================================================================
# 2. Config Validation
# =============================================================================
cat("\n2. Configuration Validation\n")
cat("   ────────────────────────\n")

# Check BO configs
tryCatch({
  source("targets_bo_configs.R", local = TRUE)
  bo <- targets_bo_configs()
  check_result(TRUE, sprintf("BO configs valid (%d configs)", length(bo)))

  # Check required fields
  required <- c("desurv_bo_bounds", "ngene_config", "ntop_config", "bo_n_init", "bo_n_iter")
  for (config_name in names(bo)) {
    missing <- setdiff(required, names(bo[[config_name]]))
    if (length(missing) > 0) {
      check_result(FALSE, sprintf("%s missing: %s", config_name, paste(missing, collapse = ", ")))
    }
  }
}, error = function(e) {
  check_result(FALSE, sprintf("BO config error: %s", e$message))
})

# Check run configs
tryCatch({
  source("targets_run_configs.R", local = TRUE)
  run <- targets_run_configs()
  check_result(TRUE, sprintf("Run configs valid (%d configs)", length(run)))
}, error = function(e) {
  check_result(FALSE, sprintf("Run config error: %s", e$message))
})

# Check val configs
tryCatch({
  source("targets_val_configs.R", local = TRUE)
  val <- targets_val_configs()
  check_result(TRUE, sprintf("Val configs valid (%d configs)", length(val)))
}, error = function(e) {
  check_result(FALSE, sprintf("Val config error: %s", e$message))
})

# =============================================================================
# 3. Pipeline Validation
# =============================================================================
cat("\n3. Pipeline Validation\n")
cat("   ───────────────────\n")

# Main pipeline
tryCatch({
  tar_validate(script = "_targets.R")
  check_result(TRUE, "Main pipeline (_targets.R) valid")
}, error = function(e) {
  check_result(FALSE, sprintf("Main pipeline error: %s", e$message))
})

# Simulation pipeline
tryCatch({
  tar_validate(script = "_targets_sims.R")
  check_result(TRUE, "Simulation pipeline (_targets_sims.R) valid")
}, error = function(e) {
  check_result(FALSE, sprintf("Simulation pipeline error: %s", e$message))
})

# =============================================================================
# 4. Debug Statement Check
# =============================================================================
cat("\n4. Debug Statement Check\n")
cat("   ─────────────────────\n")

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
r_files <- c(r_files, list.files(".", pattern = "^_targets.*\\.R$", full.names = TRUE))

debug_found <- character()
for (f in r_files) {
  content <- readLines(f, warn = FALSE)
  # Skip commented lines
  active_lines <- content[!grepl("^\\s*#", content)]
  if (any(grepl("browser\\(\\)", active_lines))) {
    debug_found <- c(debug_found, f)
  }
}

if (length(debug_found) == 0) {
  check_result(TRUE, "No active browser() statements found")
} else {
  for (f in debug_found) {
    check_result(FALSE, sprintf("browser() found in: %s", f))
  }
}

# =============================================================================
# 5. Figure Target Status
# =============================================================================
cat("\n5. Figure Target Status\n")
cat("   ────────────────────\n")

if (length(store_refs) > 0 && !is.null(store_refs[[1]])) {
  tryCatch({
    tar_config_set(store = store_refs[[1]])
    meta <- tar_meta()

    fig_targets <- meta[grepl("^fig_", meta$name), ]
    check_result(TRUE, sprintf("Found %d figure targets", nrow(fig_targets)))

    # Check for errors
    if ("error" %in% names(fig_targets)) {
      errored <- fig_targets[!is.na(fig_targets$error) & nzchar(fig_targets$error), ]
      if (nrow(errored) > 0) {
        for (i in seq_len(nrow(errored))) {
          check_result(FALSE, sprintf("Errored: %s", errored$name[i]))
        }
      } else {
        check_result(TRUE, "No errored figure targets")
      }
    }

    # Check sim targets
    sim_targets <- meta[grepl("^sim_", meta$name), ]
    if (nrow(sim_targets) > 0) {
      check_result(TRUE, sprintf("Found %d simulation targets", nrow(sim_targets)))
    } else {
      check_result(FALSE, "No simulation targets found", is_warning = TRUE)
    }

  }, error = function(e) {
    check_result(FALSE, sprintf("Cannot check targets: %s", e$message), is_warning = TRUE)
  })
} else {
  check_result(FALSE, "Cannot determine store path", is_warning = TRUE)
}

# =============================================================================
# 6. Documentation Files
# =============================================================================
cat("\n6. Documentation Check\n")
cat("   ───────────────────\n")

required_docs <- c(
  "CLAUDE.md",
  "CHANGELOG.md",
  "DEFAULTS.md",
  "CONSISTENCY_STANDARDS.md",
  "MANUSCRIPT_PIPELINE_ANALYSIS.md"
)

for (doc in required_docs) {
  if (file.exists(doc)) {
    check_result(TRUE, sprintf("%s exists", doc))
  } else {
    check_result(FALSE, sprintf("%s missing", doc), is_warning = TRUE)
  }
}

# =============================================================================
# Summary
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  Summary: %d errors, %d warnings                              ║\n", errors, warnings))
cat("╚══════════════════════════════════════════════════════════════╝\n")

if (errors > 0) {
  cat("\n⚠ Fix errors before committing!\n")
  quit(status = 1)
} else if (warnings > 0) {
  cat("\n✓ Passed with warnings. Review before pushing.\n")
} else {
  cat("\n✓ All checks passed!\n")
}
