testthat::local_edition(3)

project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
r_dir <- file.path(project_root, "R")
r_files <- list.files(r_dir, full.names = TRUE, pattern = "[.]R$")
for (path in r_files) {
  sys.source(path, envir = topenv())
}
