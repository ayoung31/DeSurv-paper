if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required to run the test suite.", call. = FALSE)
}

test_files <- list.files("R", full.names = TRUE, pattern = "[.]R$")
invisible(lapply(test_files, source))

testthat::test_dir("tests/testthat", reporter = "summary")
