#!/usr/bin/env Rscript

## Canonical repository test launcher.
## Run serially: several tests use shared ports, temporary files, compiled
## backends, or process-global state and are not safe to execute concurrently.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[startsWith(all_args, "--file=")]
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "scripts", "run_all_tests.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(repo_root)

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The `testthat` package is required to run the repository tests.")
}

cat(sprintf("Repository: %s\n", repo_root))
cat(sprintf("Test directory: %s\n", file.path(repo_root, "tests", "testthat")))
cat(sprintf("Started: %s\n", format(Sys.time())))
flush.console()

started <- proc.time()[["elapsed"]]
testthat::test_dir(
  file.path(repo_root, "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE,
  load_helpers = TRUE
)
elapsed <- proc.time()[["elapsed"]] - started
cat(sprintf("\nAll testthat tests passed. Elapsed: %.1f seconds.\n", elapsed))
