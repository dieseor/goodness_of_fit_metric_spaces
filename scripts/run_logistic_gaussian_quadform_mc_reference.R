args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_quadform_mc_script_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) >= 2L) paste(pieces[-1L], collapse = "=") else "true"
    out[[key]] <- value
  }
  out
}

parse_integer_arg <- function(x, default) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  value <- suppressWarnings(as.integer(x))
  if (is.na(value) || value < 1L) {
    stop("Expected a positive integer argument.")
  }
  value
}

`%||%` <- function(x, y) if (is.null(x)) y else x

source(resolve_quadform_mc_script_path("tests", "benchmark_logistic_gaussian_quadform_utils.R"))

named <- parse_named_args(args)
benchmark_results_csv <- named$benchmark_results_csv %||%
  file.path("tests", "benchmark_outputs", "logistic_gaussian_quadform", "20260610_full_12cores", "quadform_results.csv")
output_dir <- named$output_dir %||%
  file.path("tests", "benchmark_outputs", "logistic_gaussian_quadform", "20260610_mc_reference_12cores")
n_cores <- min(12L, parse_integer_arg(named$n_cores, 12L))
M <- parse_integer_arg(named$M, 2000000L)
chunk_size <- parse_integer_arg(named$chunk_size, 200000L)
seed <- parse_integer_arg(named$seed, 123L)

cat("Running logistic Gaussian quadratic-form Monte Carlo reference study\n")
cat(sprintf("benchmark_results_csv = %s\n", benchmark_results_csv))
cat(sprintf("output_dir = %s\n", output_dir))
cat(sprintf("n_cores = %d\n", n_cores))
cat(sprintf("M = %d\n", M))
cat(sprintf("chunk_size = %d\n", chunk_size))
cat(sprintf("seed = %d\n", seed))

result <- run_quadform_mc_reference_study(
  output_dir = output_dir,
  benchmark_results_csv = benchmark_results_csv,
  n_cores = n_cores,
  M = M,
  chunk_size = chunk_size,
  seed = seed
)

cat(sprintf("Saved MC reference study to %s\n", result$output_dir))
