args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_quadform_script_path <- function(...) {
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

parse_bool_arg <- function(x, default) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  tolower(x) %in% c("true", "1", "yes", "y")
}

source(resolve_quadform_script_path("tests", "benchmark_logistic_gaussian_quadform_utils.R"))

named <- parse_named_args(args)
output_dir <- named$output_dir %||% file.path("tests", "benchmark_outputs", "logistic_gaussian_quadform")
n_cores <- min(12L, parse_integer_arg(named$n_cores, 12L))
random_n <- parse_integer_arg(named$random_n, 240L)
random_seed <- parse_integer_arg(named$random_seed, 123L)
save_plots <- parse_bool_arg(named$save_plots, TRUE)
max_cases <- if (is.null(named$max_cases)) NULL else parse_integer_arg(named$max_cases, 1L)
stress <- parse_bool_arg(named$stress, TRUE)
realistic <- parse_bool_arg(named$realistic, TRUE)
random_cases <- parse_bool_arg(named$random_cases, TRUE)

cat("Running logistic Gaussian quadratic-form benchmark\n")
cat(sprintf("output_dir = %s\n", output_dir))
cat(sprintf("n_cores = %d\n", n_cores))
cat(sprintf("random_n = %d\n", random_n))
cat(sprintf("random_seed = %d\n", random_seed))
cat(sprintf("max_cases = %s\n", if (is.null(max_cases)) "NULL" else as.character(max_cases)))

result <- run_logistic_gaussian_quadform_benchmark(
  output_dir = output_dir,
  n_cores = n_cores,
  random_n = random_n,
  random_seed = random_seed,
  max_cases = max_cases,
  stress = stress,
  realistic = realistic,
  random_cases = random_cases,
  save_plots = save_plots
)

cat(sprintf("Saved benchmark to %s\n", result$output_dir))
