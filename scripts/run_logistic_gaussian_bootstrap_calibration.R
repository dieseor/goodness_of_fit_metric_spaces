args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_script_path <- function(...) {
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

parse_mode <- function(value) {
  if (is.null(value) || identical(value, "")) {
    return("both")
  }

  value <- tolower(value)
  if (!value %in% c("simple", "composite", "both", "smoke_simple", "smoke_composite")) {
    stop("Mode must be one of 'simple', 'composite', 'both', 'smoke_simple', 'smoke_composite'.")
  }

  value
}

parse_integer_arg <- function(value, name, default) {
  if (is.null(value) || identical(value, "")) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(sprintf("%s must be a positive integer.", name))
  }

  parsed
}

mode_arg <- if (length(args) >= 1L) args[[1]] else NULL
n_cores_arg <- if (length(args) >= 2L) args[[2]] else NULL
seed_arg <- if (length(args) >= 3L) args[[3]] else NULL
output_dir_arg <- if (length(args) >= 4L) args[[4]] else NULL

mode <- parse_mode(mode_arg)
n_cores_outer <- parse_integer_arg(n_cores_arg, "n_cores_outer", 8L)
seed <- parse_integer_arg(seed_arg, "seed", 123L)
output_dir_arg <- if (!is.null(output_dir_arg) && identical(output_dir_arg, "")) NULL else output_dir_arg

source(resolve_script_path("bootstrap", "calibration_study.R"))

cat(sprintf("Running logistic Gaussian calibration mode '%s'.\n", mode))
cat(sprintf("n_cores_outer = %d\n", n_cores_outer))
cat(sprintf("seed = %d\n", seed))

if (identical(mode, "smoke_simple")) {
  result <- run_smoke_logistic_gaussian_calibration_study(
    output_dir = output_dir_arg,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = TRUE,
    verbose = TRUE
  )
  cat(sprintf("Raw CSV: %s\n", result$raw_csv))
  cat(sprintf("Summary CSV: %s\n", result$summary_csv))
  quit(save = "no", status = 0)
}

if (identical(mode, "smoke_composite")) {
  result <- run_smoke_logistic_gaussian_composite_calibration_study(
    output_dir = output_dir_arg,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = TRUE,
    verbose = TRUE
  )
  cat(sprintf("Raw CSV: %s\n", result$raw_csv))
  cat(sprintf("Summary CSV: %s\n", result$summary_csv))
  quit(save = "no", status = 0)
}

if (identical(mode, "simple") || identical(mode, "both")) {
  output_dir_simple <- if (identical(mode, "simple")) {
    output_dir_arg
  } else if (is.null(output_dir_arg)) {
    NULL
  } else {
    file.path(output_dir_arg, "simple")
  }

  result_simple <- run_full_logistic_gaussian_calibration_study(
    output_dir = output_dir_simple,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = TRUE,
    verbose = TRUE
  )

  cat("Simple calibration completed.\n")
  cat(sprintf("Raw CSV: %s\n", result_simple$raw_csv))
  cat(sprintf("Summary CSV: %s\n", result_simple$summary_csv))
}

if (identical(mode, "composite") || identical(mode, "both")) {
  output_dir_composite <- if (identical(mode, "composite")) {
    output_dir_arg
  } else if (is.null(output_dir_arg)) {
    NULL
  } else {
    file.path(output_dir_arg, "composite")
  }

  result_composite <- run_full_logistic_gaussian_composite_calibration_study(
    output_dir = output_dir_composite,
    n_cores_outer = n_cores_outer,
    seed = seed,
    show_progress = TRUE,
    verbose = TRUE
  )

  cat("Composite calibration completed.\n")
  cat(sprintf("Raw CSV: %s\n", result_composite$raw_csv))
  cat(sprintf("Summary CSV: %s\n", result_composite$summary_csv))
}
