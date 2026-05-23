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
    return("smoke")
  }

  value <- tolower(value)
  if (!value %in% c("smoke", "full", "smoke_composite", "full_composite")) {
    stop("Mode must be one of 'smoke', 'full', 'smoke_composite', 'full_composite'.")
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
output_dir <- if (length(args) >= 4L) args[[4]] else NULL

mode <- parse_mode(mode_arg)
n_cores_outer <- parse_integer_arg(n_cores_arg, "n_cores_outer", 10L)
seed <- parse_integer_arg(seed_arg, "seed", 123L)
if (!is.null(output_dir) && identical(output_dir, "")) {
  output_dir <- NULL
}

source(resolve_script_path("bootstrap", "calibration_study.R"))

cat(sprintf("Running multiplier bootstrap calibration study in '%s' mode.\n", mode))
cat(sprintf("n_cores_outer = %d\n", n_cores_outer))
cat(sprintf("seed = %d\n", seed))
if (!is.null(output_dir)) {
  cat(sprintf("output_dir = %s\n", output_dir))
}

result <- if (identical(mode, "smoke")) {
  run_smoke_bootstrap_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed
  )
} else if (identical(mode, "smoke_composite")) {
  run_smoke_bootstrap_composite_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed
  )
} else if (identical(mode, "full_composite")) {
  run_full_bootstrap_composite_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed
  )
} else {
  run_full_bootstrap_calibration_study(
    output_dir = output_dir,
    n_cores_outer = n_cores_outer,
    seed = seed
  )
}

cat("Calibration study completed.\n")
cat(sprintf("Raw CSV: %s\n", result$raw_csv))
cat(sprintf("Summary CSV: %s\n", result$summary_csv))
