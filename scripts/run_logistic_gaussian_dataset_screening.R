#!/usr/bin/env Rscript

# This script runs the logistic Gaussian goodness-of-fit analysis for real
# compositional datasets. By default, it uses the parametric composite-null
# bootstrap: parameters are estimated from the observed data, bootstrap samples
# are generated from the fitted logistic Gaussian model, and parameters are
# re-estimated in every bootstrap sample before computing the bootstrap
# statistics.
#
# The legacy plug-in simple-null mode is still available through
# --bootstrap_mode=plugin_simple_null, but it is intended only for exploratory
# checks and should not be interpreted as a valid composite-null p-value.

resolve_logistic_gaussian_screening_runner_path <- function(...) {
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

utils_path_logistic_gaussian_screening <- resolve_logistic_gaussian_screening_runner_path(
  "real_data",
  "logistic_gaussian",
  "utils_logistic_gaussian_screening.R"
)
source(utils_path_logistic_gaussian_screening)

parse_screening_args <- function(args) {
  output <- list()

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[[1]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[key]] <- value
  }

  output
}

parse_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

parse_integer_arg <- function(x, default) {
  if (is.null(x)) {
    return(as.integer(default))
  }
  as.integer(x)
}

parse_numeric_arg <- function(x, default) {
  if (is.null(x)) {
    return(as.numeric(default))
  }
  as.numeric(x)
}

parse_optional_numeric_arg <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.numeric(x)
}

run_logistic_gaussian_dataset_screening_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- parse_screening_args(args)

  dataset_names <- if (is.null(args$datasets)) {
    default_logistic_gaussian_screening_datasets()
  } else {
    trimws(strsplit(args$datasets, ",", fixed = TRUE)[[1]])
  }

  B <- parse_integer_arg(args$B, 1000L)
  max_centers <- parse_integer_arg(args$max_centers, 100L)
  n_t <- parse_integer_arg(args$n_t, 60L)
  t_grid_tail_prob <- parse_numeric_arg(args$t_grid_tail_prob, 1e-8)
  boundary_epsilon <- parse_optional_numeric_arg(args$boundary_epsilon)
  seed <- parse_integer_arg(args$seed, 123L)
  alpha <- parse_numeric_arg(args$alpha, 0.05)
  ridge <- parse_numeric_arg(args$ridge, 1e-8)
  bootstrap_mode <- as.character(args$bootstrap_mode %||% "composite_multiplier")
  n_cores <- parse_integer_arg(args$n_cores, 1L)
  omega_grid_type <- as.character(args$omega_grid_type %||% "fixed_simplex_lattice")
  t_grid_type <- as.character(args$t_grid_type %||% "fixed_fitted_scale")
  quadform_method <- as.character(args$quadform_method %||% "hbe")
  output_dir <- as.character(args$output_dir %||% file.path("output", "calibration", "bootstrap", "logistic_gaussian", "composite"))
  make_plots <- parse_flag(args$make_plots, default = TRUE)
  run_seed_sensitivity <- parse_flag(args$run_seed_sensitivity, default = FALSE)

  if (!bootstrap_mode %in% c("composite_multiplier", "composite_parametric", "plugin_simple_null")) {
    stop(sprintf(
      "Unsupported bootstrap_mode '%s'. Use 'composite_multiplier', 'composite_parametric' or 'plugin_simple_null'.",
      bootstrap_mode
    ))
  }

  cat(logistic_gaussian_screening_warning, "\n\n", sep = "")
  cat("Datasets:", paste(dataset_names, collapse = ", "), "\n")
  cat("Number of datasets:", length(dataset_names), "\n")
  cat("B:", B, "\n")
  cat("max_centers:", max_centers, "\n")
  cat("n_t:", n_t, "\n")
  cat("t_grid_tail_prob:", t_grid_tail_prob, "\n")
  cat("boundary_epsilon:", if (is.null(boundary_epsilon)) "default(0.5/D)" else boundary_epsilon, "\n")
  cat("seed:", seed, "\n")
  cat("bootstrap_mode:", bootstrap_mode, "\n")
  cat("n_cores:", n_cores, "\n")
  cat("omega_grid_type:", omega_grid_type, "\n")
  cat("t_grid_type:", t_grid_type, "\n")
  cat("quadform_method:", quadform_method, "\n")
  cat("output_dir:", output_dir, "\n\n")
  cat("Composite bootstrap mode is the default. Use --bootstrap_mode=plugin_simple_null only for legacy exploratory checks.\n\n")

  batch_result <- run_logistic_gaussian_screening_batch(
    dataset_names = dataset_names,
    B = B,
    max_centers = max_centers,
    n_t = n_t,
    t_grid_tail_prob = t_grid_tail_prob,
    boundary_epsilon = boundary_epsilon,
    bootstrap_mode = bootstrap_mode,
    seed = seed,
    alpha = alpha,
    ridge = ridge,
    n_cores = n_cores,
    control = list(logistic_gaussian_quadform_method = quadform_method),
    omega_grid_type = omega_grid_type,
    t_grid_type = t_grid_type,
    make_plots = make_plots,
    output_dir = output_dir,
    run_seed_sensitivity = run_seed_sensitivity,
    verbose = TRUE
  )

  cat("Summary CSV:", batch_result$summary_csv, "\n")
  print(batch_result$summary)
  invisible(batch_result)
}

if (sys.nframe() == 0L) {
  run_logistic_gaussian_dataset_screening_cli()
}
