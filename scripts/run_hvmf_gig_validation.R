args <- commandArgs(trailingOnly = TRUE)

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

parse_integer_arg <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop("Arguments that must be positive integers were not valid.")
  }

  parsed
}

parse_numeric_arg <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    return(default)
  }

  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) {
    stop("Numeric arguments must be finite.")
  }

  parsed
}

output_dir <- if (length(args) >= 1L) args[[1]] else file.path("bootstrap", "results")
mle_replicates <- parse_integer_arg(if (length(args) >= 2L) args[[2]] else NULL, 1000L)
gof_outer <- parse_integer_arg(if (length(args) >= 3L) args[[3]] else NULL, 1000L)
gof_inner <- parse_integer_arg(if (length(args) >= 4L) args[[4]] else NULL, 1000L)
gof_n_cores <- parse_integer_arg(if (length(args) >= 5L) args[[5]] else NULL, 3L)
gof_n <- parse_integer_arg(if (length(args) >= 6L) args[[6]] else NULL, 100L)
gof_kappa <- parse_numeric_arg(if (length(args) >= 7L) args[[7]] else NULL, 1.5)

source(resolve_script_path("bootstrap", "hvmf_gig_study.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("MLE replicates: %d\n", mle_replicates))
cat(sprintf("GOF outer replicates: %d\n", gof_outer))
cat(sprintf("GOF inner bootstrap: %d\n", gof_inner))
cat(sprintf("GOF n_cores: %d\n", gof_n_cores))

mle_csv <- file.path(output_dir, "hvmf_gig_mle_monte_carlo_summary.csv")
mle_summary <- run_hvmf_gig_mle_monte_carlo_study(
  output_csv = mle_csv,
  n_replicates = mle_replicates
)
print(mle_summary)

gof_result <- run_hvmf_gig_simple_cvm_calibration_study(
  output_dir = output_dir,
  scenarios = list(
    list(n = gof_n, kappa = gof_kappa, M = gof_outer, B = gof_inner, n_cores = gof_n_cores)
  )
)

cat(sprintf("GOF summary CSV: %s\n", gof_result$summary_csv))
cat(sprintf("GOF log: %s\n", gof_result$log_txt))
cat(sprintf("GOF RDS: %s\n", gof_result$summary$output_rds[[1]]))