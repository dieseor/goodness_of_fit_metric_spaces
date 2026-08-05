#!/usr/bin/env Rscript

# Monte Carlo sensitivity analysis for the BHEP normality test on the SkyeAFM
# ilr data. This script deliberately uses one CPU core and does not rerun KS/CvM.
#
# It writes two CSV files:
#   1. bhep_skyeafm_mc_sensitivity.csv: p-values at nested MC sizes for
#      a = 0.5, 1, 2.
#   2. bhep_skyeafm_a_sensitivity.csv: p-values at MC = 10000 for a values
#      from 0.5 to 20.
#
# The same null simulations are reused. Thus smaller MC sizes are prefixes of
# the larger run, and the a-sensitivity uses its first 10000 replications.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

parse_csv_numeric <- function(value, name) {
  result <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(result) || any(!is.finite(result)) || any(result <= 0)) {
    stop(sprintf("`%s` must be a comma-separated list of positive numbers.", name))
  }
  result
}

parse_arguments <- function(args) {
  defaults <- list(
    results_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_dir = NULL,
    max_mc = 50000L,
    mc_prefixes = "1000,2500,5000,10000,25000,50000",
    mc_for_a = 10000L,
    a_base = "0.5,1,2",
    a_values = "0.5,1,2,3,5,7.5,10,15,20",
    alpha = 0.05,
    seed = 20260803L
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) stop("Arguments must have the form --name=value.")
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L || !parts[[1L]] %in% names(defaults)) {
      stop(sprintf("Unknown or invalid argument: %s", arg))
    }
    defaults[[parts[[1L]]]] <- parts[[2L]]
  }

  defaults$max_mc <- as.integer(defaults$max_mc)
  defaults$mc_for_a <- as.integer(defaults$mc_for_a)
  defaults$alpha <- as.numeric(defaults$alpha)
  defaults$seed <- as.integer(defaults$seed)
  defaults$mc_prefixes <- as.integer(parse_csv_numeric(defaults$mc_prefixes, "mc_prefixes"))
  defaults$a_base <- parse_csv_numeric(defaults$a_base, "a_base")
  defaults$a_values <- parse_csv_numeric(defaults$a_values, "a_values")

  if (is.null(defaults$output_dir)) {
    defaults$output_dir <- file.path(defaults$results_dir, "bhep_skyeafm_sensitivity")
  }
  if (is.na(defaults$max_mc) || defaults$max_mc < 1L ||
      is.na(defaults$mc_for_a) || defaults$mc_for_a < 1L ||
      any(defaults$mc_prefixes > defaults$max_mc) ||
      defaults$mc_for_a > defaults$max_mc ||
      is.na(defaults$alpha) || defaults$alpha <= 0 || defaults$alpha >= 1 ||
      is.na(defaults$seed) || defaults$seed < 0L) {
    stop("Incompatible Monte Carlo settings.")
  }
  if (!all(defaults$a_base %in% defaults$a_values)) {
    stop("Every value in `a_base` must also occur in `a_values`.")
  }
  defaults
}

load_skyeafm_ilr <- function(results_dir) {
  path <- file.path(results_dir, "skyeafm_results.rds")
  if (!file.exists(path)) stop(sprintf("Missing input: %s", path))
  result <- readRDS(path)
  z <- as.matrix(result$fit$Z)
  storage.mode(z) <- "double"
  if (!is.numeric(z) || any(!is.finite(z)) || nrow(z) != 23L || ncol(z) != 2L) {
    stop("SkyeAFM ilr data must be a finite 23 by 2 matrix.")
  }
  if (qr(sweep(z, 2L, colMeans(z), FUN = "-"))$rank != ncol(z)) {
    stop("SkyeAFM ilr sample covariance is singular.")
  }
  z
}

mc_p_value <- function(exceedances, mc_rep) (1L + exceedances) / (mc_rep + 1L)

main <- function() {
  if (!requireNamespace("mnt", quietly = TRUE)) {
    stop("Package `mnt` is required.")
  }
  settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
  z <- load_skyeafm_ilr(settings$results_dir)
  n <- nrow(z)
  d <- ncol(z)
  observed <- vapply(settings$a_values, function(a) as.numeric(mnt::BHEP(z, a = a)), numeric(1))
  names(observed) <- as.character(settings$a_values)
  exceedances <- setNames(integer(length(settings$a_values)), names(observed))
  a_exceedances <- setNames(integer(length(settings$a_values)), names(observed))
  prefix_exceedances <- matrix(
    0L,
    nrow = length(settings$mc_prefixes),
    ncol = length(settings$a_base),
    dimnames = list(as.character(settings$mc_prefixes), as.character(settings$a_base))
  )

  set.seed(settings$seed)
  for (replication in seq_len(settings$max_mc)) {
    null_sample <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    a_now <- if (replication <= settings$mc_for_a) settings$a_values else settings$a_base
    null_statistics <- vapply(a_now, function(a) as.numeric(mnt::BHEP(null_sample, a = a)), numeric(1))
    exceedances[as.character(a_now)] <- exceedances[as.character(a_now)] +
      as.integer(null_statistics >= observed[as.character(a_now)])

    if (replication %in% settings$mc_prefixes) {
      prefix_exceedances[as.character(replication), ] <- exceedances[as.character(settings$a_base)]
    }
    if (replication == settings$mc_for_a) {
      a_exceedances <- exceedances
    }
    if (replication %% 1000L == 0L || replication == settings$max_mc) {
      message(sprintf("Completed %d / %d null replications.", replication, settings$max_mc))
    }
  }

  mc_output <- do.call(rbind, lapply(seq_along(settings$a_base), function(j) {
    mc_rep <- settings$mc_prefixes
    count <- prefix_exceedances[, j]
    p_value <- mc_p_value(count, mc_rep)
    data.frame(
      dataset = "SkyeAFM", n = n, D = d + 1L, ilr_dimension = d,
      a = settings$a_base[[j]], mc_rep = mc_rep, test_value = observed[[as.character(settings$a_base[[j]])]],
      mc_exceedances = as.integer(count), p_value = p_value,
      mc_standard_error = sqrt(p_value * (1 - p_value) / mc_rep),
      reject = p_value <= settings$alpha, alpha = settings$alpha,
      seed = settings$seed, mnt_version = as.character(utils::packageVersion("mnt")),
      stringsAsFactors = FALSE
    )
  }))

  a_count <- a_exceedances[as.character(settings$a_values)]
  a_p_value <- mc_p_value(a_count, settings$mc_for_a)
  a_output <- data.frame(
    dataset = "SkyeAFM", n = n, D = d + 1L, ilr_dimension = d,
    a = settings$a_values, mc_rep = settings$mc_for_a, test_value = observed,
    mc_exceedances = as.integer(a_count), p_value = a_p_value,
    mc_standard_error = sqrt(a_p_value * (1 - a_p_value) / settings$mc_for_a),
    reject = a_p_value <= settings$alpha, alpha = settings$alpha,
    seed = settings$seed, mnt_version = as.character(utils::packageVersion("mnt")),
    stringsAsFactors = FALSE
  )

  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(mc_output, file.path(settings$output_dir, "bhep_skyeafm_mc_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(a_output, file.path(settings$output_dir, "bhep_skyeafm_a_sensitivity.csv"), row.names = FALSE)
  message(sprintf("Wrote sensitivity outputs to %s", settings$output_dir))
}

main()
