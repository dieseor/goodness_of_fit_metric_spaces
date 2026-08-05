#!/usr/bin/env Rscript

# Reproducible Monte Carlo BHEP p-values for the 29 compositional datasets in
# AoS/gof_metric_spaces_aos.tex. It reuses the fitted ilr matrices saved by the
# paper-table run, and deliberately does not rerun the KS/CvM bootstrap.
#
# mnt::test.BHEP() returns a statistic, a Monte Carlo critical value, and a
# decision, but no p-value. This runner calibrates the same mnt::BHEP()
# statistic under the normal null and reports the standard plus-one Monte Carlo
# p-value (1 + #{T_b >= T_obs}) / (MC.rep + 1).

# Keep the Monte Carlo run on one CPU core, including numerical backends that
# might otherwise create worker threads for matrix operations.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_arguments <- function(args) {
  defaults <- list(
    results_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_dir = NULL,
    mc_rep = 10000L,
    alpha = 0.05,
    seed = 20260803L
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("Arguments must use the form --name=value.")
    }
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !nzchar(pieces[[1L]]) || !nzchar(pieces[[2L]])) {
      stop(sprintf("Invalid argument: %s", arg))
    }
    name <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    if (!name %in% names(defaults)) {
      stop(sprintf("Unknown argument: --%s", pieces[[1L]]))
    }
    defaults[[name]] <- pieces[[2L]]
  }

  if (is.null(defaults$output_dir)) {
    defaults$output_dir <- file.path(defaults$results_dir, "bhep_pvalues")
  }
  defaults$mc_rep <- as.integer(defaults$mc_rep)
  defaults$alpha <- as.numeric(defaults$alpha)
  defaults$seed <- as.integer(defaults$seed)

  if (is.na(defaults$mc_rep) || defaults$mc_rep < 1L) {
    stop("`--mc-rep` must be a positive integer.")
  }
  if (is.na(defaults$alpha) || defaults$alpha <= 0 || defaults$alpha >= 1) {
    stop("`--alpha` must lie strictly between zero and one.")
  }
  if (is.na(defaults$seed) || defaults$seed < 0L) {
    stop("`--seed` must be a non-negative integer.")
  }
  defaults
}

table_datasets <- data.frame(
  table_dataset = c(
    "Aar_oxides", "ArcticLake", "Boxite", "ClamEast", "ClamWest",
    "HouseholdExp", "Metabolites", "Sediments", "SerumProtein", "SkyeAFM",
    "Activity10", "Activity31", "AnimalVegetation", "Bayesite", "Coxite",
    "DiagnosticProb", "Firework", "Hongite", "Hydrochem", "juraset",
    "Kongite", "PogoJump", "ShiftOperators", "Supervisor",
    "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference",
    "Yatquat_panel", "SkyeLavas"
  ),
  result_slug = c(
    "aar_oxides", "arcticlake", "boxite", "clameast", "clamwest",
    "householdexp", "metabolites", "sediments", "serumprotein", "skyeafm",
    "activity10", "activity31", "animalvegetation", "bayesite", "coxite",
    "diagnosticprob", "firework", "hongite", "hydrochem", "juraset",
    "kongite", "pogojump", "shiftoperators", "supervisor",
    "whitecells_microscopic", "whitecells_image", "yatquat_preference",
    "yatquat_panel", "skyelavasaitchison32"
  ),
  expected_n = c(
    87L, 39L, 25L, 20L, 20L, 40L, 67L, 21L, 36L, 23L,
    20L, 20L, 100L, 21L, 25L, 30L, 81L, 25L, 485L, 359L,
    25L, 28L, 27L, 18L, 30L, 30L, 40L, 40L, 32L
  ),
  expected_D = c(
    10L, 3L, 5L, 6L, 6L, 4L, 3L, 3L, 4L, 3L,
    6L, 6L, 4L, 4L, 5L, 3L, 5L, 5L, 14L, 7L,
    5L, 3L, 4L, 4L, 3L, 3L, 3L, 3L, 10L
  ),
  stringsAsFactors = FALSE
)

load_and_validate_ilr <- function(dataset, results_dir) {
  result_path <- file.path(results_dir, sprintf("%s_results.rds", dataset$result_slug))
  if (!file.exists(result_path)) {
    stop(sprintf("Missing saved result for %s: %s", dataset$table_dataset, result_path))
  }
  result <- readRDS(result_path)
  z <- result$fit$Z
  if (is.null(z)) {
    stop(sprintf("Saved result for %s does not contain `fit$Z`.", dataset$table_dataset))
  }
  z <- as.matrix(z)
  storage.mode(z) <- "double"
  n <- nrow(z)
  d <- ncol(z)
  ambient_D <- result$data_prep$D %||% (d + 1L)

  if (!is.numeric(z) || any(!is.finite(z))) {
    stop(sprintf("The ilr matrix for %s must be finite and numeric.", dataset$table_dataset))
  }
  if (!identical(n, dataset$expected_n)) {
    stop(sprintf("Unexpected n for %s: %s (expected %s).", dataset$table_dataset, n, dataset$expected_n))
  }
  if (!identical(as.integer(ambient_D), dataset$expected_D) || d != dataset$expected_D - 1L) {
    stop(sprintf("Unexpected simplex/ilr dimensions for %s.", dataset$table_dataset))
  }
  if (n < d + 1L) {
    stop(sprintf("BHEP requires n >= d + 1 for %s.", dataset$table_dataset))
  }
  if (qr(sweep(z, 2L, colMeans(z), FUN = "-"))$rank != d) {
    stop(sprintf("The ilr sample covariance is singular for %s.", dataset$table_dataset))
  }

  list(z = z, n = n, ambient_D = as.integer(ambient_D), ilr_d = d, result_path = result_path)
}

run_bhep_sensitivity <- function(settings) {
  if (!requireNamespace("mnt", quietly = TRUE)) {
    stop("Package `mnt` is required. Install CRAN package mnt (version 1.3) before running this script.")
  }

  a_values <- c(0.5, 1, 2)
  mnt_version <- as.character(utils::packageVersion("mnt"))
  rows <- vector("list", nrow(table_datasets) * length(a_values))
  row_index <- 0L

  for (dataset_index in seq_len(nrow(table_datasets))) {
    dataset <- table_datasets[dataset_index, , drop = FALSE]
    message(sprintf("[%d/%d] BHEP sensitivity: %s", dataset_index, nrow(table_datasets), dataset$table_dataset))
    input <- load_and_validate_ilr(dataset, settings$results_dir)

    for (a_index in seq_along(a_values)) {
      a <- a_values[[a_index]]
      test_seed <- settings$seed + 100L * dataset_index + a_index
      set.seed(test_seed)
      statistic <- as.numeric(mnt::BHEP(input$z, a = a))
      null_statistics <- vapply(seq_len(settings$mc_rep), function(replication) {
        null_sample <- matrix(stats::rnorm(input$n * input$ilr_d), nrow = input$n, ncol = input$ilr_d)
        as.numeric(mnt::BHEP(null_sample, a = a))
      }, numeric(1))
      exceedances <- sum(null_statistics >= statistic)
      p_value <- (1 + exceedances) / (settings$mc_rep + 1)
      reject <- p_value <= settings$alpha

      if (!is.finite(statistic) || any(!is.finite(null_statistics)) || !is.finite(p_value) || p_value <= 0 || p_value > 1) {
        stop(sprintf("Non-finite BHEP Monte Carlo output for %s at a = %s.", dataset$table_dataset, a))
      }

      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        dataset = dataset$table_dataset,
        source_dataset = if (identical(dataset$table_dataset, "SkyeLavas")) "SkyeLavasAitchison32" else dataset$table_dataset,
        n = input$n,
        D = input$ambient_D,
        ilr_dimension = input$ilr_d,
        a = a,
        test_value = statistic,
        mc_exceedances = exceedances,
        p_value = p_value,
        reject = reject,
        seed = test_seed,
        alpha = settings$alpha,
        mc_rep = settings$mc_rep,
        mnt_version = mnt_version,
        input_rds = input$result_path,
        stringsAsFactors = FALSE
      )
    }
  }

  output <- do.call(rbind, rows)
  if (nrow(output) != 87L || length(unique(output$dataset)) != 29L || !identical(sort(unique(output$a)), a_values)) {
    stop("The BHEP sensitivity output must contain 29 datasets and three tuning parameters.")
  }
  if (any(!is.finite(output$test_value)) || any(!is.finite(output$p_value)) || any(output$p_value <= 0) || any(output$p_value > 1)) {
    stop("The BHEP sensitivity output contains non-finite values.")
  }
  if (!all(output$reject == (output$p_value <= output$alpha))) {
    stop("At least one BHEP decision does not match its Monte Carlo p-value.")
  }
  output
}

main <- function() {
  settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
  if (!dir.exists(settings$results_dir)) {
    stop(sprintf("Results directory does not exist: %s", settings$results_dir))
  }

  output <- run_bhep_sensitivity(settings)
  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_path <- file.path(settings$output_dir, "bhep_pvalues.csv")
  utils::write.csv(output, output_path, row.names = FALSE)

  message(sprintf("Wrote %s BHEP results for %s datasets to %s.", nrow(output), length(unique(output$dataset)), output_path))
  invisible(output)
}

if (sys.nframe() == 0L) {
  main()
}
