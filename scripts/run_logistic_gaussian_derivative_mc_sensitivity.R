#!/usr/bin/env Rscript

# Recompute the 29 paper-table KS/CvM tests after increasing only the
# score-Monte-Carlo size used to approximate the parameter derivative dot F.
# The datasets are processed sequentially. The multiplier bootstrap itself may
# use `--n_cores` processes within each dataset.

resolve_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  hit <- candidates[file.exists(candidates) | dir.exists(candidates)]
  if (!length(hit)) stop(sprintf("Cannot find %s.", file.path(...)))
  hit[[1L]]
}

source(resolve_path("real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"))

parse_args <- function(args) {
  settings <- list(
    baseline_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_dir = NULL,
    derivative_mc_size = 10000L,
    n_cores = 5L
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) stop("Arguments must have the form --name=value.")
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[[1L]] %in% names(settings)) {
      stop(sprintf("Unknown or invalid argument: %s", arg))
    }
    settings[[pieces[[1L]]]] <- pieces[[2L]]
  }
  settings$derivative_mc_size <- as.integer(settings$derivative_mc_size)
  settings$n_cores <- as.integer(settings$n_cores)
  if (is.null(settings$output_dir)) {
    settings$output_dir <- file.path(
      dirname(settings$baseline_dir),
      sprintf("paper_table_derivative_score_mc_M%d_B5000_hbe", settings$derivative_mc_size)
    )
  }
  if (is.na(settings$derivative_mc_size) || settings$derivative_mc_size < 1L ||
      is.na(settings$n_cores) || settings$n_cores < 1L) {
    stop("Invalid settings.")
  }
  settings
}

table_datasets <- data.frame(
  display_name = c(
    "Aar_oxides", "ArcticLake", "Boxite", "ClamEast", "ClamWest",
    "HouseholdExp", "Metabolites", "Sediments", "SerumProtein", "SkyeAFM",
    "Activity10", "Activity31", "AnimalVegetation", "Bayesite", "Coxite",
    "DiagnosticProb", "Firework", "Hongite", "Hydrochem", "juraset",
    "Kongite", "PogoJump", "ShiftOperators", "Supervisor",
    "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference",
    "Yatquat_panel", "SkyeLavas"
  ),
  source_name = c(
    "Aar_oxides", "ArcticLake", "Boxite", "ClamEast", "ClamWest",
    "HouseholdExp", "Metabolites", "Sediments", "SerumProtein", "SkyeAFM",
    "Activity10", "Activity31", "AnimalVegetation", "Bayesite", "Coxite",
    "DiagnosticProb", "Firework", "Hongite", "Hydrochem", "juraset",
    "Kongite", "PogoJump", "ShiftOperators", "Supervisor",
    "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference",
    "Yatquat_panel", "SkyeLavasAitchison32"
  ),
  stored_slug = c(
    "aar_oxides", "arcticlake", "boxite", "clameast", "clamwest",
    "householdexp", "metabolites", "sediments", "serumprotein", "skyeafm",
    "activity10", "activity31", "animalvegetation", "bayesite", "coxite",
    "diagnosticprob", "firework", "hongite", "hydrochem", "juraset",
    "kongite", "pogojump", "shiftoperators", "supervisor",
    "whitecells_microscopic", "whitecells_image", "yatquat_preference",
    "yatquat_panel", "skyelavasaitchison32"
  ),
  stringsAsFactors = FALSE
)

read_baseline <- function(dir, slug) {
  path <- file.path(dir, sprintf("%s_results.rds", slug))
  if (!file.exists(path)) stop(sprintf("Missing baseline result: %s", path))
  result <- readRDS(path)
  diagnostics <- result$bootstrap$raw_result$diagnostics
  if (!identical(as.integer(result$settings$B), 5000L) ||
      !identical(result$bootstrap$mode, "composite_multiplier") ||
      !identical(result$bootstrap$bootstrap_method, "fast_multiplier") ||
      !identical(result$settings$omega_grid_type, "sample_points") ||
      !identical(result$settings$t_grid_type, "sample_distances") ||
      !identical(result$settings$quadform_method, "hbe") ||
      !identical(diagnostics$derivative_method, "score_mc") ||
      !identical(as.integer(diagnostics$derivative_mc_size), 1000L)) {
    stop(sprintf("Baseline result for %s is not the expected M=1000 paper-table run.", slug))
  }
  result
}

main <- function() {
  settings <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!dir.exists(settings$baseline_dir)) stop("Baseline directory does not exist.")
  if (.Platform$OS.type == "unix" && settings$n_cores > parallel::detectCores(logical = TRUE)) {
    stop("Requested more cores than available.")
  }

  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", nrow(table_datasets))
  for (i in seq_len(nrow(table_datasets))) {
    item <- table_datasets[i, , drop = FALSE]
    baseline <- read_baseline(settings$baseline_dir, item$stored_slug)
    message(sprintf("[%d/%d] %s: derivative score MC M=%d, B=5000", i, nrow(table_datasets), item$display_name, settings$derivative_mc_size))

    # Preserve every paper-table setting except the derivative approximation.
    control <- list(
      logistic_gaussian_quadform_method = "hbe",
      derivative_method = "score_mc",
      derivative_mc_size = settings$derivative_mc_size
    )
    rerun <- run_logistic_gaussian_screening(
      dataset_name = item$source_name,
      B = 5000L,
      max_centers = baseline$settings$max_centers,
      n_t = baseline$settings$n_t,
      t_grid_tail_prob = baseline$settings$t_grid_tail_prob,
      boundary_epsilon = baseline$settings$boundary_epsilon,
      probs_t = baseline$settings$probs_t,
      bootstrap_mode = "composite_multiplier",
      seed = baseline$settings$seed,
      alpha = baseline$settings$alpha,
      ridge = baseline$settings$ridge,
      n_cores = settings$n_cores,
      bootstrap_method = "fast_multiplier",
      control = control,
      omega_grid_type = "sample_points",
      t_grid_type = "sample_distances",
      make_plots = FALSE,
      save_outputs = TRUE,
      output_dir = settings$output_dir,
      run_seed_sensitivity = FALSE,
      verbose = TRUE
    )
    new_diagnostics <- rerun$bootstrap$raw_result$diagnostics
    if (!identical(new_diagnostics$derivative_method, "score_mc") ||
        !identical(as.integer(new_diagnostics$derivative_mc_size), settings$derivative_mc_size)) {
      stop(sprintf("Derivative configuration was not honoured for %s.", item$display_name))
    }
    rows[[i]] <- data.frame(
      dataset = item$display_name,
      source_dataset = item$source_name,
      n = rerun$data_prep$n,
      D = rerun$data_prep$D,
      B = rerun$settings$B,
      derivative_method = new_diagnostics$derivative_method,
      derivative_mc_size = new_diagnostics$derivative_mc_size,
      derivative_mc_seed = new_diagnostics$derivative_mc_seed,
      ks_pvalue_baseline = baseline$inference$ks$p_value,
      ks_pvalue_M10000 = rerun$inference$ks$p_value,
      ks_pvalue_change = rerun$inference$ks$p_value - baseline$inference$ks$p_value,
      cvm_pvalue_baseline = baseline$inference$cvm$p_value,
      cvm_pvalue_M10000 = rerun$inference$cvm$p_value,
      cvm_pvalue_change = rerun$inference$cvm$p_value - baseline$inference$cvm$p_value,
      stringsAsFactors = FALSE
    )
  }
  comparison <- do.call(rbind, rows)
  if (nrow(comparison) != 29L || any(!is.finite(as.matrix(comparison[, c(
    "ks_pvalue_baseline", "ks_pvalue_M10000", "cvm_pvalue_baseline", "cvm_pvalue_M10000"
  )])))) stop("Comparison validation failed.")
  path <- file.path(settings$output_dir, "derivative_mc_M10000_comparison.csv")
  utils::write.csv(comparison, path, row.names = FALSE)
  message(sprintf("Wrote comparison for %d datasets to %s", nrow(comparison), path))
}

main()
