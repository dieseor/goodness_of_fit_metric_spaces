#!/usr/bin/env Rscript

# Production runner for the 29 compositional datasets in the AoS table.
# Datasets are processed sequentially; each multiplier-bootstrap computation
# uses `--n_cores` workers.  Auxiliary screening profiles and univariate or
# multivariate normality diagnostics are deliberately disabled: this runner
# computes only the inferential KS/CvM statistics and bootstrap p-values.

resolve_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  hit <- candidates[file.exists(candidates) | dir.exists(candidates)]
  if (!length(hit)) stop(sprintf("Cannot find %s.", file.path(...)))
  hit[[1L]]
}

source(resolve_path("real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"))

parse_args <- function(args) {
  settings <- list(
    output_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_defaults_score_mc_M10000_B5000_cpp_shared"
    ),
    B = 5000L,
    derivative_mc_size = 10000L,
    n_cores = 5L,
    seed = 20260803L
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) stop("Arguments must have the form --name=value.")
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[[1L]] %in% names(settings)) {
      stop(sprintf("Unknown or invalid argument: %s", arg))
    }
    settings[[pieces[[1L]]]] <- pieces[[2L]]
  }
  for (name in c("B", "derivative_mc_size", "n_cores", "seed")) {
    settings[[name]] <- as.integer(settings[[name]])
  }
  if (any(!is.finite(unlist(settings[c("B", "derivative_mc_size", "n_cores", "seed")]))) ||
      settings$B < 1L || settings$derivative_mc_size < 1L || settings$n_cores < 1L) {
    stop("`B`, `derivative_mc_size`, `n_cores`, and `seed` must be positive integers.")
  }
  settings
}

paper_datasets <- data.frame(
  dataset = c(
    "Aar_oxides", "ArcticLake", "Boxite", "ClamEast", "ClamWest",
    "HouseholdExp", "Metabolites", "Sediments", "SerumProtein", "SkyeAFM",
    "Activity10", "Activity31", "AnimalVegetation", "Bayesite", "Coxite",
    "DiagnosticProb", "Firework", "Hongite", "Hydrochem", "juraset",
    "Kongite", "PogoJump", "ShiftOperators", "Supervisor",
    "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference",
    "Yatquat_panel", "SkyeLavas"
  ),
  source_dataset = c(
    "Aar_oxides", "ArcticLake", "Boxite", "ClamEast", "ClamWest",
    "HouseholdExp", "Metabolites", "Sediments", "SerumProtein", "SkyeAFM",
    "Activity10", "Activity31", "AnimalVegetation", "Bayesite", "Coxite",
    "DiagnosticProb", "Firework", "Hongite", "Hydrochem", "juraset",
    "Kongite", "PogoJump", "ShiftOperators", "Supervisor",
    "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference",
    "Yatquat_panel", "SkyeLavasAitchison32"
  ),
  stringsAsFactors = FALSE
)

scalar <- function(x, default = NA) {
  if (is.null(x) || !length(x)) default else x[[1L]]
}

validate_result <- function(result, item, settings, dataset_seed, derivative_seed,
                            expected_n = NULL) {
  diagnostics <- result$bootstrap$raw_result$diagnostics
  failures <- character()
  check <- function(ok, message) if (!isTRUE(ok)) failures <<- c(failures, message)

  check(identical(result$bootstrap$mode, "composite_multiplier"), "wrong bootstrap mode")
  check(identical(result$bootstrap$bootstrap_method, "fast_multiplier"), "wrong bootstrap method")
  check(identical(as.integer(result$settings$B), settings$B), "wrong B")
  check(identical(as.integer(result$settings$seed), dataset_seed), "wrong dataset seed")
  check(identical(result$settings$compute_auxiliary_diagnostics, FALSE),
        "auxiliary diagnostics were not disabled")
  check(identical(result$grid$omega_grid_construction, "not_computed"),
        "auxiliary screening grid was unexpectedly computed")
  check(identical(diagnostics$fast_ks_mode, "sample_points_unique_distances_streamed"),
        "KS did not use the complete sample-based grid")
  check(identical(diagnostics$fast_cvm_mode, "sample_points_unique_distances_sorted_rows"),
        "CvM did not use the complete sample-based grid")
  check(isTRUE(diagnostics$shared_sample_ks_cvm_cache),
        "KS/CvM did not share the sample-based grid preparation")
  check(identical(diagnostics$derivative_method_effective, "score_mc"), "dot F is not score Monte Carlo")
  check(identical(as.integer(diagnostics$derivative_mc_size), settings$derivative_mc_size),
        "wrong derivative Monte Carlo size")
  check(identical(as.integer(diagnostics$derivative_mc_seed), derivative_seed),
        "wrong derivative Monte Carlo seed")
  check(identical(diagnostics$fast_multiplier_backend_effective, "cpp"),
        "C++ multiplier backend was not used")
  check(identical(diagnostics$fast_multiplier_cpp_kernel_effective, "contiguous_double"),
        "wrong effective C++ kernel")
  check(isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective),
        "KS/CvM fusion was not used")
  check(identical(diagnostics$fast_multiplier_cache_corrections_requested, "auto"),
        "wrong correction-cache policy")
  check(is.logical(diagnostics$fast_multiplier_cache_corrections_effective) &&
          length(diagnostics$fast_multiplier_cache_corrections_effective) == 1L &&
          !is.na(diagnostics$fast_multiplier_cache_corrections_effective),
        "invalid effective correction-cache decision")
  check(identical(result$settings$quadform_method, "auto"), "wrong quadratic-form policy")
  check(all(is.finite(c(result$inference$ks$p_value, result$inference$cvm$p_value))),
        "non-finite p-value")
  if (!is.null(expected_n)) check(identical(as.integer(result$data_prep$n), as.integer(expected_n)), "wrong sample size")
  if (identical(item$source_dataset, "SkyeLavasAitchison32")) {
    check(identical(as.integer(result$data_prep$n), 32L), "SkyeLavas is not the 32-observation version")
  }
  if (length(failures)) stop(sprintf("Validation failed for %s: %s", item$dataset, paste(failures, collapse = "; ")))
  invisible(TRUE)
}

result_row <- function(result, item, dataset_seed, derivative_seed) {
  d <- result$bootstrap$raw_result$diagnostics
  data.frame(
    dataset = item$dataset,
    source_dataset = item$source_dataset,
    n = result$data_prep$n,
    D = result$data_prep$D,
    ilr_dimension = ncol(result$fit$Z),
    B = result$settings$B,
    alpha = result$settings$alpha,
    dataset_seed = dataset_seed,
    derivative_method = scalar(d$derivative_method_effective),
    derivative_mc_size = scalar(d$derivative_mc_size),
    derivative_mc_seed = derivative_seed,
    quadform_method = result$settings$quadform_method,
    ridge = result$settings$ridge,
    ridge_added = result$fit$ridge_added,
    omega_grid = result$settings$omega_grid_type,
    omega_grid_size = result$grid$omega_grid_size,
    bootstrap_method = scalar(d$effective_bootstrap_method),
    multiplier_backend = scalar(d$fast_multiplier_backend_effective),
    cpp_kernel = scalar(d$fast_multiplier_cpp_kernel_effective),
    fused_ks_cvm = scalar(d$fast_multiplier_fuse_ks_cvm_effective),
    correction_cache_requested = scalar(d$fast_multiplier_cache_corrections_requested),
    correction_cache_effective = scalar(d$fast_multiplier_cache_corrections_effective),
    correction_cache_bytes = scalar(d$sample_correction_cache_bytes),
    shared_correction_cache = scalar(d$shared_sample_correction_cache),
    distance_profile_backend = scalar(d$distance_profile_backend_effective),
    ks_pvalue = result$inference$ks$p_value,
    cvm_pvalue = result$inference$cvm$p_value,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  settings <- parse_args(commandArgs(trailingOnly = TRUE))
  if (.Platform$OS.type == "unix" && settings$n_cores > parallel::detectCores(logical = TRUE)) {
    stop("Requested more cores than are available.")
  }
  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      datasets = paper_datasets,
      settings = settings,
      fixed_control = list(
        derivative_method = "score_mc",
        mvnormal_quadform_method = "auto",
        fast_multiplier_backend = "cpp",
        fast_multiplier_cpp_kernel = "contiguous_double",
        fast_multiplier_fuse_ks_cvm = TRUE,
        fast_multiplier_cache_corrections = "auto",
        compute_auxiliary_diagnostics = FALSE
      ),
      created_at = Sys.time(),
      R_version = R.version.string
    ),
    file.path(settings$output_dir, "run_manifest.rds")
  )

  rows <- vector("list", nrow(paper_datasets))
  for (i in seq_len(nrow(paper_datasets))) {
    item <- paper_datasets[i, , drop = FALSE]
    dataset_seed <- as.integer(settings$seed + 100L * (i - 1L))
    derivative_seed <- as.integer(settings$seed + 50000L + i)
    prepared <- prepare_composition_dataset(item$source_dataset)
    if (!identical(prepared$status, "ok")) stop(sprintf("Could not prepare %s.", item$source_dataset))
    result_path <- file.path(
      settings$output_dir,
      sprintf("%s_results.rds", slugify_dataset_name(item$source_dataset))
    )

    if (file.exists(result_path)) {
      message(sprintf("[%d/%d] %s: reusing validated result", i, nrow(paper_datasets), item$dataset))
      result <- readRDS(result_path)
      if (!identical(result$settings$compute_auxiliary_diagnostics, FALSE)) {
        message(sprintf("[%d/%d] %s: removing legacy auxiliary artifacts", i, nrow(paper_datasets), item$dataset))
        result <- strip_auxiliary_screening_artifacts(result)
        saveRDS(result, file = result_path)
      }
    } else {
      message(sprintf("[%d/%d] %s: B=%d, M=%d, cores=%d", i, nrow(paper_datasets), item$dataset,
                      settings$B, settings$derivative_mc_size, settings$n_cores))
      result <- run_logistic_gaussian_screening(
        dataset_name = item$source_dataset,
        B = settings$B,
        # Compatibility with the current internal API; this is exactly n, not a cap.
        max_centers = prepared$n,
        bootstrap_mode = "composite_multiplier",
        seed = dataset_seed,
        alpha = 0.05,
        ridge = 1e-8,
        n_cores = settings$n_cores,
        bootstrap_method = "fast_multiplier",
        bootstrap_keep = list(
          observed_process = FALSE,
          bootstrap_statistics = TRUE,
          bootstrap_thetas = FALSE
        ),
        compute_auxiliary_diagnostics = FALSE,
        control = list(
          derivative_method = "score_mc",
          derivative_mc_size = settings$derivative_mc_size,
          derivative_mc_seed = derivative_seed,
          mvnormal_quadform_method = "auto",
          fast_multiplier_backend = "cpp",
          fast_multiplier_cpp_kernel = "contiguous_double",
          fast_multiplier_fuse_ks_cvm = TRUE,
          fast_multiplier_cache_corrections = "auto"
        ),
        omega_grid_type = "sample_points",
        t_grid_type = "sample_distances",
        make_plots = FALSE,
        save_outputs = TRUE,
        output_dir = settings$output_dir,
        run_seed_sensitivity = FALSE,
        verbose = TRUE
      )
    }
    validate_result(result, item, settings, dataset_seed, derivative_seed, expected_n = prepared$n)
    rows[[i]] <- result_row(result, item, dataset_seed, derivative_seed)
    utils::write.csv(do.call(rbind, rows[seq_len(i)]),
                     file.path(settings$output_dir, "ks_cvm_pvalues.csv"), row.names = FALSE)
  }
  results <- do.call(rbind, rows)
  if (nrow(results) != 29L || anyDuplicated(results$dataset) ||
      results$n[results$dataset == "SkyeLavas"] != 32L) {
    stop("Final validation of the 29-dataset table failed.")
  }
  message(sprintf("Completed %d datasets. Results: %s", nrow(results),
                  file.path(settings$output_dir, "ks_cvm_pvalues.csv")))
}

main()
