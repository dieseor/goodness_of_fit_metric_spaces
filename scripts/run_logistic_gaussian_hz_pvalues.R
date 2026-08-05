#!/usr/bin/env Rscript

# Monte Carlo p-values for the Henze--Zirkler (HZ) version of the BHEP test
# on the 29 ilr datasets used in AoS Table 4. Datasets are processed
# sequentially. Within each dataset, the null replications are distributed
# over `--cores` forked processes (10 by default on macOS/Linux).
#
# The HZ statistic is mnt::HZ(), namely mnt::BHEP() with
# beta_{n,d} = ((2d + 1)n/4)^(1/(d + 4)) / sqrt(2).

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_arguments <- function(args) {
  settings <- list(
    results_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_dir = NULL,
    mc_rep = 10000L,
    alpha = 0.05,
    seed = 20260803L,
    cores = 10L
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) stop("Arguments must have the form --name=value.")
    value <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    if (length(value) != 2L || !value[[1L]] %in% names(settings)) {
      stop(sprintf("Unknown or invalid argument: %s", arg))
    }
    settings[[value[[1L]]]] <- value[[2L]]
  }
  if (is.null(settings$output_dir)) {
    settings$output_dir <- file.path(settings$results_dir, "hz_pvalues")
  }
  settings$mc_rep <- as.integer(settings$mc_rep)
  settings$alpha <- as.numeric(settings$alpha)
  settings$seed <- as.integer(settings$seed)
  settings$cores <- as.integer(settings$cores)
  if (is.na(settings$mc_rep) || settings$mc_rep < 1L || settings$mc_rep >= 1000000L ||
      is.na(settings$alpha) || settings$alpha <= 0 || settings$alpha >= 1 ||
      is.na(settings$seed) || settings$seed < 0L || is.na(settings$cores) || settings$cores < 1L) {
    stop("Invalid Monte Carlo settings.")
  }
  settings
}

datasets <- data.frame(
  dataset = c(
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

load_ilr <- function(row, results_dir) {
  path <- file.path(results_dir, sprintf("%s_results.rds", row$result_slug))
  if (!file.exists(path)) stop(sprintf("Missing saved result: %s", path))
  result <- readRDS(path)
  z <- as.matrix(result$fit$Z)
  storage.mode(z) <- "double"
  n <- nrow(z)
  d <- ncol(z)
  D <- result$data_prep$D %||% (d + 1L)
  if (!is.numeric(z) || any(!is.finite(z)) || n != row$expected_n ||
      as.integer(D) != row$expected_D || d != row$expected_D - 1L || n < d + 1L) {
    stop(sprintf("Invalid ilr input for %s.", row$dataset))
  }
  if (qr(sweep(z, 2L, colMeans(z), FUN = "-"))$rank != d) {
    stop(sprintf("Singular ilr covariance for %s.", row$dataset))
  }
  list(z = z, n = n, d = d, D = as.integer(D), path = path)
}

null_hz_statistics <- function(n, d, mc_rep, seed_base, cores) {
  replicate_seeds <- seed_base + seq_len(mc_rep)
  simulate_one <- function(replication) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(replicate_seeds[[replication]])
    z_null <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    as.numeric(mnt::HZ(z_null))
  }
  if (cores == 1L) {
    return(vapply(seq_len(mc_rep), simulate_one, numeric(1)))
  }
  as.numeric(unlist(parallel::mclapply(
    seq_len(mc_rep), simulate_one,
    mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE
  ), use.names = FALSE))
}

main <- function() {
  if (!requireNamespace("mnt", quietly = TRUE)) stop("Package `mnt` is required.")
  if (.Platform$OS.type == "windows") stop("This runner requires macOS or Linux for parallel::mclapply().")
  settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
  available_cores <- parallel::detectCores(logical = TRUE)
  if (is.na(available_cores) || available_cores < settings$cores) {
    stop(sprintf("Requested %d cores, but only %s are available.", settings$cores, available_cores))
  }
  if (!dir.exists(settings$results_dir)) stop(sprintf("Results directory does not exist: %s", settings$results_dir))

  rows <- vector("list", nrow(datasets))
  for (i in seq_len(nrow(datasets))) {
    row <- datasets[i, , drop = FALSE]
    message(sprintf("[%d/%d] HZ Monte Carlo p-value: %s (using %d cores)", i, nrow(datasets), row$dataset, settings$cores))
    input <- load_ilr(row, settings$results_dir)
    beta <- (((2 * input$d + 1) * input$n / 4)^(1 / (input$d + 4))) / sqrt(2)
    observed <- as.numeric(mnt::HZ(input$z))
    dataset_seed_base <- settings$seed + 1000000L * i
    null_statistics <- null_hz_statistics(input$n, input$d, settings$mc_rep, dataset_seed_base, settings$cores)
    exceedances <- sum(null_statistics >= observed)
    p_value <- (1 + exceedances) / (settings$mc_rep + 1)
    if (!is.finite(observed) || any(!is.finite(null_statistics)) || !is.finite(p_value)) {
      stop(sprintf("Non-finite Monte Carlo output for %s.", row$dataset))
    }
    rows[[i]] <- data.frame(
      dataset = row$dataset,
      source_dataset = if (identical(row$dataset, "SkyeLavas")) "SkyeLavasAitchison32" else row$dataset,
      n = input$n, D = input$D, ilr_dimension = input$d,
      a = beta, test_value = observed, mc_exceedances = exceedances,
      p_value = p_value, reject = p_value <= settings$alpha,
      seed_base = dataset_seed_base, alpha = settings$alpha,
      mc_rep = settings$mc_rep, cores = settings$cores,
      mnt_version = as.character(utils::packageVersion("mnt")),
      input_rds = input$path, stringsAsFactors = FALSE
    )
  }
  output <- do.call(rbind, rows)
  if (nrow(output) != 29L || length(unique(output$dataset)) != 29L ||
      any(!is.finite(output$test_value)) || any(!is.finite(output$p_value)) ||
      !all(output$reject == (output$p_value <= output$alpha))) {
    stop("Output validation failed.")
  }
  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(settings$output_dir, "hz_pvalues.csv")
  utils::write.csv(output, path, row.names = FALSE)
  message(sprintf("Wrote %d HZ p-values to %s", nrow(output), path))
}

main()
