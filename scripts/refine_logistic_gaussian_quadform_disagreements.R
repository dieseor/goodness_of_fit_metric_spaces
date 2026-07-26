#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory numerical audit only.  It must not alter
# the production dispatcher or paper results without explicit approval.
# Third-route refinement for profile cells at which the default Farebrother and
# Imhof calls disagree.  Diagnostic only; it does not change production code.

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list(
    validation_csv = file.path(
      "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723",
      "quadform_real_profile_focus_validation", "quadform_exact_validation_points.csv"
    ),
    results_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_csv = file.path(
      "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723",
      "quadform_real_profile_focus_validation", "quadform_exact_refinement.csv"
    ),
    abs_threshold = 1e-6
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else ""
    if (identical(key, "validation_csv")) out$validation_csv <- value
    if (identical(key, "results_dir")) out$results_dir <- value
    if (identical(key, "output_csv")) out$output_csv <- value
    if (identical(key, "abs_threshold")) out$abs_threshold <- as.numeric(value)
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
source(file.path("bootstrap", "model_specs.R"))

points <- utils::read.csv(args$validation_csv, stringsAsFactors = FALSE)
points <- points[
  (is.finite(points$abs_farebrother_minus_imhof) & points$abs_farebrother_minus_imhof > args$abs_threshold) |
    (!is.na(points$farebrother_ifault) & points$farebrother_ifault != 0L),
  , drop = FALSE
]
if (!nrow(points)) {
  utils::write.csv(data.frame(), args$output_csv, row.names = FALSE)
  quit(status = 0L)
}

load_case <- function(dataset) {
  paths <- list.files(args$results_dir, pattern = "_results[.]rds$", full.names = TRUE)
  idx <- vapply(paths, function(path) identical(readRDS(path)$dataset_name, dataset), logical(1))
  if (sum(idx) != 1L) stop("Could not identify exactly one RDS for ", dataset)
  result <- readRDS(paths[[which(idx)]])
  raw <- result$bootstrap$raw_result
  theta <- raw$observed$theta_hat
  z <- normalize_logistic_gaussian_data(result$data_prep$X_closed)$ilr
  lambda <- as.numeric(theta$eigenvalues_full[theta$positive_idx])
  eigenvectors <- theta$eigenvectors_full[, theta$positive_idx, drop = FALSE]
  list(
    raw = raw, theta = theta, z = z, lambda = lambda,
    eigenvectors = eigenvectors
  )
}

case_cache <- list()
refine_one <- function(row) {
  dataset <- row$dataset[[1L]]
  if (is.null(case_cache[[dataset]])) case_cache[[dataset]] <<- load_case(dataset)
  case <- case_cache[[dataset]]
  i <- as.integer(row$center_index[[1L]])
  j <- as.integer(row$threshold_index[[1L]])
  q <- case$raw$observed$cvm$distance_matrix[i, j]^2
  shift <- case$theta$mu_ilr - case$z[i, ]
  delta <- as.numeric((as.numeric(shift %*% case$eigenvectors)^2) / case$lambda)
  h <- rep.int(1L, length(case$lambda))

  fare_1 <- CompQuadForm::farebrother(q, case$lambda, h, delta, maxit = 1000000L, eps = 1e-12, mode = 1L)
  fare_0 <- CompQuadForm::farebrother(q, case$lambda, h, delta, maxit = 1000000L, eps = 1e-12, mode = 0L)
  imhof_tight <- CompQuadForm::imhof(q, case$lambda, h, delta, epsabs = 1e-12, epsrel = 1e-12, limit = 1000000L)
  davies <- CompQuadForm::davies(q, case$lambda, h, delta, lim = 1000000L, acc = 1e-8)

  fare_1_cdf <- 1 - fare_1$Qq
  fare_0_cdf <- 1 - fare_0$Qq
  imhof_cdf <- 1 - imhof_tight$Qq
  davies_cdf <- 1 - davies$Qq
  reference_rule <- if (fare_1$ifault == 0L && fare_0$ifault == 0L && davies$ifault == 0L &&
                        max(abs(fare_1_cdf - fare_0_cdf), abs(fare_1_cdf - davies_cdf)) <= 1e-6) {
    "farebrother_modes_and_davies_agree"
  } else if (davies$ifault == 0L && abs(imhof_cdf - davies_cdf) <= 1e-6) {
    "imhof_and_davies_agree"
  } else {
    "unresolved"
  }

  data.frame(
    dataset = dataset,
    center_index = i,
    threshold_index = j,
    q = q,
    lambda_condition = max(case$lambda) / min(case$lambda),
    original_farebrother_cdf = row$farebrother_cdf[[1L]],
    original_farebrother_ifault = row$farebrother_ifault[[1L]],
    original_imhof_cdf = row$imhof_cdf[[1L]],
    original_imhof_abserr = row$imhof_abserr[[1L]],
    farebrother_mode1_tight_cdf = fare_1_cdf,
    farebrother_mode1_tight_ifault = fare_1$ifault,
    farebrother_mode0_tight_cdf = fare_0_cdf,
    farebrother_mode0_tight_ifault = fare_0$ifault,
    imhof_tight_cdf = imhof_cdf,
    imhof_tight_abserr = imhof_tight$abserr,
    davies_cdf = davies_cdf,
    davies_ifault = davies$ifault,
    max_fare_mode_difference = abs(fare_1_cdf - fare_0_cdf),
    fare_vs_davies_difference = abs(fare_1_cdf - davies_cdf),
    imhof_vs_davies_difference = abs(imhof_cdf - davies_cdf),
    reference_rule = reference_rule,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, lapply(split(points, seq_len(nrow(points))), refine_one))
dir.create(dirname(args$output_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, args$output_csv, row.names = FALSE)
cat("Refined ", nrow(out), " profile cells; wrote ", args$output_csv, "\n", sep = "")
