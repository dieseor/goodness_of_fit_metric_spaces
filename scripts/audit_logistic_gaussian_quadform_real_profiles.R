#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory numerical audit only.  Its outputs must
# not be used to change production dispatching or paper results until the
# mathematical selection rule has been reviewed and explicitly approved.
# Numerical audit of the weighted noncentral chi-square CDF evaluations in the
# 30 compositional paper fits.  This is a diagnostic script only: it neither
# changes the production dispatcher nor reruns a bootstrap.

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list(
    results_dir = file.path(
      "real_data", "logistic_gaussian", "screening", "fast",
      "paper_table_B5000_sampleks_fast_rerun_20260718"
    ),
    output_dir = file.path(
      "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723",
      "quadform_real_profile_diagnostics"
    ),
    sample_per_other = 120L,
    all_datasets = c(
      "Coxite", "ClamCombined", "PogoJump", "Yatquat_preference",
      "Sediments", "AnimalVegetation"
    ),
    only_datasets = NULL,
    seed = 20260724L
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else ""
    if (identical(key, "results_dir")) out$results_dir <- value
    if (identical(key, "output_dir")) out$output_dir <- value
    if (identical(key, "sample_per_other")) out$sample_per_other <- as.integer(value)
    if (identical(key, "all_datasets")) out$all_datasets <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(key, "only_datasets")) out$only_datasets <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(key, "seed")) out$seed <- as.integer(value)
  }
  if (!is.finite(out$sample_per_other) || out$sample_per_other < 1L) {
    stop("`sample_per_other` must be a positive integer.")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("bootstrap", "model_specs.R"))

safe_eval <- function(expr) {
  warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(value = value, warnings = unique(warnings))
}

quadform_exact <- function(q, lambda, delta) {
  # The true quadratic form is strictly positive a.s. for a positive-definite
  # covariance, hence its CDF at the zero-radius diagonal is exactly zero.
  # Avoid asking numerical routines to resolve this trivial boundary value.
  if (q <= 1e-12) {
    return(c(
      farebrother_cdf = 0,
      farebrother_ifault = 0,
      farebrother_status = "analytic_zero",
      farebrother_warning = "",
      imhof_cdf = 0,
      imhof_abserr = 0,
      imhof_status = "analytic_zero",
      imhof_warning = ""
    ))
  }
  fare <- safe_eval(CompQuadForm::farebrother(
    q = q, lambda = lambda, h = rep.int(1L, length(lambda)), delta = delta,
    maxit = 100000L, eps = 1e-8
  ))
  imhof <- safe_eval(CompQuadForm::imhof(
    q = q, lambda = lambda, h = rep.int(1L, length(lambda)), delta = delta,
    epsabs = 1e-8, epsrel = 1e-8, limit = 20000L
  ))

  fare_value <- fare_ifault <- NA_real_
  fare_status <- "error"
  if (!inherits(fare$value, "error")) {
    fare_value <- 1 - as.numeric(fare$value$Qq)
    fare_ifault <- as.integer(fare$value$ifault %||% NA_integer_)
    fare_status <- if (is.na(fare_ifault) || fare_ifault == 0L) "ok" else "ifault"
  }
  imhof_value <- imhof_abserr <- NA_real_
  imhof_status <- "error"
  if (!inherits(imhof$value, "error")) {
    imhof_value <- 1 - as.numeric(imhof$value$Qq)
    imhof_abserr <- as.numeric(imhof$value$abserr %||% NA_real_)
    imhof_status <- if (is.finite(imhof_value) && imhof_value >= 0 && imhof_value <= 1) "ok" else "invalid"
  }
  c(
    farebrother_cdf = fare_value,
    farebrother_ifault = fare_ifault,
    farebrother_status = fare_status,
    farebrother_warning = paste(fare$warnings, collapse = " | "),
    imhof_cdf = imhof_value,
    imhof_abserr = imhof_abserr,
    imhof_status = imhof_status,
    imhof_warning = paste(imhof$warnings, collapse = " | ")
  )
}

hbe_features <- function(lambda, delta, q) {
  kappa1 <- sum(lambda * (1 + delta))
  kappa2 <- 2 * sum(lambda^2 * (1 + 2 * delta))
  kappa3 <- 8 * sum(lambda^3 * (1 + 3 * delta))
  nu_hbe <- 8 * kappa2^3 / kappa3^2
  scale_hbe <- sqrt(2 * nu_hbe / kappa2)
  q0 <- kappa1 - 2 * kappa2^2 / kappa3
  z <- scale_hbe * (q - kappa1) + nu_hbe
  list(
    kappa1 = kappa1,
    kappa2 = kappa2,
    kappa3 = kappa3,
    nu_hbe = nu_hbe,
    scale_hbe = scale_hbe,
    q0_hbe = q0,
    z_hbe = z,
    hbe_formula = stats::pgamma(z, shape = nu_hbe / 2, scale = 2)
  )
}

choose_auto_method <- function(lambda, delta, q, kappa1) {
  cond_number <- max(lambda) / min(lambda)
  q_ratio <- q / kappa1
  max_delta <- apply(as.matrix(delta), 1L, max)
  is_ill <- (cond_number > 1e4) |
    (length(lambda) >= 5L & max_delta > 10 & q_ratio > 0.1) |
    (max_delta > 1e3 & q_ratio > 1)
  ifelse(is_ill, "imhof", "farebrother")
}

sample_validation_indices <- function(n, q, z_hbe, stored_hbe, n_sample, all_cells) {
  total <- length(q)
  if (isTRUE(all_cells)) return(seq_len(total))
  positive <- which(q > 1e-12)
  if (!length(positive)) return(seq_len(min(total, n_sample)))
  cutoff <- positive[z_hbe[positive] <= 0]
  noncutoff <- positive[z_hbe[positive] > 0]
  # Deterministic strata: HBE's artificial-support region, its immediate
  # boundary, and the remainder.  This makes the selected comparisons useful
  # for diagnosing rather than estimating an unconditional mean error.
  near_boundary <- noncutoff[order(z_hbe[noncutoff])]
  parts <- list(
    utils::head(cutoff, ceiling(n_sample / 3)),
    utils::head(near_boundary, ceiling(n_sample / 3))
  )
  used <- unique(unlist(parts, use.names = FALSE))
  remaining <- setdiff(positive, used)
  n_remaining <- max(0L, n_sample - length(used))
  if (length(remaining) && n_remaining > 0L) {
    set.seed(100000L + total + n_sample)
    parts[[length(parts) + 1L]] <- sample(remaining, min(length(remaining), n_remaining))
  }
  unique(unlist(parts, use.names = FALSE))
}

result_paths <- list.files(args$results_dir, pattern = "_results[.]rds$", full.names = TRUE)
if (!length(result_paths)) stop("No paper result RDS files found in ", args$results_dir)
if (!is.null(args$only_datasets)) {
  dataset_labels <- vapply(result_paths, function(path) readRDS(path)$dataset_name, character(1))
  result_paths <- result_paths[dataset_labels %in% args$only_datasets]
  if (!length(result_paths)) stop("No selected datasets were found.")
}

dataset_rows <- list()
point_rows <- list()
point_counter <- 0L

for (path in sort(result_paths)) {
  result <- readRDS(path)
  raw <- result$bootstrap$raw_result
  observed <- raw$observed$cvm
  theta <- raw$observed$theta_hat
  dataset <- result$dataset_name
  n <- nrow(observed$distance_matrix)
  normalized <- normalize_logistic_gaussian_data(result$data_prep$X_closed)
  lambda <- as.numeric(theta$eigenvalues_full[theta$positive_idx])
  eigenvectors <- theta$eigenvectors_full[, theta$positive_idx, drop = FALSE]
  shifts <- matrix(rep(theta$mu_ilr, each = n), nrow = n) - normalized$ilr
  nu_matrix <- shifts %*% eigenvectors
  delta_matrix <- sweep(nu_matrix^2, 2L, lambda, "/")
  q_matrix <- observed$distance_matrix^2
  q <- as.vector(q_matrix)
  center <- rep.int(seq_len(n), n)
  threshold <- rep(seq_len(n), each = n)
  delta <- delta_matrix[center, , drop = FALSE]
  kappa1 <- rowSums(sweep(1 + delta, 2L, lambda, "*"))
  kappa2 <- 2 * rowSums(sweep(1 + 2 * delta, 2L, lambda^2, "*"))
  kappa3 <- 8 * rowSums(sweep(1 + 3 * delta, 2L, lambda^3, "*"))
  nu_hbe <- 8 * kappa2^3 / kappa3^2
  scale_hbe <- sqrt(2 * nu_hbe / kappa2)
  q0 <- kappa1 - 2 * kappa2^2 / kappa3
  z <- scale_hbe * (q - kappa1) + nu_hbe
  hbe_formula <- stats::pgamma(z, shape = nu_hbe / 2, scale = 2)
  hbe_stored <- as.vector(observed$theoretical_profile)
  auto_method <- choose_auto_method(lambda, delta, q, kappa1)
  positive_q <- q > 1e-12
  cutoff <- positive_q & z <= 0
  hbe_zero <- positive_q & hbe_stored == 0
  hbe_zero_extra <- hbe_zero & !cutoff
  hbe_formula_discrepancy <- abs(hbe_formula - hbe_stored)
  effective_rank_2 <- sum(lambda)^2 / sum(lambda^2)
  all_cells <- dataset %in% args$all_datasets
  selected <- sample_validation_indices(
    n = n, q = q, z_hbe = z, stored_hbe = hbe_stored,
    n_sample = args$sample_per_other, all_cells = all_cells
  )

  dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
    dataset = dataset,
    n = n,
    ilr_dim = length(lambda),
    profile_cells = length(q),
    positive_q_cells = sum(positive_q),
    lambda_min = min(lambda),
    lambda_max = max(lambda),
    lambda_condition = max(lambda) / min(lambda),
    lambda_effective_rank_2 = effective_rank_2,
    center_delta_max_max = max(apply(delta_matrix, 1L, max)),
    center_delta_max_q95 = as.numeric(stats::quantile(apply(delta_matrix, 1L, max), 0.95, names = FALSE)),
    center_delta_sum_max = max(rowSums(delta_matrix)),
    hbe_q0_positive_centers = sum(q0[seq_len(n)] > 0),
    hbe_q0_min = min(q0[seq_len(n)]),
    hbe_q0_median = stats::median(q0[seq_len(n)]),
    hbe_q0_max = max(q0[seq_len(n)]),
    analytic_hbe_cutoff_cells = sum(cutoff),
    analytic_hbe_cutoff_fraction_positive_q = mean(cutoff[positive_q]),
    stored_hbe_zero_cells = sum(hbe_zero),
    stored_hbe_zero_fraction_positive_q = mean(hbe_zero[positive_q]),
    stored_hbe_zero_outside_analytic_cutoff = sum(hbe_zero_extra),
    max_abs_stored_vs_hbe_formula = max(hbe_formula_discrepancy, na.rm = TRUE),
    auto_farebrother_cells = sum(auto_method == "farebrother"),
    auto_imhof_cells = sum(auto_method == "imhof"),
    auto_imhof_fraction = mean(auto_method == "imhof"),
    validation_scope = if (all_cells) "all_profile_cells" else "stratified_subset",
    validation_cells = length(selected),
    stringsAsFactors = FALSE
  )

  for (idx in selected) {
    i <- center[[idx]]
    exact <- quadform_exact(q = q[[idx]], lambda = lambda, delta = delta_matrix[i, ])
    point_counter <- point_counter + 1L
    point_rows[[point_counter]] <- data.frame(
      dataset = dataset,
      validation_scope = if (all_cells) "all_profile_cells" else "stratified_subset",
      center_index = i,
      threshold_index = threshold[[idx]],
      distance = sqrt(q[[idx]]),
      q = q[[idx]],
      q_over_mean = q[[idx]] / kappa1[[idx]],
      lambda_condition = max(lambda) / min(lambda),
      delta_max = max(delta_matrix[i, ]),
      delta_sum = sum(delta_matrix[i, ]),
      kappa1 = kappa1[[idx]],
      kappa2 = kappa2[[idx]],
      kappa3 = kappa3[[idx]],
      hbe_nu = nu_hbe[[idx]],
      hbe_q0 = q0[[idx]],
      hbe_z = z[[idx]],
      hbe_analytic_cutoff = z[[idx]] <= 0,
      hbe_profile_stored = hbe_stored[[idx]],
      hbe_profile_formula = hbe_formula[[idx]],
      auto_method = auto_method[[idx]],
      farebrother_cdf = as.numeric(exact[["farebrother_cdf"]]),
      farebrother_ifault = as.integer(exact[["farebrother_ifault"]]),
      farebrother_status = as.character(exact[["farebrother_status"]]),
      farebrother_warning = as.character(exact[["farebrother_warning"]]),
      imhof_cdf = as.numeric(exact[["imhof_cdf"]]),
      imhof_abserr = as.numeric(exact[["imhof_abserr"]]),
      imhof_status = as.character(exact[["imhof_status"]]),
      imhof_warning = as.character(exact[["imhof_warning"]]),
      stringsAsFactors = FALSE
    )
  }
  message(sprintf("%s: %d validation cells (%s)", dataset, length(selected), if (all_cells) "all" else "subset"))
}

dataset_summary <- do.call(rbind, dataset_rows)
points <- do.call(rbind, point_rows)
points$abs_hbe_minus_farebrother <- abs(points$hbe_profile_stored - points$farebrother_cdf)
points$abs_hbe_minus_imhof <- abs(points$hbe_profile_stored - points$imhof_cdf)
points$abs_farebrother_minus_imhof <- abs(points$farebrother_cdf - points$imhof_cdf)
points$farebrother_imhof_exceeds_imhof_abserr <- points$abs_farebrother_minus_imhof > pmax(1e-7, 5 * points$imhof_abserr)

validation_summary <- do.call(rbind, lapply(split(points, points$dataset), function(d) {
  data.frame(
    dataset = d$dataset[[1L]],
    validation_cells_exact = nrow(d),
    farebrother_ifault_nonzero = sum(!is.na(d$farebrother_ifault) & d$farebrother_ifault != 0L),
    farebrother_error = sum(d$farebrother_status == "error"),
    imhof_non_ok = sum(!d$imhof_status %in% c("ok", "analytic_zero")),
    imhof_abserr_max = max(d$imhof_abserr, na.rm = TRUE),
    farebrother_imhof_max_abs_difference = max(d$abs_farebrother_minus_imhof, na.rm = TRUE),
    farebrother_imhof_q95_abs_difference = as.numeric(stats::quantile(d$abs_farebrother_minus_imhof, 0.95, na.rm = TRUE, names = FALSE)),
    farebrother_imhof_disagreement_beyond_reported_imhof_error = sum(d$farebrother_imhof_exceeds_imhof_abserr, na.rm = TRUE),
    hbe_farebrother_max_abs_difference = max(d$abs_hbe_minus_farebrother, na.rm = TRUE),
    hbe_farebrother_q95_abs_difference = as.numeric(stats::quantile(d$abs_hbe_minus_farebrother, 0.95, na.rm = TRUE, names = FALSE)),
    hbe_farebrother_mean_abs_difference = mean(d$abs_hbe_minus_farebrother, na.rm = TRUE),
    hbe_farebrother_cutoff_cells_tested = sum(d$hbe_analytic_cutoff),
    hbe_farebrother_max_abs_difference_outside_cutoff = if (any(!d$hbe_analytic_cutoff)) max(d$abs_hbe_minus_farebrother[!d$hbe_analytic_cutoff], na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
dataset_summary <- merge(dataset_summary, validation_summary, by = "dataset", all.x = TRUE, sort = FALSE)
dataset_summary <- dataset_summary[match(vapply(dataset_rows, `[[`, character(1), "dataset"), dataset_summary$dataset), ]

utils::write.csv(dataset_summary, file.path(args$output_dir, "dataset_quadform_regime_summary.csv"), row.names = FALSE)
utils::write.csv(points, file.path(args$output_dir, "quadform_exact_validation_points.csv"), row.names = FALSE)
saveRDS(list(dataset_summary = dataset_summary, points = points, args = args), file.path(args$output_dir, "quadform_real_profile_diagnostics.rds"))

cat("Wrote ", file.path(args$output_dir, "dataset_quadform_regime_summary.csv"), "\n", sep = "")
cat("Wrote ", file.path(args$output_dir, "quadform_exact_validation_points.csv"), "\n", sep = "")
