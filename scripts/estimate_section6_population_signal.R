#!/usr/bin/env Rscript

# Population (known-parameter) signal diagnostic for the Section 6 designs.
#
# For a null law P0 and the actual beta = 1 alternative P1, this script
# estimates, with independent Monte Carlo splits,
#
#   S = E_{Omega ~ P0}[ sup_t |F_1,Omega(t) - F_0,Omega(t)| ]
#   J = E_{Omega ~ P0, T ~ F_0,Omega}[{F_1,Omega(T) - F_0,Omega(T)}^2].
#
# The reported ``S_mean_center`` is the centre-averaged KS discrepancy.  It
# is deliberately reported alongside ``S_max_center``: the latter is a Monte
# Carlo lower-resolution proxy for the literal supremum over centres, whereas
# the former matches the integration over centres used by CvM.  J is estimated
# by a split-sample product, which removes the leading positive squared-ECDF
# bias conditional on each centre.
#
# No MLE, score derivative, bootstrap, or fitted profile is used here.  Thus
# this is a pure population diagnostic for the geometry/signal of the test.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

source("scripts/run_section6_new_scenarios.R")

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

population_signal_spec <- function(family) {
  switch(
    family,
    normal = make_mvnormal_spec(unknown_param = "both"),
    lg = make_logistic_gaussian_spec(unknown_param = "both"),
    vmf = make_vmf_spec(distance_type = "geodesic", unknown_param = "xi"),
    hvmf = make_hvmf_spec(unknown_param = "both"),
    stop("Unsupported family.")
  )
}

# In every Section 6 design P_beta = P0 + w(beta)(Q - P0).  The coefficient
# below is w(1), i.e. the discrepancy represented by an actual beta = 1 draw.
section6_beta_one_weight <- function(scenario) {
  if (scenario %in% c("normal_1_mixture", "lg_1_mixture", "vmf_1_antipodal", "hvmf_1_mixture")) {
    return(0.5)
  }
  if (scenario %in% c("normal_2_t3", "lg_2_t3", "vmf_2_projected_normal", "hvmf_2_angular")) {
    return(1)
  }
  stop(sprintf("Unknown scenario '%s'.", scenario))
}

empirical_cdf_at <- function(sorted_values, points) {
  findInterval(points, sorted_values, rightmost.closed = TRUE) / length(sorted_values)
}

centre_signal <- function(d0_a, d1_a, d0_b, d1_b, d0_eval) {
  s0_a <- sort(d0_a)
  s1_a <- sort(d1_a)
  s0_b <- sort(d0_b)
  s1_b <- sort(d1_b)

  # Averages two independent ECDF differences before taking the supremum.
  # This reduces, but cannot completely remove, the finite-MC upward bias of S.
  thresholds <- sort(unique(c(s0_a, s1_a, s0_b, s1_b)))
  delta_a <- empirical_cdf_at(s1_a, thresholds) - empirical_cdf_at(s0_a, thresholds)
  delta_b <- empirical_cdf_at(s1_b, thresholds) - empirical_cdf_at(s0_b, thresholds)
  s_hat <- max(abs(0.5 * (delta_a + delta_b)))

  # Conditional on the centre, delta_a(T) and delta_b(T) are independent and
  # have common expectation Delta(T); their product is therefore unbiased for
  # Delta(T)^2.  The independent evaluation draw implements dF_{0,omega}(T).
  delta_a_eval <- empirical_cdf_at(s1_a, d0_eval) - empirical_cdf_at(s0_a, d0_eval)
  delta_b_eval <- empirical_cdf_at(s1_b, d0_eval) - empirical_cdf_at(s0_b, d0_eval)
  j_hat <- mean(delta_a_eval * delta_b_eval)

  c(S = s_hat, J = j_hat)
}

one_population_batch <- function(design_row,
                                 n_centres,
                                 n_profile,
                                 n_eval,
                                 seed) {
  set.seed(seed)
  p0_row <- design_row
  p0_row$beta <- 0
  p1_row <- design_row
  p1_row$beta <- 1

  centres <- generate_section6_sample(transform(p0_row, n = as.integer(n_centres)))
  x0_a <- generate_section6_sample(transform(p0_row, n = as.integer(n_profile)))
  x1_a <- generate_section6_sample(transform(p1_row, n = as.integer(n_profile)))
  x0_b <- generate_section6_sample(transform(p0_row, n = as.integer(n_profile)))
  x1_b <- generate_section6_sample(transform(p1_row, n = as.integer(n_profile)))
  x0_eval <- generate_section6_sample(transform(p0_row, n = as.integer(n_eval)))

  spec <- population_signal_spec(as.character(design_row$family))
  control <- section6_control(derivative_mc_size = 1000L, derivative_seed = seed, cvm_block_size = 50L)
  d0_a <- spec$distance_matrix(x0_a, centres, control)
  d1_a <- spec$distance_matrix(x1_a, centres, control)
  d0_b <- spec$distance_matrix(x0_b, centres, control)
  d1_b <- spec$distance_matrix(x1_b, centres, control)
  d0_eval <- spec$distance_matrix(x0_eval, centres, control)

  metrics <- t(vapply(seq_len(ncol(d0_a)), function(j) {
    centre_signal(d0_a[, j], d1_a[, j], d0_b[, j], d1_b[, j], d0_eval[, j])
  }, numeric(2)))

  data.frame(
    S_mean_center = mean(metrics[, "S"]),
    S_median_center = stats::median(metrics[, "S"]),
    S_max_center = max(metrics[, "S"]),
    J = mean(metrics[, "J"]),
    J_min_center = min(metrics[, "J"]),
    stringsAsFactors = FALSE
  )
}

summarize_batches <- function(rows) {
  aggregate(
    cbind(S_mean_center, S_median_center, S_max_center, J, J_min_center) ~
      scenario + family + d + beta + beta_one_weight,
    data = rows,
    FUN = function(x) mean(x)
  ) -> estimate

  numeric_columns <- c("S_mean_center", "S_median_center", "S_max_center", "J", "J_min_center")
  for (column in numeric_columns) {
    se <- aggregate(rows[[column]],
      by = list(scenario = rows$scenario, family = rows$family, d = rows$d,
        beta = rows$beta, beta_one_weight = rows$beta_one_weight),
      FUN = function(x) stats::sd(x) / sqrt(length(x))
    )
    names(se)[ncol(se)] <- paste0(column, "_mcse")
    estimate <- merge(estimate, se,
      by = c("scenario", "family", "d", "beta", "beta_one_weight"),
      sort = FALSE
    )
  }
  estimate[order(estimate$family, estimate$scenario, estimate$d), , drop = FALSE]
}

run_population_signal <- function(dimensions = c(2L, 5L, 10L),
                                  families = section6_families,
                                  n_centres = 128L,
                                  n_profile = 8000L,
                                  n_eval = 4000L,
                                  batches = 4L,
                                  seed = 20260728L,
                                  output_dir,
                                  show_progress = TRUE) {
  dimensions <- sort(unique(as.integer(dimensions)))
  families <- intersect(section6_families, unique(as.character(families)))
  if (!length(families)) stop("No valid families requested.")
  if (any(dimensions < 2L)) stop("All dimensions must be at least two.")
  if (any(c(n_centres, n_profile, n_eval, batches) < 1L)) stop("All Monte Carlo sizes must be positive.")

  design <- do.call(rbind, lapply(families, function(family) {
    make_section6_design(family = family, dimensions = dimensions, n_values = 1L, beta_values = 1)[,
      c("scenario", "family", "alternative", "description", "d", "beta", "design_id"), drop = FALSE]
  }))
  design$design_id <- seq_len(nrow(design))
  jobs <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    cbind(design[rep(i, batches), , drop = FALSE], batch = seq_len(batches))
  }))
  rownames(jobs) <- NULL

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "batch_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  manifest_path <- file.path(output_dir, "manifest.txt")
  if (!file.exists(manifest_path)) {
    writeLines(c(
      "Known-parameter population signal diagnostic",
      "S_mean_center = E_{Omega~P0}[sup_t |F1,Omega(t)-F0,Omega(t)|] (Monte Carlo)",
      "J = E_{Omega~P0,T~F0,Omega}[(F1,Omega(T)-F0,Omega(T))^2] (split-sample Monte Carlo)",
      sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
      sprintf("families: %s", paste(families, collapse = ",")),
      sprintf("n_centres: %d", n_centres),
      sprintf("n_profile per independent split: %d", n_profile),
      sprintf("n_eval: %d", n_eval),
      sprintf("batches: %d", batches),
      sprintf("seed: %d", seed),
      "No MLE, score derivative, or bootstrap is used."
    ), manifest_path)
  }

  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else data.frame()
  job_key <- paste(jobs$scenario, jobs$d, jobs$beta, jobs$batch, sep = "|")
  done_key <- if (nrow(existing)) paste(existing$scenario, existing$d, existing$beta, existing$batch, sep = "|") else character()
  pending <- jobs[!job_key %in% done_key, , drop = FALSE]
  started <- proc.time()[["elapsed"]]

  if (isTRUE(show_progress)) {
    message(sprintf("Population signal: %d pending batches (%d already complete).", nrow(pending), nrow(jobs) - nrow(pending)))
  }
  for (i in seq_len(nrow(pending))) {
    job <- pending[i, , drop = FALSE]
    batch_seed <- section6_seed(seed, job$design_id, job$batch, stream = 91L)
    out <- one_population_batch(
      design_row = job, n_centres = n_centres, n_profile = n_profile,
      n_eval = n_eval, seed = batch_seed
    )
    out <- cbind(
      job[, c("scenario", "family", "alternative", "description", "d", "beta", "batch"), drop = FALSE],
      beta_one_weight = section6_beta_one_weight(as.character(job$scenario)),
      seed = batch_seed,
      out
    )
    existing <- rbind(existing, out)
    utils::write.csv(existing, result_path, row.names = FALSE)
    utils::write.csv(summarize_batches(existing), summary_path, row.names = FALSE)
    if (isTRUE(show_progress)) {
      elapsed <- proc.time()[["elapsed"]] - started
      done <- nrow(jobs) - nrow(pending) + i
      eta <- if (done > 0L) elapsed / i * (nrow(pending) - i) else NA_real_
      message(sprintf("completed batches: %d/%d; elapsed %.1fs; ETA %.1fs", done, nrow(jobs), elapsed, eta))
    }
  }
  invisible(list(batch_results = existing, summary = summarize_batches(existing)))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  output_dir <- args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios", "population_signal_known_parameters"
  )
  result <- run_population_signal(
    dimensions = parse_section6_csv(args$dimensions, c(2L, 5L, 10L), "integer"),
    families = parse_section6_csv(args$families, section6_families, "character"),
    n_centres = as.integer(args$n_centres %||% 128L),
    n_profile = as.integer(args$n_profile %||% 8000L),
    n_eval = as.integer(args$n_eval %||% 4000L),
    batches = as.integer(args$batches %||% 4L),
    seed = as.integer(args$seed %||% 20260728L),
    output_dir = output_dir,
    show_progress = parse_bool(args$show_progress, TRUE)
  )
  print(result$summary)
}
