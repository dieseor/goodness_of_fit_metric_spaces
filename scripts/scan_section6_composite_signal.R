#!/usr/bin/env Rscript

# Population drift diagnostic for the *composite* Section 6 tests.
#
# For an actual beta = 1 alternative P and its pseudo-true fitted-null law
# P_{theta_star}, this estimates
#
#  Sbar = E_{Omega~P} sup_t |F_{P,Omega}(t)-F_{theta_star,Omega}(t)|,
#  J    = E_{Omega,T~P} [F_{P,Omega}(T)-F_{theta_star,Omega}(T)]^2.
#
# This differs deliberately from the known-null signal screen: the centres
# and CvM thresholds are drawn from the alternative, as in the statistic, and
# the comparator is the KL/MLE projection rather than the generating null.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/scan_section6_dimension_scaled_candidates.R")

candidate_selected_default <- function() {
  c(
    "normal_bulk_mean_a0.60", "normal_t_from_t3_tailvar1.5",
    "lg_bulk_mean_a0.60", "lg_t_from_t3_tailvar1.5",
    "vmf_concentration_dquarter_c1.0", "vmf_local_location_dquarter_c0.9",
    "hvmf_local_location_dquarter_c1.1", "hvmf_local_angular_dquarter_c0.9"
  )
}

candidate_parameter_space_draw <- function(candidate, d, n, mean, Sigma = NULL) {
  z <- if (is.null(Sigma)) {
    matrix(stats::rnorm(n * d), nrow = n, ncol = d) +
      matrix(rep(mean, each = n), nrow = n, ncol = d)
  } else {
    mvtnorm::rmvnorm(n, mean = mean, sigma = Sigma)
  }
  if (candidate$family == "normal") return(z)
  logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L)
}

candidate_exact_pseudotrue <- function(candidate, d) {
  family <- as.character(candidate$family)
  rule <- as.character(candidate$rule)
  value <- as.numeric(candidate$value)

  if (family %in% c("normal", "lg")) {
    if (rule == "bulk_mean") {
      m <- value / (d - 1)^(1 / 4) * candidate_bulk_vector(d)
      return(list(
        draw = function(n) candidate_parameter_space_draw(candidate, d, n, mean = rep(0, d), Sigma = diag(d) + tcrossprod(m)),
        description = sprintf("exact Gaussian projection: mean=0; covariance=I+mm^T, ||m||=%.8f", sqrt(sum(m^2)))))
    }
    if (rule %in% c("bulk_scale", "bulk_scale_dquarter", "bulk_scale_dthird")) {
      return(list(
        draw = function(n) candidate_parameter_space_draw(candidate, d, n, mean = rep(0, d), Sigma = diag(d)),
        description = "exact Gaussian projection: mean=0; covariance=I"))
    }
    if (rule %in% c("t_df_linear", "t_df_from_t3")) {
      return(list(
        draw = function(n) candidate_parameter_space_draw(candidate, d, n, mean = rep(0, d), Sigma = diag(d)),
        description = "exact Gaussian projection: mean=0; covariance=I"))
    }
    if (rule == "t_df_from_t3_tailvar") {
      Sigma <- diag(c(1, rep(value, d - 1L)), d)
      return(list(
        draw = function(n) candidate_parameter_space_draw(candidate, d, n, mean = rep(0, d), Sigma = Sigma),
        description = sprintf("exact Gaussian projection: mean=0; covariance=diag(1,%.8f,...)", value)))
    }
  }

  if (family == "vmf") {
    q <- d
    mu0 <- candidate_unit_e(d + 1L)
    kappa <- 1.5 * d
    if (rule %in% c("concentration_mixture", "concentration_dquarter")) {
      rho <- if (rule == "concentration_mixture") value else value / d^(1 / 4)
      rbar <- 0.5 * (A_q(kappa * (1 - rho), q) + A_q(kappa * (1 + rho), q))
      kappa_star <- solve_vmf_kappa_from_rbar(rbar, q = q)
      return(list(
        draw = function(n) rotasym::r_vMF(n, mu = mu0, kappa = kappa_star),
        description = sprintf("exact vMF projection: kappa*=%.8f", kappa_star)))
    }
    if (rule %in% c("local_location", "local_location_dquarter")) {
      delta <- if (rule == "local_location") value / sqrt(d) else value / d^(1 / 4)
      mu1 <- cos(delta) * mu0 + sin(delta) * candidate_unit_e(d + 1L, 2L)
      resultant <- 0.5 * A_q(kappa, q) * (mu0 + mu1)
      rbar <- sqrt(sum(resultant^2))
      kappa_star <- solve_vmf_kappa_from_rbar(rbar, q = q)
      return(list(
        draw = function(n) rotasym::r_vMF(n, mu = resultant / rbar, kappa = kappa_star),
        description = sprintf("exact vMF projection: kappa*=%.8f; angular shift=%.8f", kappa_star, delta / 2)))
    }
    if (rule %in% c("symmetric_location_dquarter", "symmetric_location_sqrt")) {
      delta <- if (rule == "symmetric_location_dquarter") value / d^(1 / 4) else value / sqrt(d)
      rbar <- A_q(kappa, q) * cos(delta)
      kappa_star <- solve_vmf_kappa_from_rbar(rbar, q = q)
      return(list(
        draw = function(n) rotasym::r_vMF(n, mu = mu0, kappa = kappa_star),
        description = sprintf("exact symmetric vMF projection: kappa*=%.8f", kappa_star)))
    }
  }

  if (family == "hvmf" && rule %in% c("local_location", "local_location_dquarter")) {
    q <- d
    kappa <- q
    r <- if (rule == "local_location") value / sqrt(d) else value / d^(1 / 4)
    mu0 <- c(sqrt(2), candidate_unit_e(d))
    tangent <- c(0, candidate_unit_e(d, 2L))
    mu1 <- cosh(r) * mu0 + sinh(r) * tangent
    resultant <- 0.5 * hvmf_mean_resultant_ratio(q, kappa) * (mu0 + mu1)
    rbar <- sqrt(-hvmf_minkowski_inner_product(resultant, resultant))
    kappa_star <- hvmf_kappa_from_mean_resultant_ratio(q, rbar)
    mu_star <- resultant / rbar
    return(list(
      draw = function(n) rhvmf_polar(n, mu = mu_star, kappa = kappa_star),
      description = sprintf("exact HvMF projection: kappa*=%.8f; displacement=%.8f", kappa_star, r / 2)))
  }
  if (family == "hvmf" && rule %in% c("symmetric_location_dquarter", "symmetric_location_sqrt")) {
    q <- d
    kappa <- 1.5 * d
    r <- if (rule == "symmetric_location_dquarter") value / d^(1 / 4) else value / sqrt(d)
    mu0 <- c(sqrt(2), candidate_unit_e(d))
    rbar <- hvmf_mean_resultant_ratio(q, kappa) * cosh(r)
    kappa_star <- hvmf_kappa_from_mean_resultant_ratio(q, rbar)
    return(list(
      draw = function(n) rhvmf_polar(n, mu = mu0, kappa = kappa_star),
      description = sprintf("exact symmetric HvMF projection: kappa*=%.8f", kappa_star)))
  }
  NULL
}

# Only the conditional angular HvMF mixture lacks a short closed-form first
# moment in the present code.  Its MLE projection is therefore estimated once
# from a large independent draw, and the approximation is recorded.
candidate_mc_pseudotrue <- function(candidate, d, n_fit, seed) {
  set.seed(seed)
  x <- candidate_draw(candidate, d, n_fit, law = "p1")
  spec <- population_signal_spec(as.character(candidate$family))
  theta <- spec$fit_theta(x, null = list(type = "composite"), control = list())
  if (candidate$family != "hvmf") stop("MC pseudo-true projection is currently only used for HvMF.")
  list(
    draw = function(n) rhvmf_polar(n, mu = theta$mu, kappa = theta$kappa),
    description = sprintf("Monte Carlo HvMF projection (N=%d): kappa*=%.8f", n_fit, theta$kappa)
  )
}

candidate_composite_setup <- function(candidate, d, n_fit, seed) {
  exact <- candidate_exact_pseudotrue(candidate, d)
  if (!is.null(exact)) return(exact)
  candidate_mc_pseudotrue(candidate, d, n_fit = n_fit, seed = seed)
}

candidate_composite_one_batch <- function(candidate, d, setup,
                                          n_centres, n_profile, n_eval, seed) {
  set.seed(seed)
  centres <- candidate_draw(candidate, d, n_centres, "p1")
  x_theta_a <- setup$draw(n_profile)
  x_p_a <- candidate_draw(candidate, d, n_profile, "p1")
  x_theta_b <- setup$draw(n_profile)
  x_p_b <- candidate_draw(candidate, d, n_profile, "p1")
  x_p_eval <- candidate_draw(candidate, d, n_eval, "p1")
  spec <- population_signal_spec(as.character(candidate$family))
  control <- section6_control(derivative_mc_size = 1000L, derivative_seed = seed, cvm_block_size = 50L)
  d_theta_a <- spec$distance_matrix(x_theta_a, centres, control)
  d_p_a <- spec$distance_matrix(x_p_a, centres, control)
  d_theta_b <- spec$distance_matrix(x_theta_b, centres, control)
  d_p_b <- spec$distance_matrix(x_p_b, centres, control)
  d_p_eval <- spec$distance_matrix(x_p_eval, centres, control)
  values <- t(vapply(seq_len(ncol(d_theta_a)), function(j) {
    centre_signal(d_theta_a[, j], d_p_a[, j], d_theta_b[, j], d_p_b[, j], d_p_eval[, j])
  }, numeric(2)))
  data.frame(
    S_mean_center = mean(values[, "S"]), S_median_center = stats::median(values[, "S"]),
    S_max_center = max(values[, "S"]), J = mean(values[, "J"]), stringsAsFactors = FALSE
  )
}

run_composite_signal_scan <- function(candidates = candidate_selected_default(),
                                      dimensions = c(2L, 5L, 10L),
                                      n_centres = 96L,
                                      n_profile = 6000L,
                                      n_eval = 3000L,
                                      batches = 3L,
                                      pseudo_fit_size = 100000L,
                                      seed = 20260730L,
                                      output_dir,
                                      show_progress = TRUE) {
  catalog <- candidate_catalog()
  catalog <- catalog[catalog$candidate %in% candidates, , drop = FALSE]
  missing <- setdiff(candidates, catalog$candidate)
  if (length(missing)) stop(sprintf("Unknown candidate(s): %s", paste(missing, collapse = ", ")))
  dimensions <- sort(unique(as.integer(dimensions)))
  jobs <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    candidate <- catalog[i, , drop = FALSE]
    do.call(rbind, lapply(dimensions, function(d) {
      cbind(candidate[rep(1L, batches), , drop = FALSE], d = d, batch = seq_len(batches))
    }))
  }))
  rownames(jobs) <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(output_dir, "batch_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  projection_path <- file.path(output_dir, "pseudo_true_projections.csv")
  manifest_path <- file.path(output_dir, "manifest.txt")
  if (!file.exists(manifest_path)) writeLines(c(
    "Composite Section 6 population-drift diagnostic.",
    "Centres and CvM thresholds are sampled from the beta=1 alternative.",
    "The comparator is the MLE/KL pseudo-true fitted null.",
    "S is centre-averaged profile KS discrepancy; J is split-sample for the CvM integrand.",
    sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
    sprintf("n_centres: %d; n_profile: %d; n_eval: %d; batches: %d", n_centres, n_profile, n_eval, batches),
    sprintf("pseudo_fit_size for the angular HvMF approximation: %d", pseudo_fit_size),
    sprintf("seed: %d", seed)
  ), manifest_path)
  existing <- if (file.exists(result_path)) utils::read.csv(result_path, stringsAsFactors = FALSE) else data.frame()
  key <- function(x) paste(x$candidate, x$d, x$batch, sep = "|")
  pending <- jobs[!key(jobs) %in% if (nrow(existing)) key(existing) else character(), , drop = FALSE]
  if (isTRUE(show_progress)) message(sprintf("Composite candidate signal scan: %d pending batches.", nrow(pending)))
  projection_rows <- if (file.exists(projection_path)) utils::read.csv(projection_path, stringsAsFactors = FALSE) else data.frame()
  setup_cache <- list()
  for (i in seq_len(nrow(pending))) {
    job <- pending[i, , drop = FALSE]
    setup_key <- paste(job$candidate, job$d, sep = "|")
    if (is.null(setup_cache[[setup_key]])) {
      setup_seed <- section6_seed(seed, match(job$candidate, candidates), job$d, stream = 88L)
      setup_cache[[setup_key]] <- candidate_composite_setup(job, job$d, pseudo_fit_size, setup_seed)
      if (!setup_key %in% paste(projection_rows$candidate, projection_rows$d, sep = "|")) {
        projection_rows <- rbind(projection_rows, data.frame(
          candidate = job$candidate, family = job$family, d = job$d,
          projection = setup_cache[[setup_key]]$description, stringsAsFactors = FALSE
        ))
        utils::write.csv(projection_rows, projection_path, row.names = FALSE)
      }
    }
    batch_seed <- section6_seed(
      seed, match(job$candidate, candidates),
      as.integer(job$d) * 100L + as.integer(job$batch), stream = 89L
    )
    out <- candidate_composite_one_batch(
      job, job$d, setup_cache[[setup_key]], n_centres, n_profile, n_eval, batch_seed
    )
    out <- cbind(job, seed = batch_seed, out)
    existing <- rbind(existing, out)
    utils::write.csv(existing, result_path, row.names = FALSE)
    utils::write.csv(candidate_summarize(existing), summary_path, row.names = FALSE)
    if (isTRUE(show_progress)) message(sprintf("completed batches: %d/%d", nrow(existing), nrow(jobs)))
  }
  invisible(list(batch_results = existing, summary = candidate_summarize(existing)))
}

if (sys.nframe() == 0L) {
  args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
  candidates <- parse_section6_csv(args$candidates, candidate_selected_default(), "character")
  output_dir <- args$output_dir %||% file.path(
    "simulation_results", "section6_new_scenarios", "candidate_composite_signal_scan"
  )
  result <- run_composite_signal_scan(
    candidates = candidates,
    dimensions = parse_section6_csv(args$dimensions, c(2L, 5L, 10L), "integer"),
    n_centres = as.integer(args$n_centres %||% 96L),
    n_profile = as.integer(args$n_profile %||% 6000L),
    n_eval = as.integer(args$n_eval %||% 3000L),
    batches = as.integer(args$batches %||% 3L),
    pseudo_fit_size = as.integer(args$pseudo_fit_size %||% 100000L),
    seed = as.integer(args$seed %||% 20260730L),
    output_dir = output_dir,
    show_progress = parse_bool(args$show_progress, TRUE)
  )
  print(result$summary)
}
