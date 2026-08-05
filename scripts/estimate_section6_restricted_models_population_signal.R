#!/usr/bin/env Rscript

# Known-parameter population signal diagnostic for the proposed restricted
# single-spiked Normal/LG null and the current Section 6 vMF 3.1 experiment.
#
# This script deliberately measures only the population distance-profile
# separation.  It does not implement, fit, or bootstrap the restricted model.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/estimate_section6_population_signal.R")

args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
dimensions <- parse_section6_csv(args$dimensions, c(2L, 5L), "integer")
n_centres <- as.integer(args$n_centres %||% 128L)
n_profile <- as.integer(args$n_profile %||% 8000L)
n_eval <- as.integer(args$n_eval %||% 4000L)
batches <- as.integer(args$batches %||% 4L)
cores <- as.integer(args$cores %||% 1L)
seed <- as.integer(args$seed %||% 20260826L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  "population_signal_restricted_normal_lg_vmf"
)

if (any(dimensions < 2L) || any(c(n_centres, n_profile, n_eval, batches, cores) < 1L)) {
  stop("Invalid population-signal controls.")
}

single_spike_parameters <- function(d) {
  list(
    theta = 0.5 * section6_e(d),
    a = 2 * d / (d + 1),
    b = d / (d + 1),
    Sigma = diag(c(2 * d / (d + 1), rep(d / (d + 1), d - 1L)))
  )
}

draw_law <- function(family, d, n, law = c("p0", "p1")) {
  law <- match.arg(law)
  if (family %in% c("normal_single_spike", "lg_single_spike")) {
    parameter <- single_spike_parameters(d)
    z <- mvtnorm::rmvnorm(n, mean = parameter$theta, sigma = parameter$Sigma)
    if (identical(law, "p1")) {
      from_q <- stats::runif(n) < 0.5
      if (any(from_q)) {
        z[from_q, ] <- mvtnorm::rmvnorm(
          sum(from_q), mean = -parameter$theta, sigma = parameter$Sigma
        )
      }
    }
    if (family == "normal_single_spike") return(z)
    return(logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L))
  }
  if (family == "vmf_antipodal") {
    mu <- section6_e(d + 1L)
    x <- rotasym::r_vMF(n, mu = mu, kappa = d)
    if (identical(law, "p1")) {
      from_q <- stats::runif(n) < 0.5
      if (any(from_q)) x[from_q, ] <- rotasym::r_vMF(sum(from_q), mu = -mu, kappa = d)
    }
    return(x)
  }
  stop("Unsupported family.")
}

distance_spec <- function(family) {
  switch(
    family,
    normal_single_spike = population_signal_spec("normal"),
    lg_single_spike = population_signal_spec("lg"),
    vmf_antipodal = population_signal_spec("vmf"),
    stop("Unsupported family.")
  )
}

one_batch <- function(job) {
  set.seed(job$seed)
  family <- as.character(job$family)
  d <- as.integer(job$d)
  centres <- draw_law(family, d, n_centres, "p0")
  x0_a <- draw_law(family, d, n_profile, "p0")
  x1_a <- draw_law(family, d, n_profile, "p1")
  x0_b <- draw_law(family, d, n_profile, "p0")
  x1_b <- draw_law(family, d, n_profile, "p1")
  x0_eval <- draw_law(family, d, n_eval, "p0")
  spec <- distance_spec(family)
  control <- section6_control(1000L, job$seed, 50L)
  d0_a <- spec$distance_matrix(x0_a, centres, control)
  d1_a <- spec$distance_matrix(x1_a, centres, control)
  d0_b <- spec$distance_matrix(x0_b, centres, control)
  d1_b <- spec$distance_matrix(x1_b, centres, control)
  d0_eval <- spec$distance_matrix(x0_eval, centres, control)
  values <- t(vapply(seq_len(ncol(d0_a)), function(j) {
    centre_signal(d0_a[, j], d1_a[, j], d0_b[, j], d1_b[, j], d0_eval[, j])
  }, numeric(2)))
  cbind(job, data.frame(
    S_mean_center = mean(values[, "S"]),
    S_median_center = stats::median(values[, "S"]),
    J = mean(values[, "J"]),
    stringsAsFactors = FALSE
  ))
}

summarize <- function(x) {
  group <- interaction(x$family, x$d, drop = TRUE)
  pieces <- lapply(split(x, group), function(cell) {
    out <- cell[1L, c("family", "d"), drop = FALSE]
    if (out$family %in% c("normal_single_spike", "lg_single_spike")) {
      parameter <- single_spike_parameters(out$d)
      out$theta_norm <- sqrt(sum(parameter$theta^2))
      out$spike_variance <- parameter$a
      out$orthogonal_variance <- parameter$b
    } else {
      out$theta_norm <- NA_real_
      out$spike_variance <- NA_real_
      out$orthogonal_variance <- NA_real_
    }
    out$kappa <- if (out$family == "vmf_antipodal") out$d else NA_real_
    for (metric in c("S_mean_center", "S_median_center", "J")) {
      out[[metric]] <- mean(cell[[metric]])
      out[[paste0(metric, "_mcse")]] <- stats::sd(cell[[metric]]) / sqrt(nrow(cell))
    }
    out
  })
  do.call(rbind, pieces)
}

grid <- expand.grid(
  family = c("normal_single_spike", "lg_single_spike", "vmf_antipodal"),
  d = dimensions, batch = seq_len(batches),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
grid$seed <- vapply(seq_len(nrow(grid)), function(i) {
  section6_seed(seed, i, grid$batch[[i]], stream = 261L)
}, integer(1))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_dir, "batch_results.csv")
summary_path <- file.path(output_dir, "summary.csv")
manifest_path <- file.path(output_dir, "manifest.txt")
if (!file.exists(manifest_path)) writeLines(c(
  "Known-parameter population signal for proposed restricted Normal/LG and current vMF 3.1.",
  "Normal/LG P0: N_d(theta0,Sigma0), theta0=0.5e1, a_d=2d/(d+1), b_d=d/(d+1).",
  "Normal/LG P1: 0.5 N_d(theta0,Sigma0) + 0.5 N_d(-theta0,Sigma0).",
  "vMF P0: vMF(e1,kappa=d); P1: 0.5 vMF(e1,d) + 0.5 vMF(-e1,d).",
  "J is split-sample Monte Carlo for E[(F1,Omega(T)-F0,Omega(T))^2].",
  "No MLE, score derivative, Vhat, or bootstrap is used.",
  sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
  sprintf("n_centres=%d; n_profile=%d; n_eval=%d; batches=%d; seed=%d",
    n_centres, n_profile, n_eval, batches, seed)
), manifest_path)

existing <- if (file.exists(raw_path)) utils::read.csv(raw_path, stringsAsFactors = FALSE) else data.frame()
key <- function(x) paste(x$family, x$d, x$batch, sep = "|")
pending <- grid[!key(grid) %in% if (nrow(existing)) key(existing) else character(), , drop = FALSE]
message(sprintf("Restricted-model population signal: %d pending batches.", nrow(pending)))
for (start in seq.int(1L, nrow(pending), by = cores)) {
  index <- seq.int(start, min(start + cores - 1L, nrow(pending)))
  jobs <- lapply(index, function(i) pending[i, , drop = FALSE])
  values <- if (length(jobs) == 1L) list(one_batch(jobs[[1L]])) else {
    parallel::mclapply(jobs, one_batch, mc.cores = min(cores, length(jobs)), mc.preschedule = FALSE)
  }
  existing <- rbind(existing, do.call(rbind, values))
  utils::write.csv(existing, raw_path, row.names = FALSE)
  utils::write.csv(summarize(existing), summary_path, row.names = FALSE)
  message(sprintf("completed batches: %d/%d", min(start + cores - 1L, nrow(pending)), nrow(pending)))
}
print(summarize(existing), row.names = FALSE)
