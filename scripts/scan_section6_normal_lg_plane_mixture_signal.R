#!/usr/bin/env Rscript

# Known-parameter population signal scan for the two plane-mixture scenarios
# 1.1 (Normal) and 2.1 (logistic Gaussian).  This is deliberately separate
# from the bootstrap experiment: it contains no fitted parameters, score,
# derivative or bootstrap approximation.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
source("scripts/estimate_section6_population_signal.R")

args <- parse_section6_args(commandArgs(trailingOnly = TRUE))
csv_numeric <- function(name, default) {
  parse_section6_csv(args[[name]], default, "numeric")
}
csv_integer <- function(name, default) {
  parse_section6_csv(args[[name]], default, "integer")
}

dimensions <- csv_integer("dimensions", c(2L, 5L))
locations <- csv_numeric("locations", c(0.5, 0.75, 1.0))
correlations <- csv_numeric("correlations", c(0.25, 0.5, 0.75))
n_centres <- as.integer(args$n_centres %||% 128L)
n_profile <- as.integer(args$n_profile %||% 6000L)
n_eval <- as.integer(args$n_eval %||% 3000L)
batches <- as.integer(args$batches %||% 3L)
cores <- as.integer(args$cores %||% 1L)
seed <- as.integer(args$seed %||% 20260825L)
output_dir <- args$output_dir %||% file.path(
  "simulation_results", "section6_new_scenarios",
  "population_signal_normal_lg_plane_mixture_grid"
)

if (any(dimensions < 2L) || any(locations <= 0) ||
    any(correlations <= 0 | correlations >= 1) ||
    any(c(n_centres, n_profile, n_eval, batches, cores) < 1L)) {
  stop("Invalid population-signal scan controls.")
}

draw_plane_mixture <- function(family, d, n, location, correlation,
                                law = c("p0", "p1")) {
  law <- match.arg(law)
  e1 <- section6_e(d, 1L)
  Sigma_plus <- section6_sigma(d, "plus")
  Sigma_minus <- section6_sigma(d, "minus")
  Sigma_plus[1L, 2L] <- Sigma_plus[2L, 1L] <- correlation
  Sigma_minus[1L, 2L] <- Sigma_minus[2L, 1L] <- -correlation
  z <- mvtnorm::rmvnorm(n, mean = location * e1, sigma = Sigma_plus)
  if (identical(law, "p1")) {
    # This is the actual beta=1 law: P_1 = (P_0 + Q)/2.
    from_q <- stats::runif(n) < 0.5
    if (any(from_q)) {
      z[from_q, ] <- mvtnorm::rmvnorm(
        sum(from_q), mean = -location * e1, sigma = Sigma_minus
      )
    }
  }
  if (identical(family, "normal")) z else {
    logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L)
  }
}

one_batch <- function(job) {
  set.seed(job$seed)
  family <- as.character(job$family)
  d <- as.integer(job$d)
  a <- as.numeric(job$location)
  rho <- as.numeric(job$correlation)
  centres <- draw_plane_mixture(family, d, n_centres, a, rho, "p0")
  x0_a <- draw_plane_mixture(family, d, n_profile, a, rho, "p0")
  x1_a <- draw_plane_mixture(family, d, n_profile, a, rho, "p1")
  x0_b <- draw_plane_mixture(family, d, n_profile, a, rho, "p0")
  x1_b <- draw_plane_mixture(family, d, n_profile, a, rho, "p1")
  x0_eval <- draw_plane_mixture(family, d, n_eval, a, rho, "p0")
  spec <- population_signal_spec(family)
  control <- section6_control(1000L, job$seed, 50L)
  d0_a <- spec$distance_matrix(x0_a, centres, control)
  d1_a <- spec$distance_matrix(x1_a, centres, control)
  d0_b <- spec$distance_matrix(x0_b, centres, control)
  d1_b <- spec$distance_matrix(x1_b, centres, control)
  d0_eval <- spec$distance_matrix(x0_eval, centres, control)
  metrics <- t(vapply(seq_len(ncol(d0_a)), function(j) {
    centre_signal(d0_a[, j], d1_a[, j], d0_b[, j], d1_b[, j], d0_eval[, j])
  }, numeric(2)))
  cbind(job, data.frame(
    S_mean_center = mean(metrics[, "S"]),
    S_median_center = stats::median(metrics[, "S"]),
    J = mean(metrics[, "J"]),
    stringsAsFactors = FALSE
  ))
}

summarize <- function(x) {
  key <- interaction(x$family, x$d, x$location, x$correlation, drop = TRUE)
  output <- lapply(split(x, key), function(cell) {
    first <- cell[1L, c("family", "d", "location", "correlation"), drop = FALSE]
    first$null_covariance_condition_number <-
      (1 + first$correlation) / (1 - first$correlation)
    for (metric in c("S_mean_center", "S_median_center", "J")) {
      first[[metric]] <- mean(cell[[metric]])
      first[[paste0(metric, "_mcse")]] <- stats::sd(cell[[metric]]) / sqrt(nrow(cell))
    }
    first
  })
  out <- do.call(rbind, output)
  out[order(out$family, out$d, out$correlation, out$location), , drop = FALSE]
}

grid <- expand.grid(
  family = c("normal", "lg"), d = dimensions,
  location = locations, correlation = correlations,
  batch = seq_len(batches), KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$seed <- vapply(seq_len(nrow(grid)), function(i) {
  section6_seed(seed, i, grid$batch[[i]], stream = 251L)
}, integer(1))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_dir, "batch_results.csv")
summary_path <- file.path(output_dir, "summary.csv")
manifest_path <- file.path(output_dir, "manifest.txt")
if (!file.exists(manifest_path)) writeLines(c(
  "Known-parameter plane-mixture population-signal scan for Section 6 scenarios 1.1 and 2.1.",
  "P0 = N_d(a e1, I + rho(e1e2' + e2e1')); P1 = (P0 + Q)/2,",
  "Q = N_d(-a e1, I - rho(e1e2' + e2e1')).",
  "For logistic Gaussian, this construction is made in ILR coordinates.",
  "J is split-sample Monte Carlo for E[(F1,Omega(T)-F0,Omega(T))^2].",
  "No MLE, score derivative, Vhat, or bootstrap is used.",
  sprintf("dimensions: %s", paste(dimensions, collapse = ",")),
  sprintf("locations: %s", paste(locations, collapse = ",")),
  sprintf("correlations: %s", paste(correlations, collapse = ",")),
  sprintf("n_centres=%d; n_profile=%d; n_eval=%d; batches=%d; seed=%d",
    n_centres, n_profile, n_eval, batches, seed)
), manifest_path)

existing <- if (file.exists(raw_path)) utils::read.csv(raw_path, stringsAsFactors = FALSE) else data.frame()
key <- function(x) paste(x$family, x$d, x$location, x$correlation, x$batch, sep = "|")
pending <- grid[!key(grid) %in% if (nrow(existing)) key(existing) else character(), , drop = FALSE]
message(sprintf("Plane-mixture population scan: %d pending batches.", nrow(pending)))
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
