#!/usr/bin/env Rscript

# Reproducible validation of the regularised polar HvMF sampler in the four
# Section 6 null settings. The radial quantile table is regularised at
# p_max = 0.999, while the reference CDF below is the untruncated HvMF law.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(".", mustWork = TRUE)
source(file.path(repo_root, "utils.R"))

output_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  file.path(repo_root, "simulation_results", "hvmf_general_sampler_validation")
}
M <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
base_seed <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20260727L

if (!is.finite(M) || M < 1L) stop("`M` must be a strictly positive integer.")
if (!is.finite(base_seed)) stop("`base_seed` must be a finite integer.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p_max <- 0.999
distance_tolerance <- 0.015
n_values <- c(50L, 100L, 200L)
configurations <- list(
  list(q = 2L, scenario = "mixture", kappa = 2),
  list(q = 2L, scenario = "angular", kappa = 3),
  list(q = 10L, scenario = "mixture", kappa = 10),
  list(q = 10L, scenario = "angular", kappa = 15)
)

section6_mu <- function(q) {
  c(sqrt(2), 1, rep.int(0, q - 1L))
}

projection_onto_mu <- function(x, mu) {
  x[, 1L] * mu[[1L]] -
    rowSums(x[, -1L, drop = FALSE] * rep(mu[-1L], each = nrow(x)))
}

trapezoidal_cdf_table <- function(density,
                                  direct_cdf,
                                  lower,
                                  initial_upper,
                                  grid_size = 32769L,
                                  tail_tolerance = 1e-10) {
  upper <- initial_upper
  upper_mass <- direct_cdf(upper)
  while (upper_mass < 1 - tail_tolerance) {
    upper <- 2 * upper
    if (!is.finite(upper) || upper > 256) {
      stop("Could not find a finite CDF tabulation endpoint.")
    }
    upper_mass <- direct_cdf(upper)
  }

  grid <- seq(lower, upper, length.out = grid_size)
  density_values <- density(grid)
  increments <- diff(grid) *
    (density_values[-1L] + density_values[-length(grid)]) / 2
  cdf <- c(0, cumsum(increments))
  cdf <- pmin(pmax(cdf, 0), 1)

  check_grid <- seq(lower, upper, length.out = 31L)
  direct_values <- direct_cdf(check_grid)
  interpolated_values <- stats::approx(
    x = grid,
    y = cdf,
    xout = check_grid,
    rule = 2
  )$y
  approximation_error <- max(abs(direct_values - interpolated_values))

  list(
    grid = grid,
    density = density_values,
    cdf = cdf,
    upper = upper,
    upper_mass = upper_mass,
    approximation_error = approximation_error
  )
}

cdf_evaluator <- function(table) {
  force(table)
  function(x) {
    x <- as.numeric(x)
    output <- numeric(length(x))
    output[is.infinite(x) & x > 0] <- 1
    above <- is.finite(x) & x >= table$grid[[1L]]
    if (any(above)) {
      output[above] <- stats::approx(
        x = table$grid,
        y = table$cdf,
        xout = x[above],
        rule = 2
      )$y
      output[above & x >= table$upper] <- 1
    }
    pmin(pmax(output, 0), 1)
  }
}

cdf_quantile <- function(probability, table) {
  keep <- !duplicated(table$cdf)
  stats::approx(
    x = table$cdf[keep],
    y = table$grid[keep],
    xout = probability,
    rule = 2
  )$y
}

regularised_population_distance <- function(q, kappa, chi = asinh(1)) {
  quantile_table <- hvmf_build_radial_quantile_table(
    q = q,
    kappa = kappa,
    chi = chi,
    p_max = p_max,
    probability_step = 0.01,
    upper = 3
  )
  probability <- seq(0, p_max, length.out = 5001L)
  quantile <- stats::approx(
    x = quantile_table$base,
    y = quantile_table$quantile_values,
    xout = probability
  )$y
  true_cdf <- hvmf_radial_cdf(
    u = quantile,
    q = q,
    kappa = kappa,
    chi = chi
  )

  max(abs(probability / p_max - true_cdf), 1 - p_max)
}

one_sample_ks <- function(sample, cdf) {
  suppressWarnings(stats::ks.test(sample, cdf, exact = FALSE))
}

n_cells <- length(configurations) * length(n_values)
cell_alpha <- 0.05 / n_cells
calibration_lower <- stats::qbinom(cell_alpha / 2, size = M, prob = 0.05) / M
calibration_upper <- stats::qbinom(1 - cell_alpha / 2, size = M, prob = 0.05) / M

summary_records <- list()
population_records <- list()
plot_records <- list()
record_index <- 1L

for (configuration_index in seq_along(configurations)) {
  configuration <- configurations[[configuration_index]]
  q <- configuration$q
  kappa <- configuration$kappa
  scenario <- configuration$scenario
  chi <- asinh(1)
  mu <- section6_mu(q)

  cat(sprintf(
    "Preparing q = %d, scenario = %s, kappa = %g\n",
    q, scenario, kappa
  ))

  radial_density <- function(u) {
    hvmf_radial_density(u = u, q = q, kappa = kappa, chi = chi)
  }
  radial_direct_cdf <- function(u) {
    hvmf_radial_cdf(u = u, q = q, kappa = kappa, chi = chi)
  }
  radial_table <- trapezoidal_cdf_table(
    density = radial_density,
    direct_cdf = radial_direct_cdf,
    lower = 0,
    initial_upper = 3
  )
  radial_cdf <- cdf_evaluator(radial_table)

  projection_density <- function(y) {
    hvmf_mean_projection_density(y = y, q = q, kappa = kappa)
  }
  projection_direct_cdf <- function(y) {
    hvmf_mean_projection_cdf(y = y, q = q, kappa = kappa)
  }
  projection_table <- trapezoidal_cdf_table(
    density = projection_density,
    direct_cdf = projection_direct_cdf,
    lower = 1,
    initial_upper = 4
  )
  projection_cdf <- cdf_evaluator(projection_table)

  population_distance <- regularised_population_distance(q = q, kappa = kappa)
  population_records[[configuration_index]] <- data.frame(
    q = q,
    scenario = scenario,
    kappa = kappa,
    p_max = p_max,
    population_KS_distance_to_true_radial = population_distance,
    distance_tolerance = distance_tolerance,
    population_distance_pass = population_distance <= distance_tolerance,
    radial_cdf_tabulation_error = radial_table$approximation_error,
    projection_cdf_tabulation_error = projection_table$approximation_error
  )

  set.seed(base_seed + 10000L * q + 100L * configuration_index)
  n_max <- max(n_values)
  x <- rhvmf_polar(
    n = M * n_max,
    mu = mu,
    kappa = kappa,
    p_max = p_max
  )
  radial_matrix <- matrix(acosh(x[, 1L]), nrow = M, ncol = n_max, byrow = TRUE)
  projection_matrix <- matrix(
    pmax(projection_onto_mu(x, mu), 1),
    nrow = M,
    ncol = n_max,
    byrow = TRUE
  )

  exact_projection_mean <- stats::integrate(
    function(y) y * projection_density(y),
    lower = 1,
    upper = Inf,
    rel.tol = 1e-8
  )$value
  projection_tail_cutoff <- cdf_quantile(0.95, projection_table)

  for (n in n_values) {
    radial_samples <- radial_matrix[, seq_len(n), drop = FALSE]
    projection_samples <- projection_matrix[, seq_len(n), drop = FALSE]
    radial_pooled <- as.vector(t(radial_samples))
    projection_pooled <- as.vector(t(projection_samples))

    radial_gof <- lapply(seq_len(M), function(replication) {
      one_sample_ks(radial_samples[replication, ], radial_cdf)
    })
    radial_p_values <- vapply(radial_gof, function(result) result$p.value, numeric(1))
    radial_distances <- vapply(
      radial_gof,
      function(result) unname(result$statistic),
      numeric(1)
    )

    radial_pooled_ks <- one_sample_ks(radial_pooled, radial_cdf)
    projection_pooled_ks <- one_sample_ks(projection_pooled, projection_cdf)
    rejection_rate <- mean(radial_p_values < 0.05)
    N <- length(radial_pooled)
    empirical_distance_limit <-
      distance_tolerance + sqrt(log(200) / (2 * N))

    uniformity_test <- suppressWarnings(stats::ks.test(
      radial_p_values,
      "punif",
      exact = FALSE
    ))

    summary_records[[record_index]] <- data.frame(
      q = q,
      scenario = scenario,
      kappa = kappa,
      n = n,
      M = M,
      N = N,
      radial_rejection_rate_0_05 = rejection_rate,
      simultaneous_binomial_lower = calibration_lower,
      simultaneous_binomial_upper = calibration_upper,
      calibration_pass =
        rejection_rate >= calibration_lower &&
        rejection_rate <= calibration_upper,
      radial_pvalue_uniformity_p = uniformity_test$p.value,
      radial_pvalue_mean = mean(radial_p_values),
      radial_pvalue_variance = stats::var(radial_p_values),
      radial_pvalue_below_0_01 = mean(radial_p_values < 0.01),
      radial_pvalue_below_0_05 = rejection_rate,
      radial_pvalue_below_0_10 = mean(radial_p_values < 0.10),
      mean_replication_radial_KS = mean(radial_distances),
      pooled_radial_KS = unname(radial_pooled_ks$statistic),
      pooled_projection_KS = unname(projection_pooled_ks$statistic),
      empirical_distance_limit = empirical_distance_limit,
      radial_distance_pass =
        unname(radial_pooled_ks$statistic) <= empirical_distance_limit,
      projection_distance_pass =
        unname(projection_pooled_ks$statistic) <= empirical_distance_limit,
      empirical_projection_mean = mean(projection_pooled),
      exact_projection_mean = exact_projection_mean,
      empirical_projection_tail_0_05 =
        mean(projection_pooled > projection_tail_cutoff),
      exact_projection_tail = 0.05
    )

    plot_records[[record_index]] <- list(
      q = q,
      scenario = scenario,
      kappa = kappa,
      n = n,
      M = M,
      projection = projection_pooled,
      radial_p_values = radial_p_values,
      projection_KS = unname(projection_pooled_ks$statistic),
      projection_mean = mean(projection_pooled),
      exact_projection_mean = exact_projection_mean,
      projection_table = projection_table
    )

    cat(sprintf(
      "  n = %d: rejection = %.4f, radial D = %.5f, projection D = %.5f\n",
      n,
      rejection_rate,
      unname(radial_pooled_ks$statistic),
      unname(projection_pooled_ks$statistic)
    ))
    record_index <- record_index + 1L
  }
}

summary_table <- do.call(rbind, summary_records)
population_table <- do.call(rbind, population_records)

utils::write.csv(
  summary_table,
  file.path(output_dir, "hvmf_general_sampler_validation_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  population_table,
  file.path(output_dir, "hvmf_general_sampler_population_distance.csv"),
  row.names = FALSE
)

plot_projection_panel <- function(record, add_axis_labels = TRUE) {
  x_max <- cdf_quantile(0.9995, record$projection_table)
  observed_max <- max(record$projection)
  breaks <- seq(1, max(x_max, observed_max), length.out = 101L)
  y_grid <- seq(1, x_max, length.out = 1000L)
  y_density <- hvmf_mean_projection_density(
    y_grid,
    q = record$q,
    kappa = record$kappa
  )
  y_limit <- max(y_density) * 1.08

  graphics::hist(
    record$projection,
    breaks = breaks,
    probability = TRUE,
    xlim = c(1, x_max),
    ylim = c(0, y_limit),
    border = "white",
    col = "#A6CEE3",
    main = sprintf("%s, n = %d", record$scenario, record$n),
    xlab = if (add_axis_labels) "T = -<X, mu>_M" else "",
    ylab = if (add_axis_labels) "Density" else ""
  )
  graphics::lines(y_grid, y_density, col = "#D7301F", lwd = 2.5)
  graphics::mtext(
    sprintf(
      "N=%d; D=%.5f; mean=%.5f (exact %.5f)",
      length(record$projection),
      record$projection_KS,
      record$projection_mean,
      record$exact_projection_mean
    ),
    side = 3,
    line = 0.1,
    cex = 0.7
  )
}

for (record in plot_records) {
  output_path <- file.path(
    output_dir,
    sprintf(
      "projection_q%d_%s_n%d.png",
      record$q,
      record$scenario,
      record$n
    )
  )
  grDevices::png(output_path, width = 1800L, height = 1200L, res = 180L)
  graphics::par(mar = c(4.3, 4.5, 3.1, 0.8), las = 1)
  plot_projection_panel(record)
  graphics::legend(
    "topright",
    legend = "True untruncated HvMF density",
    col = "#D7301F",
    lty = 1,
    lwd = 2.5,
    bty = "n"
  )
  grDevices::dev.off()
}

for (q in c(2L, 10L)) {
  records_q <- Filter(function(record) record$q == q, plot_records)

  grDevices::png(
    file.path(output_dir, sprintf("projection_validation_q%d_combined.png", q)),
    width = 2100L,
    height = 1300L,
    res = 180L
  )
  graphics::par(
    mfrow = c(2L, 3L),
    mar = c(3.5, 3.7, 2.7, 0.5),
    oma = c(4, 5, 3.3, 0),
    las = 1
  )
  for (record in records_q) {
    plot_projection_panel(record, add_axis_labels = FALSE)
  }
  graphics::mtext(
    "T = -<X, mu>_M",
    side = 1,
    outer = TRUE,
    line = 2.1
  )
  graphics::mtext(
    "Density",
    side = 2,
    outer = TRUE,
    line = 2.5,
    las = 0
  )
  graphics::mtext(
    sprintf(
      "HvMF sampler on H^%d: projections versus the true untruncated law",
      q
    ),
    side = 3,
    outer = TRUE,
    line = 1.2,
    cex = 1.15,
    font = 2
  )
  grDevices::dev.off()

  grDevices::png(
    file.path(output_dir, sprintf("radial_ks_pvalues_q%d_combined.png", q)),
    width = 2100L,
    height = 1300L,
    res = 180L
  )
  graphics::par(
    mfrow = c(2L, 3L),
    mar = c(3.5, 3.7, 2.7, 0.5),
    oma = c(4, 5, 3.3, 0),
    las = 1
  )
  for (record in records_q) {
    graphics::hist(
      record$radial_p_values,
      breaks = seq(0, 1, by = 0.05),
      probability = TRUE,
      xlim = c(0, 1),
      ylim = c(0, 1.5),
      border = "white",
      col = "#B2DF8A",
      main = sprintf("%s, n = %d", record$scenario, record$n),
      xlab = "",
      ylab = ""
    )
    graphics::abline(h = 1, col = "#54278F", lwd = 2.5)
  }
  graphics::mtext(
    "Radial KS p-value",
    side = 1,
    outer = TRUE,
    line = 2.1
  )
  graphics::mtext(
    "Density",
    side = 2,
    outer = TRUE,
    line = 2.5,
    las = 0
  )
  graphics::mtext(
    sprintf(
      "HvMF sampler on H^%d: radial GOF p-values against the true law",
      q
    ),
    side = 3,
    outer = TRUE,
    line = 1.2,
    cex = 1.15,
    font = 2
  )
  grDevices::dev.off()
}

all_pass <-
  all(population_table$population_distance_pass) &&
  all(summary_table$calibration_pass) &&
  all(summary_table$radial_distance_pass) &&
  all(summary_table$projection_distance_pass)

cat(sprintf("Wrote validation outputs to %s\n", output_dir))
cat(sprintf("Overall acceptance: %s\n", if (all_pass) "PASS" else "FAIL"))

if (!all_pass) {
  stop("One or more HvMF sampler validation criteria failed.")
}
