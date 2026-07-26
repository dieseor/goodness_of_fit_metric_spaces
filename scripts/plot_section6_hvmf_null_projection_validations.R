#!/usr/bin/env Rscript

# Projection diagnostics for every pure-HvMF data-generating setting used in
# Section 6: the two HvMF scenarios, each at n = 50, 100, and 200.  The
# simulations persist their design, seeds, and scalar results but not X.  We
# therefore reconstruct each null sample with exactly the generator used by
# its corresponding experiment, without rerunning its bootstrap.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

repo_root <- normalizePath(".", mustWork = TRUE)
source(file.path(repo_root, "scripts", "run_hvmf_alt1_polar_power.R"))

section6_root <- file.path(repo_root, "simulation_results", "second_scenarios_power")
output_dir <- file.path(section6_root, "section6_hvmf_null_projection_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scenario_specs <- list(
  list(
    id = "hvmf_mixture",
    label = "HvMF mixture scenario",
    results_dir = file.path(section6_root, "final_hvmf_alt1_polar_kappa2_M1000_B5000"),
    generate = function(design) {
      generate_hvmf_alt1_polar_sample(n = design$n, beta = design$beta, kappa = design$kappa)
    },
    mean_direction = function(design) {
      as.numeric(strsplit(gsub(" ", "", design$null_mu), ";", fixed = TRUE)[[1L]])
    }
  ),
  list(
    id = "hvmf_angular",
    label = "HvMF angular scenario",
    results_dir = file.path(section6_root, "final_hvmf_angular_kappa2_delta_pi5_M1000_B5000"),
    generate = function(design) {
      candidate <- candidate_from_label("hvmf", as.character(design$candidate))
      generate_power_sample("hvmf", candidate = candidate, n = design$n, beta = design$beta)
    },
    mean_direction = function(design) power_mu_hvmf
  )
)

reconstruct_projection <- function(spec, design) {
  M <- as.integer(design$M)
  n <- as.integer(design$n)
  mu0 <- as.numeric(spec$mean_direction(design))
  if (length(mu0) != 3L || any(!is.finite(mu0))) {
    stop(sprintf("Invalid null mean direction for %s.", spec$id))
  }

  projection <- numeric(M * n)
  for (rep_id in seq_len(M)) {
    set.seed(power_seed(design$base_seed, design$design_id, rep_id, stream = 0L))
    x <- spec$generate(design)
    index <- ((rep_id - 1L) * n + 1L):(rep_id * n)
    # T = -<X, mu0>_M = cosh(d_H(X, mu0)).
    projection[index] <- x[, 1L] * mu0[1L] - x[, 2L] * mu0[2L] - x[, 3L] * mu0[3L]
  }
  if (any(!is.finite(projection)) || any(projection < 1 - 1e-9)) {
    stop(sprintf("Invalid reconstructed projections for %s, n = %d.", spec$id, n))
  }
  pmax(projection, 1)
}

plot_projection_histogram <- function(projection, kappa, title, subtitle, output_path, breaks, xlim, ylim) {
  y_grid <- seq(xlim[[1L]], xlim[[2L]], length.out = 1000L)
  exact_density <- hvmf_projection_density_h2(y_grid, alpha = 1, kappa = kappa)
  grDevices::png(output_path, width = 1800L, height = 1200L, res = 180L)
  graphics::par(mar = c(4.4, 4.8, 3.2, 1.1), las = 1)
  graphics::hist(
    projection,
    breaks = breaks,
    probability = TRUE,
    xlim = xlim,
    ylim = ylim,
    border = "white",
    col = "#A6CEE3",
    main = title,
    sub = subtitle,
    xlab = "T = -<X, mu_0>_M",
    ylab = "Density"
  )
  graphics::lines(y_grid, exact_density, col = "#D7301F", lwd = 3L)
  graphics::legend(
    "topright",
    legend = c("Reconstructed null samples", sprintf("Exact density: %g exp[-%g (t - 1)]", kappa, kappa)),
    fill = c("#A6CEE3", NA), border = c("white", NA),
    lty = c(NA, 1L), lwd = c(NA, 3L), col = c(NA, "#D7301F"), bty = "n"
  )
  grDevices::dev.off()
}

records <- list()
panel_index <- 1L
for (spec in scenario_specs) {
  manifest <- utils::read.csv(file.path(spec$results_dir, "manifest.csv"), stringsAsFactors = FALSE)
  designs <- manifest[manifest$beta == 0 & manifest$n %in% c(50L, 100L, 200L), , drop = FALSE]
  designs <- designs[order(designs$n), , drop = FALSE]
  if (nrow(designs) != 3L || anyDuplicated(designs$n)) {
    stop(sprintf("Expected one pure-HvMF design at n = 50, 100, 200 for %s.", spec$id))
  }

  for (row_index in seq_len(nrow(designs))) {
    design <- designs[row_index, , drop = FALSE]
    projection <- reconstruct_projection(spec, design)
    kappa <- if ("kappa" %in% names(design)) as.numeric(design$kappa) else candidate_from_label("hvmf", design$candidate)$kappa
    ks <- stats::ks.test(projection - 1, "pexp", rate = kappa)
    records[[panel_index]] <- list(
      scenario = spec$id,
      scenario_label = spec$label,
      n = as.integer(design$n),
      M = as.integer(design$M),
      kappa = kappa,
      projection = projection,
      ks_distance = unname(ks$statistic),
      mean_shifted_projection = mean(projection - 1),
      tail_probability = mean(projection - 1 > 1.5)
    )
    panel_index <- panel_index + 1L
  }
}

# All six null settings have kappa = 2, so a common scale makes the sampling
# differences and the common sampler discrepancy directly comparable.
kappa_values <- vapply(records, `[[`, numeric(1), "kappa")
if (length(unique(kappa_values)) != 1L) stop("Section 6 null settings must have a common kappa for this plot.")
kappa <- kappa_values[[1L]]
xlim <- c(1, 1 + stats::qexp(0.9995, rate = kappa))
histogram_max <- max(vapply(records, function(record) max(record$projection), numeric(1)))
breaks <- seq(1, max(xlim[[2L]], histogram_max), length.out = 101L)
ylim <- c(0, 2.1)

summary_table <- do.call(rbind, lapply(records, function(record) {
  data.frame(
    scenario = record$scenario,
    n = record$n,
    M = record$M,
    N = length(record$projection),
    kappa = record$kappa,
    mean_T_minus_1 = record$mean_shifted_projection,
    exact_mean_T_minus_1 = 1 / record$kappa,
    KS_distance = record$ks_distance,
    empirical_P_T_minus_1_gt_1_5 = record$tail_probability,
    exact_P_T_minus_1_gt_1_5 = exp(-1.5 * record$kappa)
  )
}))
utils::write.csv(summary_table, file.path(output_dir, "section6_hvmf_null_projection_summary.csv"), row.names = FALSE)

for (record in records) {
  output_path <- file.path(output_dir, sprintf("%s_n%d.png", record$scenario, record$n))
  subtitle <- sprintf(
    "M = %d, n = %d, kappa = %g; KS distance = %.5f; mean(T - 1) = %.5f (exact: %.5f)",
    record$M, record$n, record$kappa, record$ks_distance,
    record$mean_shifted_projection, 1 / record$kappa
  )
  plot_projection_histogram(
    record$projection, record$kappa, record$scenario_label, subtitle,
    output_path, breaks = breaks, xlim = xlim, ylim = ylim
  )
}

combined_path <- file.path(output_dir, "section6_hvmf_null_projection_validation_combined.png")
grDevices::png(combined_path, width = 2100L, height = 1300L, res = 180L)
graphics::par(mfrow = c(2L, 3L), mar = c(3.6, 3.8, 2.7, 0.6), las = 1, oma = c(0, 0, 3.3, 0))
y_grid <- seq(xlim[[1L]], xlim[[2L]], length.out = 1000L)
exact_density <- hvmf_projection_density_h2(y_grid, alpha = 1, kappa = kappa)
for (record in records) {
  graphics::hist(
    record$projection, breaks = breaks, probability = TRUE, xlim = xlim, ylim = ylim,
    border = "white", col = "#A6CEE3",
    main = sprintf("%s, n = %d", record$scenario_label, record$n),
    xlab = "T = -<X, mu_0>_M", ylab = "Density"
  )
  graphics::lines(y_grid, exact_density, col = "#D7301F", lwd = 2.5)
  graphics::mtext(
    sprintf("N = %d; D = %.5f; mean = %.5f", length(record$projection), record$ks_distance, record$mean_shifted_projection),
    side = 3L, line = 0.15, cex = 0.7
  )
}
graphics::mtext(
  "Section 6: mean-direction projections under the pure HvMF null (red: exact density)",
  side = 3L, outer = TRUE, line = 1.3, cex = 1.15, font = 2
)
grDevices::dev.off()

cat(sprintf("Wrote six individual figures, the combined figure, and the summary to %s\n", output_dir))
