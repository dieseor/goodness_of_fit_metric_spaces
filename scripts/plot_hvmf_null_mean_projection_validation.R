#!/usr/bin/env Rscript

# Diagnose the polar HvMF sampler used in the final first-HvMF-alternative
# experiment.  The experiment stores seeds and scalar results, rather than
# the simulated matrices.  This script reconstructs the null samples from
# those seeds and compares their mean-direction projections to the exact
# HvMF marginal law.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

repo_root <- normalizePath(".", mustWork = TRUE)
source(file.path(repo_root, "scripts", "run_hvmf_alt1_polar_power.R"))

results_dir <- file.path(
  repo_root,
  "simulation_results",
  "second_scenarios_power",
  "final_hvmf_alt1_polar_kappa2_M1000_B5000"
)
manifest_path <- file.path(results_dir, "manifest.csv")
output_dir <- file.path(results_dir, "diagnostics")
output_path <- file.path(output_dir, "hvmf_null_mean_projection_validation.png")

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
design <- manifest[manifest$n == 200L & manifest$beta == 0, , drop = FALSE]
if (nrow(design) != 1L) {
  stop("The manifest must contain exactly one null HvMF design with n = 200.")
}

mu0 <- as.numeric(strsplit(gsub(" ", "", design$null_mu), ";", fixed = TRUE)[[1L]])
if (length(mu0) != 3L || any(!is.finite(mu0))) {
  stop("Could not parse the null mean direction from the manifest.")
}

M <- as.integer(design$M)
n <- as.integer(design$n)
kappa <- as.numeric(design$kappa)
projection <- numeric(M * n)

for (rep_id in seq_len(M)) {
  set.seed(power_seed(design$base_seed, design$design_id, rep_id, stream = 0L))
  x <- generate_hvmf_alt1_polar_sample(n = n, beta = design$beta, kappa = kappa)
  index <- ((rep_id - 1L) * n + 1L):(rep_id * n)
  # T = -<X, mu0>_M = cosh(d_H(X, mu0)).
  projection[index] <- x[, 1L] * mu0[1L] - x[, 2L] * mu0[2L] - x[, 3L] * mu0[3L]
}

if (any(!is.finite(projection)) || any(projection < 1 - 1e-9)) {
  stop("The reconstructed projections are invalid.")
}
projection <- pmax(projection, 1)

# At omega = mu0 on H^2, T - 1 has exactly an Exp(kappa) distribution.
ks <- stats::ks.test(projection - 1, "pexp", rate = kappa)
y_upper <- 1 + max(stats::qexp(0.9995, rate = kappa), stats::quantile(projection - 1, 0.9995))
y_grid <- seq(1, y_upper, length.out = 1000L)
exact_density <- hvmf_projection_density_h2(y_grid, alpha = 1, kappa = kappa)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
grDevices::png(output_path, width = 1800L, height = 1200L, res = 180L)
graphics::par(mar = c(4.4, 4.8, 3.2, 1.1), las = 1)
graphics::hist(
  projection,
  breaks = 100L,
  probability = TRUE,
  xlim = c(1, y_upper),
  border = "white",
  col = "#A6CEE3",
  main = "HvMF polar sampler: projection onto the null mean direction",
  xlab = "T = -<X, mu_0>_M",
  ylab = "Density"
)
graphics::lines(y_grid, exact_density, col = "#D7301F", lwd = 3L)
graphics::legend(
  "topright",
  legend = c(
    "Reconstructed final null samples",
    sprintf("Exact HvMF density: %g exp[-%g (t - 1)]", kappa, kappa)
  ),
  fill = c("#A6CEE3", NA),
  border = c("white", NA),
  lty = c(NA, 1L),
  lwd = c(NA, 3L),
  col = c(NA, "#D7301F"),
  bty = "n"
)
graphics::mtext(
  sprintf(
    "M = %d replications, n = %d, kappa = %g;  KS distance = %.5f;  mean(T - 1) = %.5f (exact: %.5f)",
    M, n, kappa, unname(ks$statistic), mean(projection - 1), 1 / kappa
  ),
  side = 3L,
  line = 0.2,
  cex = 0.8
)
grDevices::dev.off()

cat(sprintf("Wrote %s\n", output_path))
cat(sprintf("KS distance: %.8f; mean(T - 1): %.8f; exact mean: %.8f\n", unname(ks$statistic), mean(projection - 1), 1 / kappa))
