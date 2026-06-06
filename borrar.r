#!/usr/bin/env Rscript

# Diagnostic plots and KS-grid diagnostics for the cycle-24 sunspots
# symmetric two-small-circles fit.
#
# This script deliberately uses the fitted parameters already obtained in the
# completed KS run. It does not refit the model.

suppressPackageStartupMessages({
  library(grDevices)
  library(graphics)
})

repo_dir <- "/Users/Diego/Desktop/Codigo/goodness_of_fit_metric_spaces"
setwd(repo_dir)

source("utils.R")

output_dir <- file.path("real_data", "sunspots", "output", "cycle24_sphere_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

sunspots_path <- file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv")
if (!file.exists(sunspots_path)) {
  stop("Missing input file: ", sunspots_path, "\nRun real_data/sunspots/sunspots.R first.")
}

sunspots <- utils::read.csv(sunspots_path)
X_obs <- as.matrix(sunspots[, c("x1", "x2", "x3")])
X_obs <- X_obs / sqrt(rowSums(X_obs^2))
n <- nrow(X_obs)

# Fitted parameters from the completed cycle-24 KS run.
# Do not refit here; this script is only for quick visual diagnostics.
mu_hat <- c(0.002746, 0.002837, -0.999992)
mu_hat <- mu_hat / sqrt(sum(mu_hat^2))
kappa_hat <- 25.045825
nu_hat <- 0.251825

# For display only, orient the axis so that it points to the positive z hemisphere.
# The symmetric model is invariant under mu -> -mu.
mu_plot <- if (mu_hat[3L] < 0) -mu_hat else mu_hat

message(sprintf(
  "theta_hat: mu=(%.6f, %.6f, %.6f), kappa=%.6f, nu=%.6f",
  mu_hat[1L], mu_hat[2L], mu_hat[3L], kappa_hat, nu_hat
))

X_sim <- r_sph_small_circle_symmetric_mixture2(
  n = n,
  mu = mu_hat,
  kappa = kappa_hat,
  nu = nu_hat
)
X_sim <- X_sim / sqrt(rowSums(X_sim^2))

# Save the actual point clouds used in the diagnostic plots so they can be
# inspected later without regenerating the simulation.
observed_points_csv <- file.path(output_dir, "cycle24_observed_points.csv")
simulated_points_csv <- file.path(output_dir, "cycle24_simulated_points.csv")
utils::write.csv(
  data.frame(x1 = X_obs[, 1L], x2 = X_obs[, 2L], x3 = X_obs[, 3L]),
  observed_points_csv,
  row.names = FALSE
)
utils::write.csv(
  data.frame(x1 = X_sim[, 1L], x2 = X_sim[, 2L], x3 = X_sim[, 3L]),
  simulated_points_csv,
  row.names = FALSE
)

# Save the simulated sample with the same column structure as the sunspots input
# file. This can be used as a drop-in dataset for the GOF runner once the runner
# has an --input_csv option.
sunspots_sim <- sunspots
sunspots_sim$x1 <- X_sim[, 1L]
sunspots_sim$x2 <- X_sim[, 2L]
sunspots_sim$x3 <- X_sim[, 3L]
sunspots_sim$norm_s2 <- sqrt(rowSums(X_sim^2))
sunspots_sim$theta <- atan2(X_sim[, 2L], X_sim[, 1L])
sunspots_sim$phi <- asin(pmax(-1, pmin(1, X_sim[, 3L])))
sunspots_sim$hemisphere <- ifelse(
  sunspots_sim$phi > 0,
  "N",
  ifelse(sunspots_sim$phi < 0, "S", "Equator")
)
simulated_input_csv <- file.path(output_dir, "sunspots_cycle24_s2_simulated_from_fit.csv")
utils::write.csv(sunspots_sim, simulated_input_csv, row.names = FALSE)

# Orthographic sphere plotting using base R only.
plot_sphere_points <- function(X,
                               filename,
                               title,
                               subtitle = NULL,
                               mu_axis = NULL,
                               point_cex = 0.35,
                               alpha = 0.55) {
  png(filename, width = 1400, height = 1400, res = 180)
  on.exit(dev.off(), add = TRUE)

  op <- par(mar = c(1.5, 1.5, 3.5, 1.5), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)

  plot(
    NA,
    xlim = c(-1.08, 1.08),
    ylim = c(-1.08, 1.08),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title
  )
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.35, cex = 0.78)
  }

  ang <- seq(0, 2 * pi, length.out = 500L)
  lines(cos(ang), sin(ang), lwd = 1.3)

  lon_grid <- seq(0, 2 * pi, length.out = 300L)
  for (lat_deg in seq(-60, 60, by = 30)) {
    lat <- lat_deg * pi / 180
    x <- cos(lat) * cos(lon_grid)
    z <- rep(sin(lat), length(lon_grid))
    lines(x, z, lty = 3, col = "grey75")
  }
  lat_grid <- seq(-pi / 2, pi / 2, length.out = 300L)
  for (lon_deg in seq(0, 150, by = 30)) {
    lon <- lon_deg * pi / 180
    x <- cos(lat_grid) * cos(lon)
    z <- sin(lat_grid)
    lines(x, z, lty = 3, col = "grey85")
  }

  ord <- order(X[, 2L])
  Xo <- X[ord, , drop = FALSE]

  pal <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101L)
  idx <- pmin(101L, pmax(1L, floor((Xo[, 3L] + 1) / 2 * 100) + 1L))
  cols <- grDevices::adjustcolor(pal[idx], alpha.f = alpha)

  points(Xo[, 1L], Xo[, 3L], pch = 16, cex = point_cex, col = cols)

  if (!is.null(mu_axis)) {
    arrows(0, 0, mu_axis[1L], mu_axis[3L], length = 0.08, lwd = 2.5)
    arrows(0, 0, -mu_axis[1L], -mu_axis[3L], length = 0.08, lwd = 1.5, lty = 2)
  }

  legend(
    "bottomleft",
    legend = c("south", "equator", "north"),
    pch = 16,
    col = grDevices::adjustcolor(c("#2166AC", "#F7F7F7", "#B2182B"), alpha.f = 0.9),
    bty = "n",
    cex = 0.8
  )
}

subtitle <- sprintf(
  "n=%d | fitted kappa=%.3f, nu=%.3f | axis sign shown with mu_3 >= 0",
  n, kappa_hat, nu_hat
)

observed_png <- file.path(output_dir, "cycle24_observed_sphere.png")
simulated_png <- file.path(output_dir, "cycle24_simulated_symmetric_mixture_sphere.png")

plot_sphere_points(
  X = X_obs,
  filename = observed_png,
  title = "Observed sunspot births, cycle 24",
  subtitle = subtitle,
  mu_axis = mu_plot
)

plot_sphere_points(
  X = X_sim,
  filename = simulated_png,
  title = "Simulated from fitted symmetric small-circle mixture",
  subtitle = subtitle,
  mu_axis = mu_plot
)

utils::write.csv(
  data.frame(
    n = n,
    mu_1 = mu_hat[1L],
    mu_2 = mu_hat[2L],
    mu_3 = mu_hat[3L],
    kappa = kappa_hat,
    nu = nu_hat,
    observed_png = observed_png,
    simulated_png = simulated_png,
    observed_points_csv = observed_points_csv,
    simulated_points_csv = simulated_points_csv,
    simulated_input_csv = simulated_input_csv
  ),
  file.path(output_dir, "cycle24_sphere_plot_theta.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# KS-grid diagnostic: locate the grid point where the observed sample differs
# most from the fitted model, and compare it with one simulated sample.
# -------------------------------------------------------------------------

call_dp_sym <- function(omega, t_values) {
  fn <- get("distance_profile_small_circle_symmetric_mixture2", mode = "function")
  fml <- names(formals(fn))
  has_dots <- "..." %in% fml
  args <- list()
  if ("omega" %in% fml || has_dots) args$omega <- omega
  if ("t_values" %in% fml) {
    args$t_values <- t_values
  } else if ("t" %in% fml) {
    args$t <- t_values
  } else if ("t_grid" %in% fml) {
    args$t_grid <- t_values
  } else if (has_dots) {
    args$t_values <- t_values
  } else {
    stop("Cannot infer the t argument of distance_profile_small_circle_symmetric_mixture2().")
  }
  if ("mu" %in% fml || has_dots) args$mu <- mu_hat
  if ("kappa" %in% fml || has_dots) args$kappa <- kappa_hat
  if ("nu" %in% fml || has_dots) args$nu <- nu_hat
  if ("distance_type" %in% fml || has_dots) args$distance_type <- "geodesic"
  if ("method" %in% fml || has_dots) args$method <- "legendre"
  as.numeric(do.call(fn, args))
}

empirical_profile_from_t <- function(X, omega, t_values) {
  dot_values <- as.numeric(X %*% omega)
  thresholds <- cos(t_values)
  vapply(thresholds, function(a) mean(dot_values >= a), numeric(1L))
}

message("Computing KS-grid diagnostics...")
M_value <- 60L
ks_t_points <- 200L
omega_grid <- generate_canonical_lattice(M_value, dim = 3)
if (is.list(omega_grid) && !is.null(omega_grid$omega_grid)) {
  omega_grid <- omega_grid$omega_grid
}
omega_grid <- as.matrix(omega_grid)
omega_grid <- omega_grid / sqrt(rowSums(omega_grid^2))
t_grid <- seq(1e-8, pi - 1e-8, length.out = ks_t_points)

ks_rows <- vector("list", nrow(omega_grid))
for (i in seq_len(nrow(omega_grid))) {
  omega_i <- omega_grid[i, ]
  theoretical <- call_dp_sym(omega = omega_i, t_values = t_grid)
  empirical_obs <- empirical_profile_from_t(X_obs, omega_i, t_grid)
  empirical_sim <- empirical_profile_from_t(X_sim, omega_i, t_grid)
  diff_obs <- empirical_obs - theoretical
  diff_sim <- empirical_sim - theoretical
  ks_rows[[i]] <- data.frame(
    omega_id = i,
    omega_1 = omega_i[1L],
    omega_2 = omega_i[2L],
    omega_3 = omega_i[3L],
    t = t_grid,
    theoretical = theoretical,
    empirical_observed = empirical_obs,
    diff_observed = diff_obs,
    abs_diff_observed = abs(diff_obs),
    sqrt_n_abs_diff_observed = sqrt(n) * abs(diff_obs),
    empirical_simulated = empirical_sim,
    diff_simulated = diff_sim,
    abs_diff_simulated = abs(diff_sim),
    sqrt_n_abs_diff_simulated = sqrt(n) * abs(diff_sim)
  )
}
ks_diag <- do.call(rbind, ks_rows)
ks_diag <- ks_diag[order(-ks_diag$abs_diff_observed), ]
utils::write.csv(
  ks_diag,
  file.path(output_dir, "cycle24_ks_grid_diagnostics_sorted.csv"),
  row.names = FALSE
)

ks_summary <- data.frame(
  n = n,
  n_omega = nrow(omega_grid),
  n_t = length(t_grid),
  observed_ks_unscaled = max(ks_diag$abs_diff_observed, na.rm = TRUE),
  observed_ks_sqrt_n = sqrt(n) * max(ks_diag$abs_diff_observed, na.rm = TRUE),
  simulated_ks_unscaled = max(ks_diag$abs_diff_simulated, na.rm = TRUE),
  simulated_ks_sqrt_n = sqrt(n) * max(ks_diag$abs_diff_simulated, na.rm = TRUE)
)
utils::write.csv(
  ks_summary,
  file.path(output_dir, "cycle24_ks_grid_summary.csv"),
  row.names = FALSE
)

# Marginal diagnostics along the fitted axis and along solar latitude.
z_obs_fit <- as.numeric(X_obs %*% mu_hat)
z_sim_fit <- as.numeric(X_sim %*% mu_hat)
z_obs_solar <- X_obs[, 3L]
z_sim_solar <- X_sim[, 3L]

bin_breaks <- seq(-1, 1, length.out = 81L)
make_hist_df <- function(obs, sim, breaks, label) {
  h_obs <- hist(obs, breaks = breaks, plot = FALSE)
  h_sim <- hist(sim, breaks = breaks, plot = FALSE)
  data.frame(
    coordinate = label,
    bin_left = head(breaks, -1L),
    bin_right = tail(breaks, -1L),
    bin_mid = h_obs$mids,
    observed_count = h_obs$counts,
    simulated_count = h_sim$counts,
    observed_prop = h_obs$counts / sum(h_obs$counts),
    simulated_prop = h_sim$counts / sum(h_sim$counts),
    diff_prop = h_obs$counts / sum(h_obs$counts) - h_sim$counts / sum(h_sim$counts)
  )
}

marginal_diag <- rbind(
  make_hist_df(z_obs_fit, z_sim_fit, bin_breaks, "projection_on_fitted_axis"),
  make_hist_df(z_obs_solar, z_sim_solar, bin_breaks, "solar_z_coordinate")
)
utils::write.csv(
  marginal_diag,
  file.path(output_dir, "cycle24_marginal_hist_diagnostics.csv"),
  row.names = FALSE
)

# Simple cap diagnostic at the worst KS grid point.
worst <- ks_diag[1L, ]
omega_worst <- as.numeric(worst[c("omega_1", "omega_2", "omega_3")])
t_worst <- as.numeric(worst$t)
threshold_worst <- cos(t_worst)
cap_obs <- as.numeric(X_obs %*% omega_worst >= threshold_worst)
cap_sim <- as.numeric(X_sim %*% omega_worst >= threshold_worst)
utils::write.csv(
  data.frame(
    sample = c("observed", "simulated"),
    inside_count = c(sum(cap_obs), sum(cap_sim)),
    n = n,
    inside_prop = c(mean(cap_obs), mean(cap_sim)),
    theoretical = rep(as.numeric(worst$theoretical), 2L),
    omega_1 = rep(omega_worst[1L], 2L),
    omega_2 = rep(omega_worst[2L], 2L),
    omega_3 = rep(omega_worst[3L], 2L),
    t = rep(t_worst, 2L),
    dot_threshold = rep(threshold_worst, 2L)
  ),
  file.path(output_dir, "cycle24_worst_cap_diagnostic.csv"),
  row.names = FALSE
)

message("Saved observed plot to: ", observed_png)
message("Saved simulated plot to: ", simulated_png)
message("Saved point CSVs and diagnostics in: ", output_dir)
message("Saved simulated drop-in input CSV to: ", simulated_input_csv)
message(sprintf(
  "KS diagnostic | observed sqrt(n)*D=%.4f | simulated sqrt(n)*D=%.4f",
  ks_summary$observed_ks_sqrt_n,
  ks_summary$simulated_ks_sqrt_n
))
message("Next diagnostic: run the GOF runner on simulated_input_csv after adding an --input_csv option to the runner.")