#!/usr/bin/env Rscript

# Conditional small-circle-mixture analysis for the central part of cycle 23.
# This runner deliberately stops before the GOF/bootstrap stage.

resolve_sunspots_time_varying_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

utils_path_sunspots_time_varying <- resolve_sunspots_time_varying_path("utils.R")
prep_path_sunspots_time_varying <- resolve_sunspots_time_varying_path(
  "real_data", "sunspots", "sunspots.R"
)
source(utils_path_sunspots_time_varying)

sunspots_time_varying_softplus <- function(x) {
  x <- as.numeric(x)
  pmax(x, 0) + log1p(exp(-abs(x)))
}

sunspots_time_varying_inverse_softplus <- function(x) {
  x <- as.numeric(x)
  ifelse(x > 30, x, log(expm1(x)))
}

sunspots_time_varying_normalize_hemisphere_regression <- function(hemisphere_regression = "asymmetric") {
  hemisphere_regression <- tolower(as.character(hemisphere_regression))
  if (length(hemisphere_regression) != 1L ||
      !hemisphere_regression %in% c("asymmetric", "shared")) {
    stop("`hemisphere_regression` must be either 'asymmetric' or 'shared'.")
  }
  hemisphere_regression
}

sunspots_time_varying_validate_theta <- function(theta,
                                                  nu_eps = 1e-6,
                                                  c_min = 1e-8,
                                                  c_max = 1e6,
                                                  hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  if (identical(hemisphere_regression, "shared") &&
      all(c("a", "b", "c") %in% names(theta)) &&
      !all(c("a_N", "b_N", "a_S", "b_S") %in% names(theta))) {
    theta <- c(list(a_N = theta$a, b_N = theta$b, a_S = theta$a, b_S = theta$b), theta["c"])
  }
  required <- c("a_N", "b_N", "a_S", "b_S", "c")
  if (!is.list(theta) || !all(required %in% names(theta))) {
    stop("`theta` must contain a_N, b_N, a_S, b_S, and c.")
  }

  values <- vapply(theta[required], as.numeric, numeric(1L))
  if (any(!is.finite(values))) {
    stop("All time-varying mixture parameters must be finite scalars.")
  }
  if (!is.finite(nu_eps) || nu_eps <= 0 || nu_eps >= 0.5) {
    stop("`nu_eps` must lie in (0, 0.5).")
  }
  if (!is.finite(c_min) || !is.finite(c_max) || c_min <= 0 || c_max <= c_min) {
    stop("`c_min` and `c_max` must satisfy 0 < c_min < c_max.")
  }

  nu_N_end <- values[["a_N"]] + values[["b_N"]]
  nu_S_end <- values[["a_S"]] + values[["b_S"]]
  if (values[["b_N"]] >= 0 || values[["b_S"]] >= 0) {
    stop("Both latitude regressions must be strictly decreasing: b_N < 0 and b_S < 0.")
  }
  if (identical(hemisphere_regression, "shared") &&
      (!isTRUE(all.equal(values[["a_N"]], values[["a_S"]], tolerance = 0)) ||
       !isTRUE(all.equal(values[["b_N"]], values[["b_S"]], tolerance = 0)))) {
    stop("The shared hemisphere regression requires a_N = a_S and b_N = b_S.")
  }
  if (any(c(values[["a_N"]], nu_N_end, values[["a_S"]], nu_S_end) <= nu_eps) ||
      any(c(values[["a_N"]], nu_N_end, values[["a_S"]], nu_S_end) >= 1 - nu_eps)) {
    stop("Both linear latitude paths must remain in (nu_eps, 1 - nu_eps) on [0, 1].")
  }
  if (values[["c"]] < c_min || values[["c"]] > c_max) {
    stop("`c` is outside the configured admissible interval.")
  }

  c(values, list(
    nu_N_end = nu_N_end,
    nu_S_end = nu_S_end,
    hemisphere_regression = hemisphere_regression,
    n_parameters = if (identical(hemisphere_regression, "shared")) 3L else 5L
  ))
}

sunspots_time_varying_unpack_par <- function(par,
                                             nu_eps = 1e-6,
                                             c_min = 1e-8,
                                             c_max = 1e6,
                                             hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  par <- as.numeric(par)
  expected_length <- if (identical(hemisphere_regression, "shared")) 3L else 5L
  if (length(par) != expected_length || any(!is.finite(par))) {
    stop(sprintf("The unconstrained parameter vector must contain %d finite values.", expected_length))
  }

  nu_upper <- 1 - nu_eps
  decode_path <- function(start_raw, end_raw) {
    # BFGS can evaluate very large finite values, for which plogis rounds to 0 or 1.
    prob_eps <- sqrt(.Machine$double.eps)
    start_prob <- pmin(pmax(stats::plogis(start_raw), prob_eps), 1 - prob_eps)
    end_prob <- pmin(pmax(stats::plogis(end_raw), prob_eps), 1 - prob_eps)
    start <- nu_eps + (nu_upper - nu_eps) * start_prob
    end <- nu_eps + (start - nu_eps) * end_prob
    c(a = start, b = end - start)
  }

  north <- decode_path(par[[1L]], par[[2L]])
  south <- if (identical(hemisphere_regression, "shared")) north else decode_path(par[[3L]], par[[4L]])
  c_index <- if (identical(hemisphere_regression, "shared")) 3L else 5L
  theta <- list(
    a_N = north[["a"]],
    b_N = north[["b"]],
    a_S = south[["a"]],
    b_S = south[["b"]],
    c = min(max(sunspots_time_varying_softplus(par[[c_index]]), c_min), c_max)
  )
  sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
}

sunspots_time_varying_pack_theta <- function(theta,
                                             nu_eps = 1e-6,
                                             c_min = 1e-8,
                                             c_max = 1e6,
                                             hemisphere_regression = "asymmetric") {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  theta <- sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  encode_path <- function(start, end) {
    start_prob <- (start - nu_eps) / (1 - 2 * nu_eps)
    end_prob <- (end - nu_eps) / (start - nu_eps)
    c(stats::qlogis(start_prob), stats::qlogis(end_prob))
  }
  if (identical(hemisphere_regression, "shared")) {
    return(c(
      encode_path(theta$a_N, theta$nu_N_end),
      sunspots_time_varying_inverse_softplus(theta$c)
    ))
  }
  c(encode_path(theta$a_N, theta$nu_N_end), encode_path(theta$a_S, theta$nu_S_end),
    sunspots_time_varying_inverse_softplus(theta$c))
}

sunspots_time_varying_nu <- function(u, theta) {
  u <- as.numeric(u)
  if (any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) {
    stop("`u` must contain finite values in (0, 1).")
  }
  theta <- sunspots_time_varying_validate_theta(theta)
  list(
    north = theta$a_N + theta$b_N * u,
    south = theta$a_S + theta$b_S * u
  )
}

# The normalization is evaluated by the existing vectorized small-circle implementation.
sunspots_time_varying_component_log_axis_density <- function(z, c_value, nu) {
  z <- as.numeric(z)
  nu <- as.numeric(nu)
  if (length(z) != length(nu)) {
    stop("`z` and `nu` must have the same length.")
  }
  log_norm <- small_circle_log_norm_constant(kappa = c_value, nu = nu)
  -log(2) - log_norm - c_value * (z - nu)^2
}

sunspots_time_varying_log_density <- function(x, u, theta) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  u <- as.numeric(u)
  if (length(u) != nrow(x)) {
    stop("`u` must have one entry per row of `x`.")
  }
  theta <- sunspots_time_varying_validate_theta(theta)
  nu <- sunspots_time_varying_nu(u, theta)
  z <- pmin(pmax(x[, 3L], -1), 1)

  log_north <- sunspots_time_varying_component_log_axis_density(z, theta$c, nu$north)
  log_south <- sunspots_time_varying_component_log_axis_density(-z, theta$c, nu$south)
  -log(2 * pi) + rotational_logsumexp2(log(0.5) + log_north, log(0.5) + log_south)
}

sunspots_time_varying_loglik <- function(x, u, theta) {
  sum(sunspots_time_varying_log_density(x = x, u = u, theta = theta))
}

sunspots_time_varying_initial_theta <- function(x, u, nu_eps = 1e-6) {
  z <- pmin(pmax(as.numeric(x[, 3L]), -1), 1)
  early <- u <= stats::quantile(u, probs = 1 / 3, names = FALSE)
  late <- u >= stats::quantile(u, probs = 2 / 3, names = FALSE)
  path_start <- function(values, index_early, index_late) {
    fallback <- stats::median(values, na.rm = TRUE)
    start <- stats::median(values[index_early], na.rm = TRUE)
    end <- stats::median(values[index_late], na.rm = TRUE)
    if (!is.finite(start)) start <- fallback
    if (!is.finite(end)) end <- fallback * 0.8
    start <- min(max(start, 0.08), 0.85)
    end <- min(max(end, nu_eps * 2), start - 0.01)
    if (end <= nu_eps) end <- max(nu_eps * 2, start * 0.8)
    c(a = start, b = end - start)
  }

  north_values <- z[z >= 0]
  south_values <- -z[z < 0]
  north <- path_start(z[z >= 0], early[z >= 0], late[z >= 0])
  south <- path_start(-z[z < 0], early[z < 0], late[z < 0])
  fitted_mean <- ifelse(z >= 0, north[["a"]] + north[["b"]] * u,
                        -(south[["a"]] + south[["b"]] * u))
  variance <- stats::var(z - fitted_mean)
  c_value <- min(max(1 / (2 * max(variance, 0.01)), 1), 200)

  list(a_N = north[["a"]], b_N = north[["b"]],
       a_S = south[["a"]], b_S = south[["b"]], c = c_value)
}

sunspots_time_varying_shared_theta <- function(theta,
                                                nu_eps = 1e-6,
                                                c_min = 1e-8,
                                                c_max = 1e6) {
  theta <- sunspots_time_varying_validate_theta(
    theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max
  )
  list(
    a_N = mean(c(theta$a_N, theta$a_S)),
    b_N = mean(c(theta$b_N, theta$b_S)),
    a_S = mean(c(theta$a_N, theta$a_S)),
    b_S = mean(c(theta$b_N, theta$b_S)),
    c = theta$c
  )
}

fit_sunspots_time_varying_asymmetric_mixture <- function(x,
                                                          u,
                                                          hemisphere_regression = "asymmetric",
                                                          control = list()) {
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  u <- as.numeric(u)
  if (length(u) != nrow(x)) stop("`u` must have one entry per row of `x`.")
  if (any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) stop("`u` must lie in (0, 1).")

  nu_eps <- as.numeric(control$nu_eps %||% 1e-6)
  c_min <- as.numeric(control$c_min %||% 1e-8)
  c_max <- as.numeric(control$c_max %||% 1e6)
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  optim_control <- control$optim_control %||% list(maxit = 500L, reltol = 1e-10)
  base_start <- control$start_theta %||%
    sunspots_time_varying_initial_theta(x = x, u = u, nu_eps = nu_eps)
  if (identical(hemisphere_regression, "shared")) {
    base_start <- sunspots_time_varying_shared_theta(
      base_start, nu_eps = nu_eps, c_min = c_min, c_max = c_max
    )
  }
  base_start <- sunspots_time_varying_validate_theta(
    base_start, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )

  start_thetas <- list(
    base_start,
    modifyList(base_start, list(c = 10)),
    modifyList(base_start, list(c = 30)),
    modifyList(base_start, list(c = 60))
  )
  start_thetas <- lapply(start_thetas, sunspots_time_varying_validate_theta,
                         nu_eps = nu_eps, c_min = c_min, c_max = c_max,
                         hemisphere_regression = hemisphere_regression)

  objective <- function(par) {
    theta <- sunspots_time_varying_unpack_par(
      par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
      hemisphere_regression = hemisphere_regression
    )
    value <- -sunspots_time_varying_loglik(x = x, u = u, theta = theta)
    if (is.finite(value)) value else .Machine$double.xmax / 100
  }

  fits <- lapply(start_thetas, function(start_theta) {
    par0 <- sunspots_time_varying_pack_theta(
      start_theta, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
      hemisphere_regression = hemisphere_regression
    )
    try(stats::optim(par = par0, fn = objective, method = "BFGS", control = optim_control), silent = TRUE)
  })
  fits <- Filter(function(fit) !inherits(fit, "try-error") && is.finite(fit$value), fits)
  if (length(fits) == 0L) stop("All time-varying mixture optimizations failed.")

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "value"))]]
  theta_hat <- sunspots_time_varying_unpack_par(
    best$par, nu_eps = nu_eps, c_min = c_min, c_max = c_max,
    hemisphere_regression = hemisphere_regression
  )
  c(theta_hat, list(
    loglik = -best$value,
    opt = best,
    n_starts = length(start_thetas),
    n_successful_starts = length(fits),
    start_theta = base_start
  ))
}

sunspots_time_varying_axis_density <- function(z, u, theta) {
  z <- as.numeric(z)
  u <- as.numeric(u)
  if (length(z) != length(u)) stop("`z` and `u` must have the same length.")
  nu <- sunspots_time_varying_nu(u, theta)
  exp(rotational_logsumexp2(
    log(0.5) + sunspots_time_varying_component_log_axis_density(z, theta$c, nu$north),
    log(0.5) + sunspots_time_varying_component_log_axis_density(-z, theta$c, nu$south)
  ))
}

sunspots_time_varying_conditional_pit <- function(z, u, theta) {
  z <- as.numeric(z)
  nu <- sunspots_time_varying_nu(u, theta)
  mapply(function(z_i, nu_N_i, nu_S_i) {
    0.5 * small_circle_axis_cdf(z_i, kappa = theta$c, nu = nu_N_i) +
      0.5 * (1 - small_circle_axis_cdf(-z_i, kappa = theta$c, nu = nu_S_i))
  }, z, nu$north, nu$south)
}

prepare_sunspots_cycle23_time_varying_data <- function(input_csv,
                                                        start_date = "1997-06-01",
                                                        end_date = "2006-01-01") {
  if (!file.exists(input_csv)) source(prep_path_sunspots_time_varying)
  if (!file.exists(input_csv)) stop(sprintf("Input CSV not found: %s", input_csv))

  df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  required <- c("cycle", "date", "x1", "x2", "x3")
  if (!all(required %in% names(df))) {
    stop(sprintf("Input CSV is missing: %s", paste(setdiff(required, names(df)), collapse = ", ")))
  }
  if (!all(df$cycle == 23L)) stop("The time-varying runner is restricted to cycle 23 data.")

  dates <- as.POSIXct(df$date, tz = "UTC")
  start_time <- as.POSIXct(as.Date(start_date), tz = "UTC")
  end_time <- as.POSIXct(as.Date(end_date), tz = "UTC")
  if (is.na(start_time) || is.na(end_time) || start_time >= end_time) {
    stop("`start_date` and `end_date` must define a non-empty chronological interval.")
  }
  keep <- !is.na(dates) & dates >= start_time & dates < end_time
  retained <- df[keep, , drop = FALSE]
  retained$date <- dates[keep]
  retained <- retained[order(retained$date, retained$NOAA), , drop = FALSE]
  if (nrow(retained) == 0L) stop("The requested date interval retains no observations.")

  retained$time_rank_mid <- rank(as.numeric(retained$date), ties.method = "average")
  retained$u <- (retained$time_rank_mid - 0.5) / nrow(retained)
  retained
}

plot_sunspots_time_varying_diagnostics <- function(data, theta, output_dir, n_blocks = 4L) {
  u_grid <- seq(0, 1, length.out = 401L)
  u_grid <- pmin(pmax(u_grid, 1e-8), 1 - 1e-8)
  nu_grid <- sunspots_time_varying_nu(u_grid, theta)
  latitude <- asin(pmin(pmax(data$x3, -1), 1)) * 180 / pi

  migration_path <- file.path(output_dir, "cycle23_time_varying_latitude_paths.png")
  grDevices::png(migration_path, width = 1400, height = 900, res = 140)
  plot(data$u, latitude, pch = 16, cex = 0.35,
       col = grDevices::adjustcolor(ifelse(latitude >= 0, "#c43c39", "#2b6cb0"), alpha.f = 0.22),
       xlab = "Uniformized time U", ylab = "Latitude (degrees)",
       main = "Cycle 23: fitted north and south small-circle paths")
  lines(u_grid, asin(nu_grid$north) * 180 / pi, col = "#8b0000", lwd = 3)
  lines(u_grid, -asin(nu_grid$south) * 180 / pi, col = "#003f7f", lwd = 3)
  legend("topright", legend = c("North observations", "South observations", "Fitted north band", "Fitted south band"),
         col = c("#c43c39", "#2b6cb0", "#8b0000", "#003f7f"),
         pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3), bty = "n")
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  breaks <- quantile(data$u, probs = seq(0, 1, length.out = n_blocks + 1L), names = FALSE)
  data$time_block <- cut(data$u, breaks = unique(breaks), include.lowest = TRUE)
  density_path <- file.path(output_dir, "cycle23_time_varying_axial_densities.png")
  grDevices::png(density_path, width = 1400, height = 1000, res = 140)
  old_par <- par(mfrow = grDevices::n2mfrow(n_blocks), mar = c(4, 4, 3, 1))
  z_grid <- seq(-1, 1, length.out = 501L)
  for (block in levels(data$time_block)) {
    block_data <- data[data$time_block == block, , drop = FALSE]
    empirical <- stats::density(block_data$x3, from = -1, to = 1, n = length(z_grid))
    fitted <- rowMeans(vapply(block_data$u, function(u_value) {
      sunspots_time_varying_axis_density(z_grid, rep(u_value, length(z_grid)), theta)
    }, numeric(length(z_grid))))
    plot(empirical, xlim = c(-1, 1), ylim = c(0, max(empirical$y, fitted) * 1.05),
         lwd = 2, col = "black", xlab = "z = x3", ylab = "Axial density",
         main = sprintf("U in %s (n = %d)", block, nrow(block_data)))
    lines(z_grid, fitted, col = "#1f78b4", lwd = 2)
    legend("topright", legend = c("Kernel estimate", "Fitted conditional average"),
           col = c("black", "#1f78b4"), lwd = 2, bty = "n", cex = 0.85)
    grid(col = "#e5e5e5")
  }
  par(old_par)
  grDevices::dev.off()

  pit <- sunspots_time_varying_conditional_pit(data$x3, data$u, theta)
  pit_path <- file.path(output_dir, "cycle23_time_varying_conditional_pit.png")
  grDevices::png(pit_path, width = 1400, height = 700, res = 140)
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(pit, breaks = 20, freq = FALSE, col = "#9ecae1", border = "white",
       xlab = "Conditional PIT", main = "Conditional PIT histogram", xlim = c(0, 1))
  abline(h = 1, col = "#8b0000", lwd = 2)
  qqplot(stats::ppoints(length(pit)), sort(pit), xlim = c(0, 1), ylim = c(0, 1),
         pch = 16, cex = 0.45, col = "#2b6cb0", xlab = "Uniform quantiles",
         ylab = "Conditional PIT quantiles", main = "Conditional PIT Q-Q plot")
  abline(0, 1, col = "#8b0000", lwd = 2)
  par(old_par)
  grDevices::dev.off()

  list(migration_path = migration_path, density_path = density_path, pit_path = pit_path, pit = pit)
}

run_sunspots_cycle23_time_varying_asymmetric_mixture <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle23_time_varying_asymmetric_mixture"),
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    n_blocks = 4L,
    hemisphere_regression = "asymmetric",
    control = list()) {
  hemisphere_regression <- sunspots_time_varying_normalize_hemisphere_regression(hemisphere_regression)
  data <- prepare_sunspots_cycle23_time_varying_data(
    input_csv = input_csv, start_date = start_date, end_date = end_date
  )
  x <- as.matrix(data[, c("x1", "x2", "x3")])
  theta_hat <- fit_sunspots_time_varying_asymmetric_mixture(
    x = x, u = data$u, hemisphere_regression = hemisphere_regression, control = control
  )
  log_density <- sunspots_time_varying_log_density(x, data$u, theta_hat)
  loglik <- sum(log_density)
  n <- nrow(data)
  aic <- 2 * theta_hat$n_parameters - 2 * loglik
  bic <- theta_hat$n_parameters * log(n) - 2 * loglik

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- plot_sunspots_time_varying_diagnostics(data, theta_hat, output_dir, n_blocks = as.integer(n_blocks))
  data$conditional_log_density <- log_density
  data$conditional_pit <- diagnostics$pit
  utils::write.csv(data, file.path(output_dir, "cycle23_time_varying_retained_data.csv"), row.names = FALSE)
  utils::write.csv(data.frame(
    a_N = theta_hat$a_N, b_N = theta_hat$b_N, nu_N_start = theta_hat$a_N, nu_N_end = theta_hat$nu_N_end,
    a_S = theta_hat$a_S, b_S = theta_hat$b_S, nu_S_start = theta_hat$a_S, nu_S_end = theta_hat$nu_S_end,
    c = theta_hat$c, hemisphere_regression = hemisphere_regression,
    n_parameters = theta_hat$n_parameters, loglik = loglik, aic = aic, bic = bic,
    convergence = theta_hat$opt$convergence, optim_message = theta_hat$opt$message %||% "",
    n = n, start_date = start_date, end_date_exclusive = end_date,
    stringsAsFactors = FALSE
  ), file.path(output_dir, "cycle23_time_varying_fit.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  invisible(list(data = data, theta_hat = theta_hat, loglik = loglik, aic = aic, bic = bic,
                 diagnostics = diagnostics, output_dir = output_dir))
}

parse_sunspots_time_varying_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat("Options: --input_csv=PATH --output_dir=PATH --start_date=YYYY-MM-DD --end_date=YYYY-MM-DD --n_blocks=INTEGER --hemisphere_regression=asymmetric|shared\n")
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key %in% c("input_csv", "output_dir", "start_date", "end_date", "hemisphere_regression")) out[[key]] <- value
    if (key == "n_blocks") out[[key]] <- as.integer(value)
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_time_varying_asymmetric_mixture, parse_sunspots_time_varying_args())
}
