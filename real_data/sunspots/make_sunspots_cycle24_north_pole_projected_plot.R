#!/usr/bin/env Rscript

resolve_sunspots_plot_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

runner_path_sunspots_plot <- resolve_sunspots_plot_path(
  "real_data", "sunspots", "run_sunspots_cycle24_small_circle_symmetric_mixture_gof.R"
)
source(runner_path_sunspots_plot)

make_sunspots_cycle24_north_pole_projected_plot <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    theta_csv = file.path("real_data", "sunspots", "output", "cycle24_small_circle_symmetric_mixture", "sunspots_cycle24_theta_hat.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle24_small_circle_symmetric_mixture"),
    l_max = 200L,
    quad_n = 400L,
    tol = 1e-10) {
  sunspots_df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  theta_df <- utils::read.csv(theta_csv, stringsAsFactors = FALSE)

  x <- as.matrix(sunspots_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)
  theta_hat <- list(
    mu = c(theta_df$mu_1[[1L]], theta_df$mu_2[[1L]], theta_df$mu_3[[1L]]),
    kappa = theta_df$kappa_hat[[1L]],
    nu = theta_df$nu_hat[[1L]]
  )

  north_pole <- c(0, 0, 1)
  z <- pmin(pmax(as.numeric(x %*% north_pole), -1), 1)
  projected_fit <- sunspots_projected_gof(
    z = z,
    cdf_fun = function(z_grid) {
      1 - distance_profile_small_circle_symmetric_mixture2(
        omega = north_pole,
        t_values = acos(pmin(pmax(z_grid, -1), 1)),
        mu = theta_hat$mu,
        kappa = theta_hat$kappa,
        nu = theta_hat$nu,
        distance_type = "geodesic",
        method = "legendre",
        l_max = as.integer(l_max),
        quad_n = as.integer(quad_n),
        tol = as.numeric(tol)
      )
    }
  )

  cdf_grid_df <- data.frame(
    z_grid = projected_fit$z_grid,
    cdf = projected_fit$fitted_cdf,
    empirical_cdf = projected_fit$empirical_cdf,
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    cdf_grid_df,
    file.path(output_dir, "sunspots_cycle24_north_pole_projected_cdf_grid.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(z = z),
    file.path(output_dir, "sunspots_cycle24_north_pole_projected_data.csv"),
    row.names = FALSE
  )

  plot_sunspots_projected_diagnostics(
    z = z,
    cdf_grid_df = cdf_grid_df,
    output_dir = output_dir
  )
}

if (sys.nframe() == 0L) {
  make_sunspots_cycle24_north_pole_projected_plot()
}
