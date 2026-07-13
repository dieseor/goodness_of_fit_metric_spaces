library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "sunspots", "run_sunspots_cycle23_time_varying_asymmetric_mixture.R"))

test_that("time-varying parameter map enforces global positivity and migration", {
  theta <- sunspots_time_varying_unpack_par(c(0.4, -0.7, 0.2, 0.3, 3))
  expect_lt(theta$b_N, 0)
  expect_lt(theta$b_S, 0)
  expect_gt(theta$a_N + theta$b_N, 0)
  expect_gt(theta$a_S + theta$b_S, 0)
  expect_lt(theta$a_N, 1)
  expect_lt(theta$a_S, 1)
})

test_that("shared hemisphere regression has three free parameters and the intended symmetric paths", {
  theta <- list(a_N = 0.54, b_N = -0.18, a_S = 0.54, b_S = -0.18, c = 17)
  par <- sunspots_time_varying_pack_theta(theta, hemisphere_regression = "shared")
  recovered <- sunspots_time_varying_unpack_par(par, hemisphere_regression = "shared")

  expect_length(par, 3L)
  expect_identical(recovered$n_parameters, 3L)
  expect_equal(recovered$a_N, recovered$a_S, tolerance = 0)
  expect_equal(recovered$b_N, recovered$b_S, tolerance = 0)
  expect_error(
    sunspots_time_varying_validate_theta(
      list(a_N = 0.54, b_N = -0.18, a_S = 0.50, b_S = -0.14, c = 17),
      hemisphere_regression = "shared"
    ),
    "requires a_N = a_S"
  )
})

test_that("conditional axial density and PIT use the existing small-circle law", {
  theta <- list(a_N = 0.55, b_N = -0.20, a_S = 0.50, b_S = -0.15, c = 18)
  u_value <- 0.4
  axial_integral <- integrate(
    function(z) sunspots_time_varying_axis_density(z, rep(u_value, length(z)), theta),
    lower = -1, upper = 1
  )$value
  expect_equal(axial_integral, 1, tolerance = 1e-8)

  z <- seq(-0.9, 0.9, length.out = 25)
  pit <- sunspots_time_varying_conditional_pit(z, rep(u_value, length(z)), theta)
  expect_true(all(diff(pit) >= -1e-12))
  expect_true(all(pit >= 0 & pit <= 1))
})

test_that("time ranks are recomputed after the manual date restriction", {
  temporary_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temporary_csv), add = TRUE)
  dates <- as.POSIXct(c("1997-05-31", "1997-06-02", "1997-06-02", "2005-12-31", "2006-01-01"), tz = "UTC")
  example_data <- data.frame(
    cycle = 23L, date = format(dates, tz = "UTC", usetz = TRUE), NOAA = seq_along(dates),
    x1 = 0, x2 = sqrt(1 - c(0.2, 0.3, 0.4, 0.5, 0.6)^2), x3 = c(0.2, 0.3, 0.4, 0.5, 0.6)
  )
  utils::write.csv(example_data, temporary_csv, row.names = FALSE)
  retained <- prepare_sunspots_cycle23_time_varying_data(
    temporary_csv, start_date = "1997-06-01", end_date = "2006-01-01"
  )
  expect_equal(nrow(retained), 3L)
  expect_equal(retained$time_rank_mid, c(1.5, 1.5, 3))
  expect_equal(retained$u, c(1 / 3, 1 / 3, 5 / 6))
})
