library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "wind", "run_risoe_jensen_like_all.R"))

test_that("make_risoe_all_window_configs uses only complete years and exposes max available", {
  complete_years <- c(1996:2001, 2003:2007)
  configs <- make_risoe_all_window_configs(complete_years)

  expect_length(configs, 3L)
  expect_identical(configs[[1L]]$years, c(1996:2001, 2003:2006))
  expect_identical(configs[[2L]]$years, c(1997:2001, 2003:2007))
  expect_identical(configs[[3L]]$years, complete_years)
  expect_identical(configs[[1L]]$n_years_nominal, 10L)
  expect_identical(configs[[3L]]$n_years_nominal, 11L)
})

test_that("augment_hvmf_case_df preserves Jensen-like scaling columns", {
  datetime <- as.POSIXct(c("2001-11-03 11:55:00", "2001-12-03 11:55:00"), tz = "UTC")
  df_case <- data.frame(
    datetime = datetime,
    year = c(2001L, 2001L),
    month = c(11L, 12L),
    day = c(3L, 3L),
    hour = c(11L, 11L),
    minute = c(55L, 55L),
    height_m = c(77, 77),
    speed = c(5, 10),
    direction_deg = c(180, 270),
    speed_mean_height = c(7.5, 7.5),
    speed_scaled = c(2 / 3, 4 / 3),
    angle_rad = c(pi, 3 * pi / 2),
    x0 = c(cosh(2 / 3), cosh(4 / 3)),
    x1 = c(-sinh(2 / 3), 0),
    x2 = c(0, -sinh(4 / 3)),
    minkowski_norm = c(-1, -1),
    stringsAsFactors = FALSE
  )

  augmented <- augment_hvmf_case_df(
    df_case = df_case,
    dataset_id = "demo",
    pattern_name = "set12",
    window_id = "window_x",
    source_file = "risoe_m_all.nc"
  )

  expect_identical(augmented$dataset_id, c("demo", "demo"))
  expect_identical(augmented$pattern, c("set12", "set12"))
  expect_identical(augmented$window_id, c("window_x", "window_x"))
  expect_equal(unique(augmented$speed_mean), 7.5)
  expect_true(all(abs(augmented$minkowski_norm + 1) < 1e-12))
})
