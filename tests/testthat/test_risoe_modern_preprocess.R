library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "wind", "preprocess_risoe_modern_hvmf.R"))

test_that("select_noon_nov_dec keeps one November-December record per day and breaks noon ties early", {
  df <- data.frame(
    datetime = as.POSIXct(
      c(
        "2000-10-31 11:55:00",
        "2000-11-01 11:55:00",
        "2000-11-01 12:05:00",
        "2000-11-02 12:15:00",
        "2000-12-01 11:45:00",
        "2000-12-01 12:25:00"
      ),
      tz = "UTC"
    ),
    ws77 = c(1, 2, 3, 4, 5, 6),
    wd77 = c(10, 20, 30, 40, 50, 60),
    ws125 = c(1, 2, 3, 4, 5, 6),
    wd125 = c(10, 20, 30, 40, 50, 60)
  )

  selected <- select_noon_nov_dec(df, tie_break = "earliest", fixed_tz = "UTC")

  expect_equal(nrow(selected), 3L)
  expect_equal(format(selected$datetime[[1L]], "%Y-%m-%d %H:%M:%S", tz = "UTC"), "2000-11-01 11:55:00")
  expect_equal(format(selected$datetime[[2L]], "%Y-%m-%d %H:%M:%S", tz = "UTC"), "2000-11-02 12:15:00")
  expect_equal(format(selected$datetime[[3L]], "%Y-%m-%d %H:%M:%S", tz = "UTC"), "2000-12-01 11:45:00")
  expect_true(all(extract_calendar_fields(selected$datetime, fixed_tz = "UTC")$month %in% c(11L, 12L)))
  expect_equal(anyDuplicated(as.Date(selected$datetime, tz = "UTC")), 0L)
})

test_that("build_hvmf_wind_set normalizes speeds, wraps directions, and stays on H^2", {
  df <- data.frame(
    datetime = as.POSIXct(
      c(
        "2000-11-01 11:55:00",
        "2000-11-02 12:15:00",
        "2000-12-01 11:45:00",
        "2000-12-02 11:35:00"
      ),
      tz = "UTC"
    ),
    ws77 = c(2, 4, 0, 5),
    wd77 = c(30, 390, 90, Inf)
  )

  result <- build_hvmf_wind_set(df, speed_col = "ws77", direction_col = "wd77", height_m = 77, fixed_tz = "UTC")

  expect_equal(nrow(result), 2L)
  expect_equal(attr(result, "dropped_days"), 2L)
  expect_equal(unique(result$height_m), 77)
  expect_equal(unique(result$speed_mean_height), 3, tolerance = 1e-12)
  expect_equal(result$speed_scaled, c(2 / 3, 4 / 3), tolerance = 1e-12)
  expect_equal(result$direction_deg, c(30, 30), tolerance = 1e-12)

  expected_x0 <- cosh(result$speed_scaled)
  expected_x1 <- sinh(result$speed_scaled) * cos(result$angle_rad)
  expected_x2 <- sinh(result$speed_scaled) * sin(result$angle_rad)

  expect_equal(result$x0, expected_x0, tolerance = 1e-12)
  expect_equal(result$x1, expected_x1, tolerance = 1e-12)
  expect_equal(result$x2, expected_x2, tolerance = 1e-12)
  expect_true(all(abs(result$minkowski_norm + 1) < 1e-12))
  expect_true(all(result$month %in% c(11L, 12L)))
})
