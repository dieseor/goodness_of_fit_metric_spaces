library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("wind", "run_risoe_jensen_like_screening.R"))

test_that("month scope and day pattern filters reproduce Jensen-like subsets", {
  df <- data.frame(
    datetime = as.POSIXct(
      c(
        "2000-11-01 11:55:00",
        "2000-11-03 11:55:00",
        "2000-11-07 11:55:00",
        "2000-12-05 11:55:00",
        "2000-12-30 11:55:00"
      ),
      tz = "UTC"
    ),
    year = c(2000L, 2000L, 2000L, 2000L, 2000L),
    month = c(11L, 11L, 11L, 12L, 12L),
    day = c(1L, 3L, 7L, 5L, 30L),
    hour = rep(11L, 5L),
    minute = rep(55L, 5L),
    ws77 = c(1, 2, 3, 4, 5),
    wd77 = c(10, 20, 30, 40, 50)
  )

  nov_only <- filter_month_scope(df, "nov")
  expect_equal(nrow(nov_only), 3L)

  pattern_set12 <- filter_day_pattern(df, c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L))
  expect_equal(pattern_set12$day, c(3L, 7L, 30L))

  built <- build_jensen_like_hvmf_set(
    selected_df = df,
    speed_col = "ws77",
    direction_col = "wd77",
    height_m = 77,
    day_pattern = c(3L, 7L, 11L, 15L, 19L, 23L, 27L, 30L),
    month_scope = "nov",
    fixed_tz = "UTC"
  )

  expect_equal(nrow(built), 2L)
  expect_true(all(built$month == 11L))
  expect_equal(built$day, c(3L, 7L))
})
