source(file.path("..", "..", "wind", "plot_risoe_clean_temporal_heatmaps.R"))

test_that("apply_height_cleaning enforces conservative cutoffs", {
  df <- data.frame(
    datetime = as.POSIXct(c(
      "2007-07-31 12:00:00",
      "2007-08-01 12:00:00",
      "2004-11-30 12:00:00",
      "2004-12-01 12:00:00"
    ), tz = "UTC"),
    stringsAsFactors = FALSE
  )

  out_77 <- apply_height_cleaning(df, height_m = 77L)
  expect_equal(
    as.character(as.Date(out_77$datetime, tz = "UTC")),
    c("2007-07-31", "2004-11-30", "2004-12-01")
  )

  out_125 <- apply_height_cleaning(df, height_m = 125L)
  expect_equal(
    as.character(as.Date(out_125$datetime, tz = "UTC")),
    "2004-11-30"
  )
})

test_that("build_temporal_scale anchors October-November-December correctly", {
  scale_oct <- build_temporal_scale(c(10L))
  expect_equal(scale_oct$breaks, 1L)
  expect_equal(scale_oct$labels, "1 Oct")

  scale_nov_dec <- build_temporal_scale(c(11L, 12L))
  expect_equal(scale_nov_dec$breaks, c(1L, 31L))
  expect_equal(scale_nov_dec$labels, c("1 Nov", "1 Dec"))

  scale_all <- build_temporal_scale(c(10L, 11L, 12L))
  expect_equal(scale_all$breaks, c(1L, 32L, 62L))
  expect_equal(scale_all$labels, c("1 Oct", "1 Nov", "1 Dec"))
})
