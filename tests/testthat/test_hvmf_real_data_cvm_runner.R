library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "wind", "run_hvmf_real_data_cvm.R"))

test_that("coerce_flag_values handles logical-like inputs robustly", {
  expect_equal(coerce_flag_values(c("True", "False", "1", "0", "", "no")), c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(coerce_flag_values(c(TRUE, NA, FALSE), default = FALSE), c(TRUE, FALSE, FALSE))
  expect_error(coerce_flag_values("maybe"), "Could not coerce flag values")
})

test_that("speed-direction mapping lands on H^2 and passes verification", {
  coords <- build_h2_from_speed_direction(
    speed_scaled = c(0, 1, 1.5),
    direction_deg = c(0, 90, 180)
  )

  verified <- verify_h2_matrix(coords, tol = 1e-10)

  expect_equal(nrow(verified$data_matrix), 3L)
  expect_true(all(abs(verified$minkowski_norm + 1) < 1e-12))
  expect_true(verified$max_minkowski_error < 1e-12)
})

test_that("prepare_jensen_hvmf_dataset filters excluded rows and builds coordinates", {
  tmp_csv <- tempfile(fileext = ".csv")
  df <- data.frame(
    speed_scaled = c(0.5, 1.0, 1.5),
    direction_deg = c(0, 90, 180),
    excluded_manual = c("False", "True", "False"),
    user_marked_added = c("False", "False", "True"),
    stringsAsFactors = FALSE
  )
  write.csv(df, tmp_csv, row.names = FALSE)
  on.exit(unlink(tmp_csv), add = TRUE)

  prepared <- prepare_jensen_hvmf_dataset(tmp_csv, tol = 1e-10)

  expect_equal(prepared$n, 2L)
  expect_true(prepared$max_minkowski_error < 1e-12)
  expect_match(prepared$notes, "included 1 manually added points; excluded_manual removed 1 rows")
})
