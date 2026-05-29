library(testthat)

skip_if_not_installed("compositions")

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("logistic_gaussian_screening", "utils_logistic_gaussian_screening.R"))

test_that("prepare_composition_dataset returns closed positive data", {
  prepared <- prepare_composition_dataset("SkyeAFM")

  expect_equal(prepared$n, 23L)
  expect_equal(prepared$D, 3L)
  expect_false(prepared$has_zeros)
  expect_false(prepared$has_missing)
  expect_equal(prepared$status, "ok")
  expect_equal(prepared$source_package, "compositions")
  expect_equal(prepared$source_dataset_name, "SkyeAFM")
  expect_equal(unname(rowSums(prepared$X_closed)), rep(1, prepared$n), tolerance = 1e-12)
})

test_that("logistic Gaussian composite screening runs on SkyeAFM in smoke size", {
  result <- run_logistic_gaussian_screening(
    dataset_name = "SkyeAFM",
    B = 3,
    max_centers = 23,
    n_t = 8,
    bootstrap_mode = "composite_multiplier",
    seed = 2026,
    make_plots = FALSE,
    save_outputs = FALSE,
    verbose = FALSE
  )

  expect_equal(result$screening_type, "logistic_gaussian_composite_or_screening")
  expect_equal(result$bootstrap$mode, "composite_multiplier")
  expect_equal(result$settings$bootstrap_mode, "composite_multiplier")
  expect_equal(result$bootstrap$engine, "multiplier_bootstrap_logistic_gaussian")
  expect_true(result$inference$ks$p_value >= 0 && result$inference$ks$p_value <= 1)
  expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
  expect_true(result$classification$diagnosis %in% c("promising", "borderline", "reject", "problematic"))
  expect_equal(dim(result$observed$empirical_profile), c(23, length(result$grid$t_grid)))
})

test_that("default dataset registry includes expanded composite candidates", {
  datasets <- default_logistic_gaussian_screening_datasets(
    include_external = TRUE
  )

  expect_true(all(c(
    "SkyeAFM",
    "SkyeLavasComplete",
    "AarMajorOxides",
    "Sediments",
    "HouseholdExp",
    "ClamEast",
    "ClamWest",
    "ClamCombined",
    "FerrettiGut",
    "FerrettiOral",
    "Shi2015",
    "HongKongBudgetsA",
    "HongKongBudgetsB",
    "HongKongBudgetsCombined"
  ) %in% datasets))
})

test_that("not found datasets return structured not_found results", {
  result <- run_logistic_gaussian_screening(
    dataset_name = "SkyeLavasComplete",
    B = 2,
    max_centers = 10,
    n_t = 5,
    bootstrap_mode = "composite_multiplier",
    seed = 2026,
    make_plots = FALSE,
    save_outputs = FALSE,
    verbose = FALSE
  )

  if (identical(result$data_prep$status, "not_found")) {
    expect_equal(result$classification$diagnosis, "not_found")
    expect_true(is.na(result$inference$ks$p_value))
    expect_true(is.na(result$inference$cvm$p_value))
  } else {
    expect_equal(result$data_prep$status, "ok")
    expect_true(result$inference$ks$p_value >= 0 && result$inference$ks$p_value <= 1)
    expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
  }
})

test_that("ClamCombined preparation keeps site labels", {
  prepared <- prepare_composition_dataset("ClamCombined")

  expect_equal(prepared$status, "ok")
  expect_equal(prepared$n, 40L)
  expect_equal(prepared$D, 6L)
  expect_equal(prepared$group_variable, "site")
  expect_equal(as.integer(table(prepared$group)[c("ClamEast", "ClamWest")]), c(20L, 20L))
  expect_equal(unname(rowSums(prepared$X_closed)), rep(1, prepared$n), tolerance = 1e-12)
})
