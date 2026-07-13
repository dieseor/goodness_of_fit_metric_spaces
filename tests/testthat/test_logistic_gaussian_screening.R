library(testthat)

skip_if_not_installed("compositions")

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"))

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
  expect_equal(result$settings$quadform_method, "hbe")
  expect_equal(result$bootstrap$engine, "multiplier_bootstrap_logistic_gaussian")
  expect_true(result$inference$ks$p_value >= 0 && result$inference$ks$p_value <= 1)
  expect_true(result$inference$cvm$p_value >= 0 && result$inference$cvm$p_value <= 1)
  expect_true(result$classification$diagnosis %in% c("promising", "borderline", "reject", "problematic"))
  expect_equal(
    dim(result$observed$empirical_profile),
    c(nrow(result$grid$omega), length(result$grid$t_grid))
  )
})

test_that("the final registry contains the 28 specified compositions and preserved externals", {
  expected_parts <- list(
    Aar_oxides = c("SiO2", "TiO2", "Al2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "Fe2O3t"),
    Activity10 = c("teac", "cons", "admi", "rese", "wake", "slee"),
    Activity31 = c("teac", "cons", "admi", "rese", "wake", "slee"),
    AnimalVegetation = c("disc", "spick", "din", "spin"),
    ArcticLake = c("sand", "silt", "clay"),
    Bayesite = c("A", "B", "C", "D"),
    Boxite = c("A", "B", "C", "D", "E"),
    ClamEast = c("dl", "dm", "ds", "ll", "lm", "ls"),
    ClamWest = c("dl", "dm", "ds", "ll", "lm", "ls"),
    Coxite = c("A", "B", "C", "D", "E"),
    DiagnosticProb = c("A", "B", "C"),
    Firework = c("a", "b", "c", "d", "e"),
    Hongite = c("A", "B", "C", "D", "E"),
    HouseholdExp = c("Housing", "Food", "Other", "Services"),
    Hydrochem = c("H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3", "PO4", "SO4", "HCO3", "TOC"),
    juraset = c("Cd", "Cu", "Pb", "Co", "Cr", "Ni", "Zn"),
    Kongite = c("A", "B", "C", "D", "E"),
    Metabolites = c("met1", "met2", "met3"),
    PogoJump = c("yat", "yee", "sam"),
    Sediments = c("sand", "silt", "clay"),
    SerumProtein = c("a", "b", "c", "d"),
    ShiftOperators = c("A", "B", "C", "D"),
    SkyeAFM = c("A", "F", "M"),
    Supervisor = c("C", "D", "E", "F"),
    WhiteCells_microscopic = c("mG", "mL", "mM"),
    WhiteCells_image = c("iG", "iL", "iM"),
    Yatquat_preference = c("prFL", "prSK", "prST"),
    Yatquat_panel = c("paFL", "paSK", "paST")
  )
  expected_entries <- stats::setNames(names(expected_parts), names(expected_parts))
  expected_entries[c("Aar_oxides", "WhiteCells_microscopic", "WhiteCells_image", "Yatquat_preference", "Yatquat_panel")] <-
    c("Aar", "WhiteCells", "WhiteCells", "Yatquat", "Yatquat")

  expect_identical(names(composition_registry), names(expected_parts))
  expect_equal(length(composition_registry), 28L)
  expect_equal(anyDuplicated(names(composition_registry)), 0L)
  expect_identical(lapply(composition_registry, `[[`, "parts"), expected_parts)
  expect_identical(
    vapply(composition_registry, `[[`, character(1), "data_entry"),
    expected_entries
  )
  expect_identical(
    vapply(composition_registry, `[[`, character(1), "object_name"),
    expected_entries
  )
  expect_identical(
    default_logistic_gaussian_screening_datasets(),
    c(names(expected_parts), external_logistic_gaussian_screening_datasets())
  )

  registry <- logistic_gaussian_screening_dataset_registry()
  expect_true(all(names(expected_parts) %in% names(registry)))
  expect_false(any(c("Aar", "AarMajorOxides", "WhiteCells", "Blood23", "Glacial", "Skulls", "jura259") %in% names(registry)))
  expect_true(all(external_logistic_gaussian_screening_datasets() %in% names(registry)))
})

test_that("every canonical dataset loads its exact selected components", {
  for (dataset_name in names(composition_registry)) {
    prepared <- prepare_composition_dataset(dataset_name)
    expect_equal(prepared$status, "ok", info = dataset_name)
    expect_identical(prepared$component_names, composition_registry[[dataset_name]]$parts, info = dataset_name)
    expect_false(prepared$has_zeros, info = dataset_name)
    expect_false(prepared$has_missing, info = dataset_name)
    expect_equal(prepared$n_duplicate_rows, 0L, info = dataset_name)
    expect_equal(
      unname(rowSums(prepared$X_closed)),
      rep(1, prepared$n),
      tolerance = 1e-12,
      info = dataset_name
    )
  }
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
