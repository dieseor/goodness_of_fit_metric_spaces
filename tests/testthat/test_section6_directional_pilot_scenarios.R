library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

section6_env <- new.env(parent = globalenv())
sys.source("scripts/run_section6_new_scenarios.R", envir = section6_env)

vmf_env <- new.env(parent = globalenv())
sys.source("scripts/run_vmf_antipodal_fixed_kappa_pilot.R", envir = vmf_env)

test_that("the new projected-normal vMF pilot has concentration d and mean norm sqrt(d)", {
  design <- vmf_env$fixed_kappa_design(
    dimensions = c(2L, 5L), n_values = 50L, beta_values = 0.5,
    kappa = 1, scenario_type = "projected_normal_sqrt_d"
  )
  expect_equal(design$kappa, c(2, 5))
  expect_equal(design$projected_normal_mean_norm, sqrt(c(2, 5)))
  expect_identical(unique(design$scenario), "vmf_2_projected_normal_sqrt_d")

  set.seed(1)
  x <- vmf_env$generate_fixed_kappa_antipodal(design[1L, , drop = FALSE])
  expect_equal(ncol(x), 3L)
  expect_lt(max(abs(rowSums(x^2) - 1)), 1e-12)
})

test_that("the new HvMF pilots preserve the hyperboloid constraint", {
  catalog <- section6_env$section6_scenario_catalog()
  scenarios <- c(
    "hvmf_1_dimension_scaled_location_mixture",
    "hvmf_1_radial_c_inv_sqrt2",
    "hvmf_1_radial_c_half",
    "hvmf_2_angular_sqrt_d_concentration"
  )
  expect_true(all(vapply(scenarios, function(s) isTRUE(catalog[[s]]$experimental), logical(1))))

  for (scenario in scenarios) {
    design <- section6_env$make_section6_design(
      family = "hvmf", dimensions = c(2L, 5L), n_values = 40L,
      beta_values = 0.5, scenarios = scenario
    )
    for (i in seq_len(nrow(design))) {
      set.seed(i)
      x <- section6_env$generate_section6_sample(design[i, , drop = FALSE])
      expect_true(all(x[, 1L] > 0))
      expect_lt(max(abs(-x[, 1L]^2 + rowSums(x[, -1L, drop = FALSE]^2) + 1)), 1e-8)
    }
  }

  d <- 5
  a <- sqrt(2 / d)
  mu1 <- c(sqrt(1 + a^2), a * section6_env$section6_e(d, index = 2L))
  expect_equal(-mu1[[1L]]^2 + sum(mu1[-1L]^2), -1, tolerance = 1e-14)
})


test_that("radial HvMF pilot locations have the requested hyperbolic shifts", {
  catalog <- section6_env$section6_scenario_catalog()
  c_values <- c(
    hvmf_1_radial_c_inv_sqrt2 = 1 / sqrt(2),
    hvmf_1_radial_c_half = 1 / 2
  )

  expect_true(all(vapply(
    names(c_values),
    function(s) isTRUE(catalog[[s]]$experimental),
    logical(1)
  )))

  for (d in c(2L, 5L)) {
    mu0 <- c(sqrt(2), section6_env$section6_e(d))

    for (scenario in names(c_values)) {
      c_shift <- unname(c_values[[scenario]])
      mu1 <- section6_env$section6_hvmf_radial_mu1(d, c_shift)

      expect_equal(
        -mu1[[1L]]^2 + sum(mu1[-1L]^2),
        -1,
        tolerance = 1e-13
      )

      minkowski <- -mu0[[1L]] * mu1[[1L]] +
        sum(mu0[-1L] * mu1[-1L])

      expect_equal(
        acosh(-minkowski),
        c_shift / sqrt(d),
        tolerance = 1e-12
      )
    }
  }
})
