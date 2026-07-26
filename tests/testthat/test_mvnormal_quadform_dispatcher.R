# EN DUDA (2026-07-26): tests for the shared MVN/logistic-Gaussian dispatcher.
library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "model_specs.R"))

with_mock_quadform_backend <- function(mock, code) {
  backend_environment <- environment(mvnormal_quadform_run_farebrother)
  original <- get("mvnormal_quadform_backend_call", envir = backend_environment, inherits = FALSE)
  assign("mvnormal_quadform_backend_call", mock, envir = backend_environment)
  on.exit(assign("mvnormal_quadform_backend_call", original, envir = backend_environment), add = TRUE)
  force(code)
}

test_that("generic controls default to auto and take precedence over LG aliases", {
  expect_identical(mvnormal_quadform_settings(list())$method, "auto")
  expect_identical(
    mvnormal_quadform_settings(list(
      mvnormal_quadform_method = "davies",
      logistic_gaussian_quadform_method = "hbe"
    ))$method,
    "davies"
  )
})

test_that("the shared evaluator is used by both MVN and logistic-Gaussian adapters", {
  mvn_theta <- normalize_mvnormal_theta(list(mu = c(0, 0), Sigma = diag(2)))
  mvn_probability <- make_mvnormal_spec()$profile_eval(
    omega = c(0, 0), t = 1,
    theta = mvn_theta,
    control = list(mvnormal_quadform_method = "auto")
  )
  expect_equal(mvn_probability, stats::pchisq(1, df = 2), tolerance = 1e-12)

  lg_theta <- normalize_logistic_gaussian_theta(list(mu_ilr = 0, Sigma_ilr = matrix(1, 1, 1)))
  lg_probability <- make_logistic_gaussian_spec()$profile_eval(
    omega = c(0.5, 0.5), t = 1,
    theta = lg_theta,
    control = list(mvnormal_quadform_method = "auto")
  )
  expect_equal(lg_probability, stats::pchisq(1, df = 1), tolerance = 1e-12)
})

test_that("both mathematical pre-screens bypass Farebrother", {
  calls <- character()
  mock <- function(method, ...) {
    calls <<- c(calls, method)
    switch(method,
      farebrother = list(Qq = 0.5, ifault = 0L),
      imhof = list(Qq = 0.4, abserr = 1e-6),
      davies = list(Qq = 0.4, ifault = 0L)
    )
  }
  with_mock_quadform_backend(mock, {
    kappa <- mvnormal_quadform_cdf_diagnostics(1, c(1, 1e-5), c(1, 1), c(0, 0))
    expect_true(kappa$details[[1]]$pre_screen)
    expect_false("farebrother" %in% calls)

    calls <<- character()
    a0 <- mvnormal_quadform_cdf_diagnostics(1600, c(1, 0.9), c(1, 1), c(800, 800))
    expect_true(a0$details[[1]]$pre_screen)
    expect_false("farebrother" %in% calls)
  })
})

test_that("a nonzero Farebrother ifault enters the Imhof rescue chain", {
  calls <- character()
  mock <- function(method, ...) {
    calls <<- c(calls, method)
    switch(method,
      farebrother = list(Qq = 0.5, ifault = 1L),
      imhof = list(Qq = 0.4, abserr = 1e-5),
      davies = list(Qq = 0.4, ifault = 0L)
    )
  }
  with_mock_quadform_backend(mock, {
    result <- mvnormal_quadform_cdf_diagnostics(1, c(1, 0.9), c(1, 1), c(0, 0))
    expect_identical(result$details[[1]]$status, "imhof_usual")
    expect_equal(result$probability, 0.6)
    expect_identical(calls, c("farebrother", "imhof"))
  })
})

test_that("Imhof clips only an excursion covered by its reported error", {
  calls <- character()
  mock <- function(method, ...) {
    calls <<- c(calls, method)
    switch(method,
      farebrother = list(Qq = 0.5, ifault = 4L),
      imhof = list(Qq = 1 + 5e-5, abserr = 1e-4),
      davies = list(Qq = 0.5, ifault = 0L)
    )
  }
  with_mock_quadform_backend(mock, {
    result <- mvnormal_quadform_cdf_diagnostics(1, c(1, 0.9), c(1, 1), c(0, 0))
    expect_identical(result$details[[1]]$status, "imhof_usual")
    expect_identical(result$probability, 0)
  })
  expect_true(is.na(mvnormal_quadform_safe_clip(-2e-4, 1e-4)))
})

test_that("Davies is used only after both Imhof attempts fail its error criterion", {
  calls <- character()
  mock <- function(method, ...) {
    calls <<- c(calls, method)
    switch(method,
      farebrother = list(Qq = 0.5, ifault = 4L),
      imhof = list(Qq = 0.4, abserr = 2e-4),
      davies = list(Qq = 0.25, ifault = 0L)
    )
  }
  with_mock_quadform_backend(mock, {
    result <- mvnormal_quadform_cdf_diagnostics(1, c(1, 0.9), c(1, 1), c(0, 0))
    expect_identical(result$details[[1]]$status, "davies")
    expect_equal(result$probability, 0.75)
    expect_identical(calls, c("farebrother", "imhof", "imhof", "davies"))
  })
})

test_that("terminal MC is reproducible, restores RNG, and is shared across radii", {
  mock <- function(method, ...) {
    switch(method,
      farebrother = list(Qq = 0.5, ifault = 4L),
      imhof = list(Qq = 0.4, abserr = 1),
      davies = list(Qq = 0.4, ifault = 1L)
    )
  }
  control <- list(
    mvnormal_quadform_mc_abs_error = 0.1,
    mvnormal_quadform_mc_batch_size = 1000L,
    mvnormal_quadform_mc_seed = 901L,
    mvnormal_quadform_mc_label = "test shared radii"
  )
  set.seed(902)
  rng_before <- .Random.seed
  dispatcher_environment <- environment(mvnormal_quadform_cdf_diagnostics)
  original_mc <- get("mvnormal_quadform_mc_cdf", envir = dispatcher_environment, inherits = FALSE)
  mc_calls <- 0L
  assign("mvnormal_quadform_mc_cdf", function(...) {
    mc_calls <<- mc_calls + 1L
    original_mc(...)
  }, envir = dispatcher_environment)
  on.exit(assign("mvnormal_quadform_mc_cdf", original_mc, envir = dispatcher_environment), add = TRUE)
  with_mock_quadform_backend(mock, {
    first <- mvnormal_quadform_cdf_diagnostics(c(0.8, 1.4), c(1, 0.8), c(1, 1), c(0, 0), control)
    expect_identical(vapply(first$details, `[[`, character(1), "status"), c("mc", "mc"))
    expect_identical(mc_calls, 1L)
    second <- mvnormal_quadform_cdf_diagnostics(c(0.8, 1.4), c(1, 0.8), c(1, 1), c(0, 0), control)
    expect_equal(first$probability, second$probability)
  })
  expect_identical(.Random.seed, rng_before)
})
