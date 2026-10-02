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
