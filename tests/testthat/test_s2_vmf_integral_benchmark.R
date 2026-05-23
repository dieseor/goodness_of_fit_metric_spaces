utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
source(utils_path)

project_root <- if (file.exists("utils.R")) {
  getwd()
} else {
  normalizePath(file.path("..", ".."))
}
old_wd <- setwd(project_root)
on.exit(setwd(old_wd), add = TRUE)

source(file.path("convergence_empirical_process", "gaussian_process_s2_integral_benchmark.R"))

test_that("sphere distance thresholds are converted correctly", {
  expect_equal(sphere_distance_to_dot_threshold(0, "chordal"), 1)
  expect_equal(sphere_distance_to_dot_threshold(2, "chordal"), -1)
  expect_equal(sphere_distance_to_dot_threshold(0, "geodesic"), 1)
  expect_equal(sphere_distance_to_dot_threshold(pi, "geodesic"), -1)
})

test_that("S2 exact joint probability is symmetric", {
  mu <- c(0, 0, 1)
  omega1 <- s2_point_from_angles(0.2, 0.9)
  omega2 <- s2_point_from_angles(1.1, 0.8)
  prob12 <- joint_probability_vmf_s2_simple_integral(
    omega1 = omega1,
    t1 = 0.85,
    omega2 = omega2,
    t2 = 0.95,
    mu = mu,
    kappa = 2,
    distance_type = "geodesic",
    rel.tol_outer = 1e-6,
    abs.tol_outer = 1e-8,
    rel.tol_inner = 1e-7,
    abs.tol_inner = 1e-9
  )
  prob21 <- joint_probability_vmf_s2_simple_integral(
    omega1 = omega2,
    t1 = 0.95,
    omega2 = omega1,
    t2 = 0.85,
    mu = mu,
    kappa = 2,
    distance_type = "geodesic",
    rel.tol_outer = 1e-6,
    abs.tol_outer = 1e-8,
    rel.tol_inner = 1e-7,
    abs.tol_inner = 1e-9
  )
  expect_equal(prob12, prob21, tolerance = 1e-7)
})

test_that("S2 same-center intersection reduces to the smaller cap", {
  mu <- c(0, 0, 1)
  omega <- s2_point_from_angles(0.4, 0.7)
  t1 <- 0.70
  t2 <- 1.05

  joint_prob <- joint_probability_vmf_s2_simple_integral(
    omega1 = omega,
    t1 = t1,
    omega2 = omega,
    t2 = t2,
    mu = mu,
    kappa = 1.5,
    distance_type = "geodesic",
    rel.tol_outer = 1e-6,
    abs.tol_outer = 1e-8,
    rel.tol_inner = 1e-7,
    abs.tol_inner = 1e-9
  )

  marginal_prob <- distance_profile_vmf_s2_integral(
    omega = omega,
    mu = mu,
    kappa = 1.5,
    t_values = min(t1, t2),
    distance_type = "geodesic",
    rel.tol = 1e-7,
    abs.tol = 1e-9
  )

  expect_equal(joint_prob, marginal_prob, tolerance = 1e-7)
})

test_that("S2 exact covariance matrix smoke test returns finite symmetric matrix", {
  mu <- c(0, 0, 1)
  omega_grid <- generate_canonical_lattice(3, dim = 3)
  t_grid <- c(0.6, 1.2)

  sigma_exact <- cov_vmf_s2_simple_integral(
    omega_grid = omega_grid,
    t_grid = t_grid,
    mu = mu,
    kappa = 0.5,
    distance_type = "geodesic",
    rel.tol_outer = 1e-5,
    abs.tol_outer = 1e-7,
    rel.tol_inner = 1e-6,
    abs.tol_inner = 1e-8,
    subdivisions_outer = 120L,
    subdivisions_inner = 120L
  )

  expect_equal(dim(sigma_exact), c(6, 6))
  expect_true(all(is.finite(sigma_exact)))
  expect_equal(sigma_exact, t(sigma_exact), tolerance = 1e-8)
})

test_that("S2 closed-form special cases agree with the integral route", {
  mu <- c(0, 0, 1)
  cases <- default_s2_closed_form_cases(mu = mu, distance_type = "geodesic")
  target_case <- cases[[4]]
  truth_prob <- closed_form_joint_probability_vmf_s2_special(
    case = target_case,
    kappa = 2,
    distance_type = "geodesic"
  )
  integral_prob <- joint_probability_vmf_s2_simple_integral(
    omega1 = target_case$omega1,
    t1 = target_case$t1,
    omega2 = target_case$omega2,
    t2 = target_case$t2,
    mu = mu,
    kappa = 2,
    distance_type = "geodesic",
    rel.tol_outer = 1e-7,
    abs.tol_outer = 1e-9,
    rel.tol_inner = 1e-8,
    abs.tol_inner = 1e-10
  )
  expect_equal(integral_prob, truth_prob, tolerance = 1e-7)
})

test_that("cov_vmf accepts the integral S2 branch for the simple null", {
  source("convergence_empirical_process/gaussian_process_vmf.R")
  mu <- c(0, 0, 1)
  omega_grid <- generate_canonical_lattice(3, dim = 3)
  t_grid <- c(0.6, 1.2)

  sigma_exact <- cov_vmf(
    omega_grid = omega_grid,
    t_grid = t_grid,
    mu = mu,
    kappa = 0.5,
    distance_type = "geodesic",
    n_cores = 1,
    h0 = "simple",
    cov_method = "integral_s2_simple"
  )

  expect_equal(dim(sigma_exact), c(6, 6))
  expect_true(all(is.finite(sigma_exact)))
  expect_equal(sigma_exact, t(sigma_exact), tolerance = 1e-8)
})
