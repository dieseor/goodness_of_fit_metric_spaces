library(testthat)

utils_path <- if (file.exists("utils.R")) {
  "utils.R"
} else {
  file.path("..", "..", "utils.R")
}
suppressWarnings(source(utils_path, local = FALSE))

model_specs_path <- if (file.exists(file.path("bootstrap", "model_specs.R"))) {
  file.path("bootstrap", "model_specs.R")
} else {
  file.path("..", "..", "bootstrap", "model_specs.R")
}
suppressWarnings(source(model_specs_path, local = FALSE))

test_that("vMF S2 Legendre coefficients recover the uniform case", {
  coeffs <- vmf_s2_legendre_coefficients(kappa = 0, l_max = 8L)

  expect_equal(coeffs[[1L]], 1, tolerance = 1e-15)
  expect_true(max(abs(coeffs[-1L])) < 1e-15)
})

test_that("vMF S2 Legendre scalar profile matches the exact integral closely", {
  mu <- c(0, 0, 1)
  omega <- c(sqrt(1 - 0.4^2), 0, 0.4)
  t_values <- seq(0, pi, length.out = 41)

  for (cfg in list(list(kappa = 0.5, l_max = 10L),
                   list(kappa = 2, l_max = 10L),
                   list(kappa = 10, l_max = 20L))) {
    exact <- theoretical_distance_profile_vmf_s2_fast(
      omega = omega,
      mu = mu,
      kappa = cfg$kappa,
      t_values = t_values,
      distance_type = "geodesic"
    )
    legendre <- distance_profile_vmf_s2_legendre(
      omega = omega,
      mu = mu,
      kappa = cfg$kappa,
      t_values = t_values,
      distance_type = "geodesic",
      l_max = cfg$l_max
    )

    expect_lt(max(abs(exact - legendre)), 5e-6)
  }
})

test_that("vMF S2 Legendre grid profile matches the tabulated route closely", {
  set.seed(123)
  mu <- c(0, 0, 1)
  omega_grid <- matrix(rnorm(30 * 3), ncol = 3)
  omega_grid <- omega_grid / sqrt(rowSums(omega_grid^2))
  t_grid <- seq(0, pi, length.out = 51)

  tabulated <- distance_profile_vmf_s2_grid(
    omega_grid = omega_grid,
    mu = mu,
    kappa = 10,
    t_grid = t_grid,
    distance_type = "geodesic",
    n_u = 4097L
  )
  legendre <- distance_profile_vmf_s2_legendre_grid(
    omega_grid = omega_grid,
    mu = mu,
    kappa = 10,
    t_grid = t_grid,
    distance_type = "geodesic",
    l_max = 20L
  )

  expect_lt(max(abs(tabulated - legendre)), 5e-6)
})

test_that("vMF model spec Legendre matrix evaluators agree with current routes", {
  if (!requireNamespace("rotasym", quietly = TRUE)) {
    skip("rotasym not installed")
  }

  set.seed(99)
  x <- rotasym::r_vMF(10, mu = c(1, 0, 0), kappa = 2)
  x <- normalize_vmf_data(x)
  theta <- list(mu = c(1, 0, 0), kappa = 2)
  spec <- make_vmf_spec(distance_type = "geodesic", unknown_param = "xi")
  omega_grid <- generate_canonical_lattice(8, dim = 3)
  t_grid <- seq(0.1, 2.8, length.out = 9)
  distance_matrix <- acos(pmax(pmin(x %*% t(x), 1), -1))

  tabulated_grid <- spec$profile_matrix_eval(
    omega_grid = omega_grid,
    t_grid = t_grid,
    theta = theta,
    control = list(vmf_profile_method = "tabulated", vmf_profile_n_u = 4097L)
  )
  legendre_grid <- spec$profile_matrix_eval(
    omega_grid = omega_grid,
    t_grid = t_grid,
    theta = theta,
    control = list(vmf_profile_method = "legendre", vmf_profile_l_max = 10L)
  )

  tabulated_cvm <- spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(vmf_profile_method = "tabulated", vmf_profile_n_u = 4097L)
  )
  legendre_cvm <- spec$sample_profile_matrix_eval(
    data = x,
    distance_matrix = distance_matrix,
    theta = theta,
    control = list(vmf_profile_method = "legendre", vmf_profile_l_max = 10L)
  )

  expect_lt(max(abs(tabulated_grid - legendre_grid)), 5e-6)
  expect_lt(max(abs(tabulated_cvm - legendre_cvm)), 5e-6)
})
