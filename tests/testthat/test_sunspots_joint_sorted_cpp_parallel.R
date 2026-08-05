library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path(
  "real_data",
  "sunspots",
  "sunspots_cycle23_joint_time_space.R"
))
source(file.path("bootstrap", "model_specs.R"))
source(file.path("bootstrap", "multiplier_bootstrap.R"))
source(file.path(
  "bootstrap",
  "sunspots_joint_time_space_model_spec.R"
))

sorted_cpp_fixture <- function(n = 16L, seed = 20260805L) {
  control <- list(
    hemisphere_regression = "asymmetric",
    time_quad_n = 8L,
    profile_l_max = 20L,
    profile_quad_n = 80L,
    center_block_size = 4L,
    ks_block_size = 4L,
    cvm_block_size = 4L,
    distance_profile_backend = "cpp"
  )
  eta <- sunspots_joint_time_canonicalize_eta(
    list(
      weight1 = 0.42,
      alpha1 = 3.5,
      beta1 = 8.5,
      alpha2 = 8,
      beta2 = 3.2
    ),
    control = control
  )
  fit <- list(
    eta_hat = eta,
    theta_hat = list(
      a_N = 0.59,
      b_N = -0.21,
      a_S = 0.51,
      b_S = -0.15,
      c = 16
    )
  )
  par <- sunspots_joint_pack_par(fit, control = control)
  simulated <- sunspots_joint_with_seed(
    seed,
    sample_sunspots_joint_time_space(
      n,
      par,
      control = control
    )
  )
  z <- cbind(simulated$x, simulated$s)
  spec <- make_sunspots_joint_time_space_spec("asymmetric")
  distances <- spec$distance_matrix(z, z, control)
  sorted_distances <- sort_distance_matrix_rows(
    distances
  )$sorted_distance_matrix

  list(
    control = control,
    fit = fit,
    data = simulated,
    z = z,
    spec = spec,
    sorted_distances = sorted_distances
  )
}

test_that("sorted joint C++ kernel matches the generic C++ and R kernels", {
  skip_if_not_installed("Rcpp")
  fixture <- sorted_cpp_fixture()
  prepared <- spec_sample_profile_sorted_prepare(
    spec = fixture$spec,
    data = fixture$z,
    sorted_distance_matrix = fixture$sorted_distances,
    theta = fixture$fit,
    control = fixture$control
  )
  rows <- 1:4
  radii <- fixture$sorted_distances[rows, , drop = FALSE]

  optimized <- sunspots_joint_profile_block_sorted(
    radii = radii,
    rho = prepared$rho[rows],
    center_s = prepared$center_s[rows],
    time_nodes = prepared$time_nodes,
    time_weights = prepared$time_weights,
    coefficients = prepared$coefficients,
    backend = "cpp"
  )
  generic_cpp <- sunspots_joint_profile_block(
    radii = radii,
    rho = prepared$rho[rows],
    center_s = prepared$center_s[rows],
    time_nodes = prepared$time_nodes,
    time_weights = prepared$time_weights,
    coefficients = prepared$coefficients,
    backend = "cpp"
  )
  reference_r <- sunspots_joint_profile_block(
    radii = radii,
    rho = prepared$rho[rows],
    center_s = prepared$center_s[rows],
    time_nodes = prepared$time_nodes,
    time_weights = prepared$time_weights,
    coefficients = prepared$coefficients,
    backend = "r"
  )

  expect_equal(optimized, generic_cpp, tolerance = 1e-12)
  expect_equal(optimized, reference_r, tolerance = 1e-12)

  unsorted <- radii
  unsorted[1L, 1:2] <- rev(unsorted[1L, 1:2])
  if (unsorted[1L, 1L] == unsorted[1L, 2L]) {
    unsorted[1L, 1L] <- unsorted[1L, 2L] + 0.1
  }
  expect_error(
    sunspots_joint_profile_block_sorted(
      radii = unsorted,
      rho = prepared$rho[rows],
      center_s = prepared$center_s[rows],
      time_nodes = prepared$time_nodes,
      time_weights = prepared$time_weights,
      coefficients = prepared$coefficients,
      backend = "cpp"
    ),
    regexp = "sorted"
  )
})

test_that("parallel lightweight observed statistics equal serial statistics", {
  skip_if_not_installed("Rcpp")
  skip_on_os("windows")

  fixture <- sorted_cpp_fixture(n = 20L, seed = 20260806L)
  serial_control <- utils::modifyList(
    fixture$control,
    list(observed_profile_n_cores = 1L)
  )
  parallel_control <- utils::modifyList(
    fixture$control,
    list(observed_profile_n_cores = 2L)
  )

  serial <- compute_sample_ks_cvm_observed_stats_light(
    spec = fixture$spec,
    normalized_data = spec_normalize_data(
      fixture$spec,
      fixture$z,
      serial_control
    ),
    sorted_distance_matrix = fixture$sorted_distances,
    theta = fixture$fit,
    control = serial_control
  )
  parallel <- compute_sample_ks_cvm_observed_stats_light(
    spec = fixture$spec,
    normalized_data = spec_normalize_data(
      fixture$spec,
      fixture$z,
      parallel_control
    ),
    sorted_distance_matrix = fixture$sorted_distances,
    theta = fixture$fit,
    control = parallel_control
  )

  expect_equal(parallel$ks, serial$ks, tolerance = 1e-12)
  expect_equal(parallel$cvm, serial$cvm, tolerance = 1e-12)
})


test_that("joint sunspots specs are certified for the generic C++ backend", {
  expect_true(distance_profile_cpp_supports_spec(
    make_sunspots_joint_time_space_spec("asymmetric")$name
  ))
  expect_true(distance_profile_cpp_supports_spec(
    make_sunspots_joint_time_space_spec("shared")$name
  ))
})
