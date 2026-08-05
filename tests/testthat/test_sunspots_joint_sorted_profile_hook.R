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
  "real_data",
  "sunspots",
  "run_sunspots_cycle23_time_varying_asymmetric_mixture_gof.R"
))

test_that(
  "sunspots sorted-profile hook is exact and prepares coefficients once",
  {
    control <- list(
      hemisphere_regression = "asymmetric",
      time_quad_n = 8L,
      profile_l_max = 12L,
      profile_quad_n = 60L,
      center_block_size = 3L,
      ks_block_size = 3L,
      cvm_block_size = 3L,
      distance_profile_backend = "r"
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
      20260805L,
      sample_sunspots_joint_time_space(
        12L,
        par,
        control = control
      )
    )

    z <- cbind(simulated$x, simulated$s)
    spec <- make_sunspots_joint_time_space_spec("asymmetric")
    normalized <- spec_normalize_data(spec, z, control)
    distances <- spec$distance_matrix(z, z, control)
    sorted_distances <- sort_distance_matrix_rows(
      distances
    )$sorted_distance_matrix

    prepared <- spec_sample_profile_sorted_prepare(
      spec = spec,
      data = normalized,
      sorted_distance_matrix = sorted_distances,
      theta = fit,
      control = control
    )
    row_indices <- c(1L, 4L, 7L)
    block_values <- spec_sample_profile_sorted_block_eval(
      spec = spec,
      data = normalized,
      sorted_distance_matrix = sorted_distances,
      theta = fit,
      row_indices = row_indices,
      prepared = prepared,
      control = control
    )

    reference_values <- matrix(
      0,
      nrow = length(row_indices),
      ncol = nrow(sorted_distances)
    )
    for (k in seq_along(row_indices)) {
      i <- row_indices[[k]]
      omega_i <- spec_observation_at_normalized(
        spec,
        normalized,
        i,
        control
      )
      reference_values[k, ] <- spec$profile_eval(
        omega_i,
        sorted_distances[i, ],
        fit,
        control
      )
    }
    expect_equal(
      block_values,
      reference_values,
      tolerance = 1e-12
    )

    original_coefficients <-
      sunspots_joint_conditional_legendre_coefficients
    coefficient_calls <- 0L
    assign(
      "sunspots_joint_conditional_legendre_coefficients",
      function(...) {
        coefficient_calls <<- coefficient_calls + 1L
        original_coefficients(...)
      },
      envir = .GlobalEnv
    )
    on.exit(
      assign(
        "sunspots_joint_conditional_legendre_coefficients",
        original_coefficients,
        envir = .GlobalEnv
      ),
      add = TRUE
    )

    hook_stats <- compute_sample_ks_cvm_observed_stats_light(
      spec = spec,
      normalized_data = normalized,
      sorted_distance_matrix = sorted_distances,
      theta = fit,
      control = control
    )
    expect_equal(coefficient_calls, 1L)

    assign(
      "sunspots_joint_conditional_legendre_coefficients",
      original_coefficients,
      envir = .GlobalEnv
    )

    fallback_spec <- spec
    fallback_spec$sample_profile_sorted_prepare <- NULL
    fallback_spec$sample_profile_sorted_block_eval <- NULL
    fallback_stats <- compute_sample_ks_cvm_observed_stats_light(
      spec = fallback_spec,
      normalized_data = normalized,
      sorted_distance_matrix = sorted_distances,
      theta = fit,
      control = control
    )

    expect_equal(hook_stats$ks, fallback_stats$ks, tolerance = 1e-12)
    expect_equal(hook_stats$cvm, fallback_stats$cvm, tolerance = 1e-12)
  }
)
