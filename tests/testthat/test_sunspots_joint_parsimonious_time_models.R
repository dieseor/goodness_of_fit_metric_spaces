library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path(
  "real_data",
  "sunspots",
  "sunspots_cycle23_joint_time_space.R"
))
source(file.path(
  "real_data",
  "sunspots",
  "sunspots_cycle23_joint_time_models_parsimonious.R"
))
source(file.path(
  "bootstrap",
  "sunspots_joint_time_space_parsimonious_model_spec.R"
))

parsimonious_test_theta <- function(shared = TRUE) {
  if (isTRUE(shared)) {
    return(list(
      a_N = 0.48,
      b_N = -0.2,
      a_S = 0.48,
      b_S = -0.2,
      c = 18,
      n_parameters = 3L
    ))
  }
  list(
    a_N = 0.52,
    b_N = -0.22,
    a_S = 0.45,
    b_S = -0.16,
    c = 18,
    n_parameters = 5L
  )
}

parsimonious_test_eta <- function(time_model) {
  if (identical(time_model, "beta")) {
    return(sunspots_joint_parsimonious_time_validate_eta(
      list(
        time_model = "beta",
        alpha = 3.2,
        beta = 5.4
      ),
      time_model = "beta"
    ))
  }
  sunspots_joint_parsimonious_time_validate_eta(
    list(
      time_model = "uniform_beta",
      beta_weight = 0.62,
      alpha = 3.2,
      beta = 5.4
    ),
    time_model = "uniform_beta"
  )
}

test_that("parsimonious temporal densities, CDFs, and quadratures are valid", {
  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    integral <- integrate(
      function(s) {
        sunspots_joint_parsimonious_time_density(
          s,
          eta,
          time_model = time_model
        )
      },
      lower = 0,
      upper = 1,
      subdivisions = 2000L,
      rel.tol = 1e-10
    )$value
    expect_equal(integral, 1, tolerance = 1e-9)

    expect_equal(
      sunspots_joint_parsimonious_time_cdf(
        c(-1, 0, 1, 2),
        eta,
        time_model = time_model
      ),
      c(0, 0, 1, 1),
      tolerance = 1e-14
    )

    quadrature <- sunspots_joint_parsimonious_time_quadrature(
      eta,
      n_nodes = 64L,
      time_model = time_model
    )
    expect_equal(sum(quadrature$weights), 1, tolerance = 1e-12)
    expect_lt(quadrature$mass_error, 1e-12)

    high_precision <- integrate(
      function(s) {
        exp(s) *
          sunspots_joint_parsimonious_time_density(
            s,
            eta,
            time_model = time_model
          )
      },
      lower = 0,
      upper = 1,
      subdivisions = 2000L,
      rel.tol = 1e-11
    )$value
    expect_equal(
      sum(quadrature$weights * exp(quadrature$nodes)),
      high_precision,
      tolerance = 1e-10
    )
  }
})

test_that("parsimonious temporal pack and unpack are inverse operations", {
  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    par <- sunspots_joint_parsimonious_time_pack_eta(
      eta,
      time_model = time_model
    )
    recovered <- sunspots_joint_parsimonious_time_unpack_eta(
      par,
      time_model = time_model
    )
    expect_equal(recovered$alpha, eta$alpha, tolerance = 1e-14)
    expect_equal(recovered$beta, eta$beta, tolerance = 1e-14)
    expect_equal(
      recovered$beta_weight,
      eta$beta_weight,
      tolerance = 1e-14
    )
  }
})

test_that("analytic parsimonious temporal scores match finite differences", {
  set.seed(20260806)
  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    s <- sample_sunspots_joint_parsimonious_time(
      80L,
      eta,
      time_model = time_model
    )
    par <- sunspots_joint_parsimonious_time_pack_eta(
      eta,
      time_model = time_model
    )
    analytic <- sunspots_joint_parsimonious_time_score_matrix(
      s,
      par,
      time_model = time_model
    )
    step <- 1e-6
    numeric <- vapply(seq_along(par), function(index) {
      plus <- par
      minus <- par
      plus[[index]] <- plus[[index]] + step
      minus[[index]] <- minus[[index]] - step
      plus_eta <- sunspots_joint_parsimonious_time_unpack_eta(
        plus,
        time_model = time_model
      )
      minus_eta <- sunspots_joint_parsimonious_time_unpack_eta(
        minus,
        time_model = time_model
      )
      (
        sunspots_joint_parsimonious_time_log_density(
          s,
          plus_eta,
          time_model = time_model
        ) -
          sunspots_joint_parsimonious_time_log_density(
            s,
            minus_eta,
            time_model = time_model
          )
      ) / (2 * step)
    }, numeric(length(s)))
    expect_equal(
      as.numeric(analytic),
      as.numeric(numeric),
      tolerance = 2e-5
    )
  }
})

test_that("parsimonious temporal scores have approximately zero mean", {
  set.seed(20260807)
  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    s <- sample_sunspots_joint_parsimonious_time(
      30000L,
      eta,
      time_model = time_model
    )
    par <- sunspots_joint_parsimonious_time_pack_eta(
      eta,
      time_model = time_model
    )
    score <- sunspots_joint_parsimonious_time_score_matrix(
      s,
      par,
      time_model = time_model
    )
    expect_true(all(abs(colMeans(score)) < 0.035))
  }
})

test_that("parsimonious temporal MLEs are finite on simulated samples", {
  set.seed(20260808)
  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    s <- sample_sunspots_joint_parsimonious_time(
      1000L,
      eta,
      time_model = time_model
    )
    fit <- suppressWarnings(
      fit_sunspots_joint_parsimonious_time(
        s,
        time_model = time_model,
        control = list(
          parsimonious_time_n_starts =
            if (identical(time_model, "beta")) 3L else 6L,
          parsimonious_time_nelder_mead_control =
            list(maxit = 1500L, reltol = 1e-9)
        )
      )
    )
    expect_true(is.finite(fit$loglik))
    expect_gt(fit$alpha, 0)
    expect_gt(fit$beta, 0)
    expect_gte(fit$n_successful_starts, 1L)
    expect_equal(fit$n_parameters, if (
      identical(time_model, "beta")
    ) 2L else 3L)
  }
})

test_that("uniform-beta nonidentification and mixture boundaries are detected", {
  unidentified <- list(
    time_model = "uniform_beta",
    beta_weight = 0.5,
    alpha = 1,
    beta = 1
  )
  flags <- sunspots_joint_parsimonious_time_boundary_flags(
    unidentified,
    time_model = "uniform_beta"
  )
  expect_true(flags$identification)
  expect_false(
    sunspots_joint_parsimonious_time_fast_regular(
      unidentified,
      time_model = "uniform_beta"
    )
  )

  bounds <- sunspots_joint_parsimonious_time_control()
  boundary <- list(
    time_model = "uniform_beta",
    beta_weight = bounds$weight_eps,
    alpha = 2,
    beta = 3
  )
  flags <- sunspots_joint_parsimonious_time_boundary_flags(
    boundary,
    time_model = "uniform_beta"
  )
  expect_true(flags$weight)
  expect_false(
    sunspots_joint_parsimonious_time_fast_regular(
      boundary,
      time_model = "uniform_beta"
    )
  )
})

test_that("parsimonious joint scores have the expected dimensions", {
  for (time_model in c("beta", "uniform_beta")) {
    for (shared in c(TRUE, FALSE)) {
      eta <- parsimonious_test_eta(time_model)
      fit <- list(
        time_model = time_model,
        eta_hat = eta,
        theta_hat = parsimonious_test_theta(shared)
      )
      control <- list(
        time_model = time_model,
        hemisphere_regression = if (shared) "shared" else "asymmetric"
      )
      par <- sunspots_joint_parsimonious_pack_par(fit, control)
      set.seed(if (shared) 20260809 else 20260810)
      sample <- sample_sunspots_joint_time_space_parsimonious(
        35L,
        par,
        control
      )
      score <- sunspots_joint_parsimonious_score_matrix(
        sample,
        par,
        control
      )
      expected <- sunspots_joint_parsimonious_time_n_parameters(
        time_model
      ) + if (shared) 3L else 5L
      expect_identical(dim(score), c(35L, expected))
      expect_true(all(is.finite(score)))
    }
  }
})

test_that("parsimonious temporal profiles agree between R and C++", {
  skip_if_not_installed("Rcpp")
  theta <- parsimonious_test_theta(TRUE)
  set.seed(20260811)
  radii <- matrix(runif(42L), nrow = 6L)
  rho <- runif(6L, -1, 1)
  center_s <- runif(6L)

  for (time_model in c("beta", "uniform_beta")) {
    eta <- parsimonious_test_eta(time_model)
    quadrature <- sunspots_joint_parsimonious_time_quadrature(
      eta,
      n_nodes = 18L,
      time_model = time_model
    )
    coefficients <- sunspots_joint_conditional_legendre_coefficients(
      theta,
      quadrature$nodes,
      l_max = 40L,
      quad_n = 180L
    )
    reference <- sunspots_joint_profile_block(
      radii,
      rho,
      center_s,
      quadrature$nodes,
      quadrature$weights,
      coefficients,
      backend = "r"
    )
    compiled <- sunspots_joint_profile_block(
      radii,
      rho,
      center_s,
      quadrature$nodes,
      quadrature$weights,
      coefficients,
      backend = "cpp"
    )
    expect_equal(compiled, reference, tolerance = 1e-12)
  }
})

test_that("parsimonious joint spec names are certified for C++", {
  expected <- c(
    "sunspots_joint_time_space_beta_shared",
    "sunspots_joint_time_space_beta_asymmetric",
    "sunspots_joint_time_space_uniform_beta_shared",
    "sunspots_joint_time_space_uniform_beta_asymmetric"
  )
  expect_true(all(vapply(
    expected,
    distance_profile_cpp_supports_spec,
    logical(1L)
  )))

  expect_identical(
    make_sunspots_joint_time_space_parsimonious_spec(
      "beta", "shared"
    )$name,
    expected[[1L]]
  )
  expect_identical(
    make_sunspots_joint_time_space_parsimonious_spec(
      "uniform_beta", "asymmetric"
    )$name,
    expected[[4L]]
  )
})
