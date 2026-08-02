library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

test_that("multivariate-Normal MLE equals the defining sample moments in every audited dimension", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28100L + d)
    mu <- 0.5 * c(1, rep.int(0, d - 1L))
    Sigma <- diag(d)
    Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
    x <- mvtnorm::rmvnorm(50L, mean = mu, sigma = Sigma)

    theta <- fit_mvnormal_theta(
      x, null = list(type = "composite"), unknown_param = "both"
    )
    mu_manual <- colMeans(x)
    centered <- sweep(x, 2L, mu_manual, FUN = "-")
    Sigma_manual <- crossprod(centered) / nrow(x)

    expect_equal(theta$mu, mu_manual, tolerance = 1e-12)
    expect_equal(theta$Sigma, Sigma_manual, tolerance = 1e-12)
  }
})

test_that("weighted Normal and logistic-Gaussian MLEs use the weighted likelihood moments", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28150L + d)
    mu <- 0.5 * c(1, rep.int(0, d - 1L))
    Sigma <- diag(d)
    Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
    weights <- stats::runif(50L, min = 0.1, max = 2)
    probability_weights <- weights / sum(weights)
    x <- mvtnorm::rmvnorm(50L, mean = mu, sigma = Sigma)
    theta_normal <- fit_mvnormal_theta(
      x, weights = weights, null = list(type = "composite"), unknown_param = "both"
    )
    mu_manual <- colSums(x * probability_weights)
    centered <- sweep(x, 2L, mu_manual, FUN = "-")
    Sigma_manual <- crossprod(centered * sqrt(probability_weights))
    expect_equal(theta_normal$mu, mu_manual, tolerance = 1e-12)
    expect_equal(theta_normal$Sigma, Sigma_manual, tolerance = 1e-12)

    simplex <- logistic_gaussian_ilr_to_simplex(x, ambient_dim = d + 1L)
    theta_lg <- fit_logistic_gaussian_theta(
      simplex, weights = weights, null = list(type = "composite"), unknown_param = "both"
    )
    z <- logistic_gaussian_ilr_matrix(simplex)
    mu_manual_lg <- colSums(z * probability_weights)
    centered_lg <- sweep(z, 2L, mu_manual_lg, FUN = "-")
    Sigma_manual_lg <- crossprod(centered_lg * sqrt(probability_weights))
    expect_equal(theta_lg$mu_ilr, mu_manual_lg, tolerance = 1e-12)
    expect_equal(theta_lg$Sigma_ilr, Sigma_manual_lg, tolerance = 1e-12)
  }
})

test_that("the Gaussian fast branch uses fitted influence coordinates, not known parameters", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28250L + d)
    Sigma <- diag(d)
    Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
    x <- mvtnorm::rmvnorm(60L, mean = 0.5 * c(1, rep.int(0, d - 1L)), sigma = Sigma)
    theta <- fit_mvnormal_theta(
      x, null = list(type = "composite"), unknown_param = "both"
    )
    centered <- sweep(x, 2L, theta$mu, FUN = "-")
    Sigma_inv <- solve(theta$Sigma)
    p_sigma <- d * (d + 1L) / 2L
    score <- cbind(
      centered %*% t(Sigma_inv),
      t(vapply(seq_len(nrow(x)), function(i) {
        rr <- centered[i, , drop = FALSE]
        fast_multiplier_sym_score_to_vech(
          0.5 * (Sigma_inv %*% crossprod(rr) %*% Sigma_inv - Sigma_inv)
        )
      }, numeric(p_sigma)))
    )
    influence <- cbind(
      centered,
      t(vapply(seq_len(nrow(x)), function(i) {
        rr <- centered[i, , drop = FALSE]
        fast_multiplier_vech(crossprod(rr) - theta$Sigma)
      }, numeric(p_sigma)))
    )
    idx <- which(lower.tri(theta$Sigma, diag = TRUE), arr.ind = TRUE)
    covariance_influence <- matrix(0, d + p_sigma, d + p_sigma)
    covariance_influence[seq_len(d), seq_len(d)] <- theta$Sigma
    for (a in seq_len(p_sigma)) {
      ii <- idx[a, 1L]
      jj <- idx[a, 2L]
      for (b in seq_len(p_sigma)) {
        kk <- idx[b, 1L]
        ll <- idx[b, 2L]
        covariance_influence[d + a, d + b] <-
          theta$Sigma[ii, kk] * theta$Sigma[jj, ll] +
          theta$Sigma[ii, ll] * theta$Sigma[jj, kk]
      }
    }
    expect_equal(influence, score %*% covariance_influence, tolerance = 1e-10)
  }
})

test_that("logistic-Gaussian MLE is exactly the Normal MLE in ilr coordinates", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28200L + d)
    mu <- 0.5 * c(1, rep.int(0, d - 1L))
    Sigma <- diag(d)
    Sigma[1L, 2L] <- Sigma[2L, 1L] <- 0.75
    simplex <- rlogistic_gaussian_simplex(
      50L, mu_ilr = mu, Sigma_ilr = Sigma
    )
    z <- logistic_gaussian_ilr_matrix(simplex)
    theta <- fit_logistic_gaussian_theta(
      simplex, null = list(type = "composite"), unknown_param = "both"
    )
    mu_manual <- colMeans(z)
    centered <- sweep(z, 2L, mu_manual, FUN = "-")
    Sigma_manual <- crossprod(centered) / nrow(z)

    expect_equal(theta$mu_ilr, mu_manual, tolerance = 1e-12)
    expect_equal(theta$Sigma_ilr, Sigma_manual, tolerance = 1e-12)
    expect_equal(
      logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L), simplex,
      tolerance = 1e-12
    )
  }
})

test_that("unweighted vMF MLE solves the resultant equation in every audited dimension", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28300L + d)
    mu <- c(1, rep.int(0, d))
    x <- normalize_vmf_data(rotasym::r_vMF(100L, mu = mu, kappa = 1.5 * d))
    theta <- fit_vmf_theta(x, null = list(type = "composite"), unknown_param = "xi")
    resultant <- colMeans(x)
    r_bar <- sqrt(sum(resultant^2))
    theta_direct <- normalize_vmf_theta(list(
      mu = resultant / r_bar,
      kappa = solve_vmf_kappa_from_rbar(r_bar, q = d)
    ))
    scores <- t(vapply(seq_len(nrow(x)), function(i) {
      psi_xi(x[i, ], theta$xi, d)
    }, numeric(d + 1L)))

    # movMF is the production unweighted solver.  Its internal numerical
    # tolerance is looser than uniroot's, hence the deliberately numerical
    # (rather than machine-precision) comparison.
    expect_lt(max(abs(theta$xi - theta_direct$xi)), 1e-4)
    expect_lt(abs(A_q(theta$kappa, d) - r_bar), 1e-6)
    expect_lt(max(abs(colMeans(scores))), 1e-6)
  }
})

test_that("weighted vMF MLE solves the weighted resultant equation in every audited dimension", {
  for (d in c(2L, 5L, 10L)) {
    set.seed(28350L + d)
    mu <- c(1, rep.int(0, d))
    x <- normalize_vmf_data(rotasym::r_vMF(100L, mu = mu, kappa = 1.5 * d))
    weights <- stats::runif(nrow(x), min = 0.1, max = 2)
    theta <- fit_vmf_theta(
      x, weights = weights, null = list(type = "composite"), unknown_param = "xi"
    )
    probability_weights <- weights / sum(weights)
    resultant <- colSums(x * probability_weights)
    r_bar <- sqrt(sum(resultant^2))
    scores <- t(vapply(seq_len(nrow(x)), function(i) {
      psi_xi(x[i, ], theta$xi, d)
    }, numeric(d + 1L)))

    expect_equal(theta$mu, resultant / r_bar, tolerance = 1e-12)
    expect_equal(A_q(theta$kappa, d), r_bar, tolerance = 1e-10)
    expect_lt(max(abs(colSums(scores * probability_weights))), 1e-10)
  }
})

test_that("vMF score derivative agrees with the Bessel-ratio identity", {
  for (d in c(2L, 5L, 10L)) {
    for (kappa in c(d, 1.5 * d)) {
      A <- A_q(kappa, d)
      identity <- 1 - A^2 - d * A / kappa
      h <- 1e-5 * max(1, kappa)
      finite_difference <- (A_q(kappa + h, d) - A_q(kappa - h, d)) / (2 * h)
      expect_equal(A_q_prime(kappa, d), identity, tolerance = 1e-10)
      expect_equal(A_q_prime(kappa, d), finite_difference, tolerance = 1e-8)
    }
  }
})
