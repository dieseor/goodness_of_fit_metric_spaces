
library(rotasym)
library(sphunif)
library(testthat)
library(rgl)
library(viridis)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(progressr)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(tibble)

# Set working directory to save stuff
setwd("On the spherical cardioid distribution and its goodness-of-fit/")

# Load functions
source("cardioid-source.R")

## Checks

checks <- FALSE
if (checks) {

test_that("C_k^{(d-1)/2}(1) = d_{k,d} tau^{-1}_{k,d}", {

  for (p in 2:4) {
    d <- p - 1
    for (k in 0:4) {
      C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
      expect_equal(C1, sphunif::d_p_k(p = p, k = k) / tau_k_d(k = k, d = d))
    }
  }

})

test_that("Density integrates one", {

  for (p in 2:4) {
    for (r in seq(-1, 1, by = 0.5)) {
      for (k in 1:3) {
        expect_equal(int_sph_MC(d_sph_car, p = p, mu = c(rep(0, p - 1), 1),
                                rho = r, k = k, M = 1e4),
                     1, tolerance = 5e-2)
      }
    }
  }

})

test_that("Simulation methods are equivalent", {

  n <- 1e3
  for (rho in seq(-1, 1, by = 0.5)) {
    for (p in 2:4) {
      for (k in 1:3) {
        mu <- c(1, rep(0, p - 1))
        X_1 <- r_sph_car(n = n, mu = mu, rho = rho, k = k,
                         rejection = TRUE, odd_trick = FALSE, cubic = FALSE)
        X_2 <- r_sph_car(n = n, mu = mu, rho = rho, k = k,
                         rejection = FALSE, odd_trick = FALSE, cubic = FALSE)
        expect_gt(ks.test(X_1[, 1], X_2[, 1])$p.value, 0.01)
        expect_gt(ks.test(X_1[, 2], X_2[, 2])$p.value, 0.01)
        plot(ecdf(X_1[, 1]), main = paste("p =", p, ", rho =", rho, ", k =", k))
        lines(ecdf(X_2[, 1]), col = 2)
      }
    }
  }

})

test_that("Odd k simulation is correct", {

  n <- 1e4
  for (rho in seq(-1, 1, by = 0.5)) {
    for (p in 2:4) {
      for (k in c(1, 3, 5)) {
        mu <- c(1, rep(0, p - 1))
        X_1 <- r_sph_car(n = n, mu = mu, rho = rho, k = k,
                         rejection = FALSE, odd_trick = TRUE, cubic = TRUE)
        X_2 <- r_sph_car(n = n, mu = mu, rho = rho, k = k,
                         rejection = TRUE, odd_trick = FALSE, cubic = FALSE)
        expect_gt(ks.test(X_1[, 1], X_2[, 1])$p.value, 0.01)
        expect_gt(ks.test(X_1[, 2], X_2[, 2])$p.value, 0.01)
        plot(ecdf(X_1[, 1]), main = paste("p =", p, ", rho =", rho, ", k =", k))
        lines(ecdf(X_2[, 1]), col = 2)
      }
    }
  }

})

test_that("Convolution on the circle gives uniform density", {

  n <- 1e4
  rho_1 <- 0.5
  rho_2 <- -0.3
  k_1 <- 1
  k_2 <- 3
  mu_1 <- pi
  mu_2 <- pi / 2
  X_1 <- r_sph_car(n = n, mu = c(cos(mu_1), sin(mu_1)), rho = rho_1, k = k_1)
  X_2 <- r_sph_car(n = n, mu = c(cos(mu_2), sin(mu_2)), rho = rho_2, k = k_2)
  Theta_1 <- atan2(X_1[, 2], X_1[, 1])
  Theta_2 <- atan2(X_2[, 2], X_2[, 1])
  Theta <- (Theta_1 + Theta_2) %% (2 * pi)
  expect_gt(sphunif::unif_test(Theta, type = "Watson")$p.value, 0.1)

})

test_that("Convolution on the circle for non-uniform density", {

  n <- 1e4
  rho_1 <- 0.9
  rho_2 <- -0.9
  k_1 <- k_2 <- 3
  mu_1 <- pi
  mu_2 <- pi / 2
  X_1 <- r_sph_car(n = n, mu = c(cos(mu_1), sin(mu_1)), rho = rho_1, k = k_1)
  X_2 <- r_sph_car(n = n, mu = c(cos(mu_2), sin(mu_2)), rho = rho_2, k = k_2)
  X_3 <- r_sph_car(n = n, mu = c(cos(mu_1 + mu_2), sin(mu_1 + mu_2)),
                   rho = rho_1 * rho_2 / 2, k = k_1)
  Theta_1 <- atan2(X_1[, 2], X_1[, 1])
  Theta_2 <- atan2(X_2[, 2], X_2[, 1])
  Theta <- (Theta_1 + Theta_2) %% (2 * pi)
  Theta_3 <- atan2(X_3[, 2], X_3[, 1]) %% (2 * pi)
  expect_gt(ks.test(Theta, Theta_3)$p.value, 0.01)

})

test_that("Correct integration of projected density", {

  for (p in 3:5) {
    for (r in seq(-1, 1, by = 0.5)) {
      for (k in 1:3) {
        expect_equal(integrate(function(x) d_proj_car(x, rho = r, k = k, p = p),
                               lower = -1, upper = 1)$value,
                     1, tolerance = 1e-2)
      }
    }
  }

})

test_that("Projected quantile function inverts cdf", {

  for (p in 2:4) {
    for (r in seq(-1, 1, by = 0.5)) {
      for (k in 1:3) {
        x <- seq(-1, 1, l = 10)
        u <- seq(0 + 1e-15, 1 - 1e-15, l = 10)
        expect_equal(
          p_proj_car(x = q_proj_car(u = u, rho = r, k = k, p = p),
                     rho = r, k = k, p = p),
          u, tolerance = 1e-4)
        expect_equal(
          q_proj_car(u = p_proj_car(x = x, rho = r, k = k, p = p),
                     rho = r, k = k, p = p),
          x, tolerance = 1e-4)
      }
    }
  }

})

test_that("Correct integration of projection gamma", {

  for (p in 2:4) {
    for (r in seq(-1, 1, by = 0.5)) {
      for (k in 1:3) {
        mu <- c(rep(0, p - 1), 1)
        gamma <- c(1, rep(0, p - 1))
        expect_equal(integrate(function(x)
          d_proj_car_gamma(x, rho = r, k = k, p = p, mu = mu, gamma = gamma),
          lower = -1, upper = 1)$value,
          1, tolerance = 1e-2)
      }
    }
  }

})

test_that("Sampling and projecting on mu agrees with projected cdf", {

  n <- 1e3
  for (p in 2:4) {
    for (rho in seq(-1, 1, by = 0.5)) {
      for (k in 1:3) {
        mu <- c(rep(0, p - 1), 1)
        X <- r_sph_car(n = n, mu = mu, rho = rho, k = k, rejection = FALSE,
                       odd_trick = FALSE, cubic = FALSE)
        projs <- X %*% mu
        expect_gt(ks.test(projs, p_proj_car, rho = rho, k = k, p = p)$p.value,
                  0.01)
      }
    }
  }

})

test_that("Theoretical vs. empirical moments for k = m = 1", {

  set.seed(1234)
  M <- 5e4
  rho <- 0.5
  k <- 1
  for (p in 2:5) {
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    m1_emp <- mom_emp_sph_car(mu = mu, rho = rho, k = k, m = k, M = M)
    m1_theo <- mom_sph_car(mu = mu, rho = rho, k = k, m = k)
    expect_equal(m1_emp, m1_theo, tolerance = 0.01)
  }

})

test_that("Theoretical vs. empirical moments for k = m = 2", {

  set.seed(1234)
  M <- 5e4
  rho <- 0.5
  k <- 2
  for (p in 2:5) {
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    m1_emp <- mom_emp_sph_car(mu = mu, rho = rho, k = k, m = k, M = M)
    m1_theo <- mom_sph_car(mu = mu, rho = rho, k = k, m = k)
    expect_equal(m1_emp, m1_theo, tolerance = 0.01)
  }

})

test_that("Theoretical vs. empirical moments for k = m = 3", {

  set.seed(1234)
  M <- 5e4
  rho <- 0.75
  k <- 3
  for (p in 2:5) {
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    m3_emp <- mom_emp_sph_car(mu = mu, rho = rho, k = k, m = k, M = M)
    m3_theo <- mom_sph_car(mu = mu, rho = rho, k = k, m = k)
    expect_equal(m3_emp, m3_theo, tolerance = 0.01)
  }

})

test_that("Theoretical vs. empirical moments for k = m = 4", {

  set.seed(42)
  M <- 5e4
  rho <- 0.5
  k <- 4
  for (p in 3:5) {
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    m4_emp <- mom_emp_sph_car(mu = mu, rho = rho, k = k, m = k, M = M)
    m4_theo <- mom_sph_car(mu = mu, rho = rho, k = k, m = k)
    expect_equal(m4_emp, m4_theo, tolerance = 0.01)
  }

})

test_that("Theoretical vs. empirical variance for k = 1", {

  set.seed(1234)
  M <- 5e4
  rho <- 0.5
  for (p in 2:5) {
    mu <- c(rep(0, p - 1), 1)
    X <- r_sph_car(n = M, mu = mu, rho = rho, k = 1)
    var1_emp <- cov(X)
    var1_theo <- diag(rep(1 / p, p)) - (rho / p)^2 * tcrossprod(mu)
    expect_equal(var1_emp, var1_theo, tolerance = 0.01)
  }

})

test_that("MM estimator works for k = 1", {

  set.seed(1235)
  for (p in 2:5) {
    n <- 1e3
    rho <- 0.5
    k <- 1
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    X <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
    mm <- mm_sph_car(X = X, k = k)
    mm_mu <- mm_sph_car(X = X, k = k, mu = mu)
    expect_equal(mm$k, k)
    expect_equal(mm$rho, rho, tolerance = 0.1)
    expect_lt(1 - mm$mu %*% mu, 0.1)
    expect_equal(mm_mu$rho, rho, tolerance = 0.1)
  }

})

test_that("MM estimator works for k = 2", {

  set.seed(894)
  for (p in 2:5) {
    n <- 1e3
    rho <- 0.9
    k <- 2
    mu <- rotasym::r_unif_sphere(n = 2, p = p)[1, ]
    X <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
    mm <- mm_sph_car(X = X, k = k)
    mm_mu <- mm_sph_car(X = X, k = k, mu = mu)
    expect_equal(mm$k, k)
    expect_equal(mm$rho, rho, tolerance = 0.1)
    expect_lt(1 - abs(mm$mu %*% mu), 0.1)
    expect_equal(mm_mu$rho, rho, tolerance = 0.1)
  }

})

test_that("I_3(k) integral", {

  for (p in 3:5) {
    d <- p - 1
    for (k in 1:4) {
      I3_int <- integrate(function(t)
        Gegen_polyn(theta = acos(t), k = k, p = p)^3 *
        (1 - t^2)^(d / 2 - 1), lower = -1, upper = 1)$value
      I3_lhs <- (1 + 2 * k / (d - 1)) * I3_int /
        sphunif::Gegen_coefs(k = k, p = p, only_const = TRUE)
      I3_rhs <- (2 * k + d - 1)^2 / ((3 * k + d - 1) * (d - 1)) *
        factorial(k) / (factorial(k / 2)^3) *
        gamma((d + k - 1) / 2)^3 * gamma(d + 3 * k / 2 - 1) /
        (gamma(d + k - 1) * gamma((d - 1) / 2)^2 *
          gamma((d + 3 * k - 1) / 2))
      expect_equal(I3_lhs, ifelse(k %% 2 == 0, I3_rhs, 0), tolerance = 1e-6)
      if (k == 2) {
        expect_equal(I3_rhs, 2 * (d - 1) * (d + 3) / (d + 5), tolerance = 1e-6)
      }
    }
  }

})

test_that("Asymptotic variance for MM and GM coincide for k = 1,2", {

  rho <- 0.8
  for (p in 2:5) {
    expect_equal(1 / avar_mm_rho(rho = rho, k = 1, p = p),
                 1 / avar_gm_rho(rho = rho, k = 1, p = p))
    expect_equal(1 / avar_mm_rho(rho = rho, k = 2, p = p),
                 1 / avar_gm_rho(rho = rho, k = 2, p = p))
  }

})

test_that("Asymptotic variance for ML when d = 1 and k = 1, 2", {

  rho <- 0.6
  expect_equal(avar_mle_rho(rho = rho, k = 1, p = 2, exact = FALSE),
               avar_mle_rho(rho = rho, k = 1, p = 2, exact = TRUE))
  expect_equal(avar_mle_rho(rho = rho, k = 2, p = 2, exact = FALSE),
               avar_mle_rho(rho = rho, k = 2, p = 2, exact = TRUE))
  expect_equal(avar_mle_rho(rho = rho, k = 3, p = 2, exact = FALSE),
               avar_mle_rho(rho = rho, k = 3, p = 2, exact = TRUE))
  expect_equal(avar_mle_mu(rho = rho, k = 1, p = 2, exact = FALSE),
               avar_mle_mu(rho = rho, k = 1, p = 2, exact = TRUE))
  expect_equal(avar_mle_mu(rho = rho, k = 2, p = 2, exact = FALSE),
               avar_mle_mu(rho = rho, k = 2, p = 2, exact = TRUE))
  expect_equal(avar_mle_mu(rho = rho, k = 3, p = 2, exact = FALSE),
               avar_mle_mu(rho = rho, k = 3, p = 2, exact = TRUE))

})

test_that("Asymptotic variance for ML when d = 2 and k = 1, 2", {

  rho <- 0.7
  expect_equal(avar_mle_rho(rho = rho, k = 1, p = 3, exact = FALSE),
               avar_mle_rho(rho = rho, k = 1, p = 3, exact = TRUE))
  expect_equal(avar_mle_rho(rho = rho, k = 2, p = 3, exact = FALSE),
               avar_mle_rho(rho = rho, k = 2, p = 3, exact = TRUE))
  expect_equal(avar_mle_mu(rho = rho, k = 1, p = 3, exact = FALSE),
               avar_mle_mu(rho = rho, k = 1, p = 3, exact = TRUE))
  expect_equal(avar_mle_mu(rho = rho, k = 2, p = 3, exact = FALSE),
               avar_mle_mu(rho = rho, k = 2, p = 3, exact = TRUE))

})

test_that("Mgf expression works", {

  M <- 1e5
  for (p in 2:4) {
    for (k in 1:2) {
      for (rho in seq(-1, 1, by = 0.5)) {
        mu <- drop(rotasym::r_unif_sphere(n = 1, p = p))
        t <- matrix(rnorm(5 * p, sd = 2), nrow = 5)
        mgf_emp <- rowMeans(exp(t %*%
          t(r_sph_car(n = M, mu = mu, rho = rho, k = k))))
        mgf_theo <- sph_car_mgf(t = t, mu = mu, rho = rho, k = k)
        expect_equal(mgf_emp, mgf_theo, tolerance = 0.1)
      }
    }
  }

})

test_that("Cf expression works", {

  M <- 1e5
  for (p in 2:4) {
    for (k in 1:2) {
      for (rho in seq(-1, 1, by = 0.5)) {
        mu <- drop(rotasym::r_unif_sphere(n = 1, p = p))
        t <- matrix(rnorm(5 * p, sd = 2), nrow = 5)
        cf_emp <- rowMeans(exp(1.0i * t %*%
          t(r_sph_car(n = M, mu = mu, rho = rho, k = k))))
        cf_theo <- sph_car_cf(t = t, mu = mu, rho = rho, k = k)
        expect_equal(cf_emp, cf_theo, tolerance = 0.1)
      }
    }
  }

})

test_that("Gof statistic versions (Monte Carlo, weights, CvM) coincide
           for p = 2", {

  p <- 2
  n <- 3
  N <- 1e4
  mu <- c(rotasym::r_unif_sphere(n = 1, p = p))
  rho <- 0.75
  for (k in 1:3) {

    theta <- c(k, rho, mu)
    gammas <- rotasym::r_unif_sphere(n = N, p = p)
    X <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
    stat_mc <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_mc", gammas = gammas)
    stat_w <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_w", gammas = gammas)
    stat_cvm <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_cvm", gammas = gammas)
    expect_equal(stat_mc, stat_w)
    expect_equal(stat_mc, stat_cvm, tolerance = 1e-2)

  }

})

test_that("Gof statistic versions (Monte Carlo, weights, CvM) coincide
           for p = 3", {

  p <- 3
  n <- 3
  N <- 1e4
  mu <- c(rotasym::r_unif_sphere(n = 1, p = p))
  rho <- 0.15
  for (k in 1:2) {

    theta <- c(k, rho, mu)
    gammas <- rotasym::r_unif_sphere(n = N, p = p)
    X <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
    stat_mc <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_mc", gammas = gammas)
    stat_w <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_w", gammas = gammas)
    stat_cvm <- sph_car_stat_gof(X = X, theta = theta, weight = "CvM",
                                type_stat = "unif_cvm", gammas = gammas)
    expect_equal(stat_mc, stat_w)
    expect_equal(stat_mc, stat_cvm, tolerance = 1e-2)

  }

})

# p <- 5
# d <- p - 1
# rho <- 0.2
# k <- 3
# I3_int <- integrate(function(t) {
#   Gegen_polyn(theta = acos(t), k = k, p = p)^3 * (1 - t^2)^(d / 2 - 1)
#   }, lower = -1, upper = 1, abs.tol = 1e-10)$value
# I3_lhs <- (1 + 2 * k / (d - 1)) * I3_int /
#         sphunif::Gegen_coefs(k = k, p = p, only_const = TRUE)
# I3_rhs <- (2 * k + d - 1)^2 / ((3 * k + d - 1) * (d - 1)) *
#         factorial(k) / (factorial(k / 2)^3) *
#         gamma((d + k - 1) / 2)^3 * gamma(d + 3 * k / 2 - 1) /
#         (gamma(d + k - 1) * gamma((d - 1) / 2)^2 *
#           gamma((d + 3 * k - 1) / 2))
# X <- r_sph_car(n=1e5,k=k,mu=c(1,rep(0,p-1)),rho=rho,rejection=TRUE)
# (1 + 2 * k / (d - 1))^2 * var(drop(Gegen_polyn(theta=drop(acos(X[,1])),k=k,p=p)))
# d_p_k(p = p, k = k) - rho *  ifelse(k %% 2 == 0, I3_rhs, 0) - rho^2

# mean(drop(Gegen_polyn(theta = drop(acos(X[, 1])), k = k, p = p)))
# (1 + 2 * k / (d - 1))^(-1) * rho

# mean(drop(Gegen_polyn(theta = drop(acos(X[, 1])), k = k, p = p))^2)
# (1 + 2 * k / (d - 1))^(-2) * (d_p_k(p = p, k = k) + rho * I3_lhs)

# var(drop(Gegen_polyn(theta=drop(acos(X[,1])),k=k,p=p)))
# (1 + 2 * k / (d - 1))^(-2) * (d_p_k(p=p,k=k) + rho * I3_lhs-rho^2)

# var(stats_mm_rho_mu)
# (d_p_k(p=p,k=k) + rho * I3_lhs-rho^2)
}

## Plots for the paper

plots <- FALSE
if (plots) {

## Plot for d = 1
{

# Plot parameters
r_min <- 0.5
n_rug <- 200
tick_len <- 0.15
angle_pad <- 0.5

# Location mu
th_mu <- pi / 2
mu <- c(cos(th_mu), sin(th_mu))

# Grid of angles
th <- seq(0, 2 * pi, length.out = 400)
x <- cbind(cos(th), sin(th))

# Panels to display: combinations of (rho, k)
panels <- data.frame(
  rho = c(1, 1, -1, 1, 1, -1),
  k = c(1, 2, 2, 3, 4, 4)
)

# Explicit facet order to force the desired 2x3 layout:
# (k=1, rho=1) (k=2, rho=1) (k=2, rho=-1)
# (k=3, rho=1) (k=4, rho=1) (k=4, rho=-1)
panels$panel_raw <- paste0("k == ", panels$k, "*','~~rho == ", panels$rho)
panels$panel <- factor(
  panels$panel_raw,
  levels = c(
    "k == 1*','~~rho == 1",
    "k == 2*','~~rho == 1",
    "k == 2*','~~rho == -1",
    "k == 3*','~~rho == 1",
    "k == 4*','~~rho == 1",
    "k == 4*','~~rho == -1"
  )
)

# Evaluate the density on the angle grid for each panel
df <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  dens <- c(d_sph_car(x, mu = mu, rho = panels$rho[i], k = panels$k[i]))
  data.frame(th = th, dens = dens, panel = panels$panel[i])
}))

# Radius of the plot
r_max <- max(df$dens)
outer_r <- r_min + r_max

# Convert density to polar radius
df$r <- r_min + df$dens

# Simulate sample for rug points and color them by their density values
set.seed(1)
rug_df <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  xs <- r_sph_car(n = n_rug, mu = mu, rho = panels$rho[i], k = panels$k[i])
  th_s <- atan2(xs[, 2], xs[, 1])
  th_s <- ifelse(th_s < 0, th_s + 2 * pi, th_s)
  dens_s <- c(d_sph_car(xs, mu = mu, rho = panels$rho[i], k = panels$k[i]))
  data.frame(th = th_s, dens = dens_s, panel = panels$panel[i])
}))

# Rug ticks are drawn radially inward
rug_df$y <- r_min
rug_df$yend <- pmax(0, r_min - tick_len * r_max)

# Radial indicator at th_mu
mu_df <- panels
mu_df$panel <- panels$panel
mu_df$th <- th_mu
mu_df$y0 <- r_min
mu_df$dens_mu <- vapply(seq_len(nrow(panels)), function(i) {
  c(d_sph_car(matrix(mu, nrow = 1), mu = mu, rho = panels$rho[i],
              k = panels$k[i]))
}, numeric(1))
mu_df$y1 <- r_min + mu_df$dens_mu

# Angular grid lines every 2*pi/8
x_breaks <- seq(0, 2 * pi, by = 2 * pi / 8)

# Radial grid lines (shared across facets)
y_breaks <- r_min + pretty(c(0, r_max), n = 4)

# Labels for the radial rings
ring_lab_df <- data.frame(
  th = 5 * pi / 4,
  r = y_breaks,
  lab = formatC(y_breaks - r_min, format = "f", digits = 1)
)

# Angle labels for the main rays, placed *outside* the outer ring
angle_df <- data.frame(
  th = c(0, pi / 2, pi, 3 * pi / 2),
  r = outer_r + angle_pad * r_max,
  lab = c("0", "pi/2", "pi", "3*pi/2")
)

# Create the plot
plot_d1 <- ggplot(df, aes(x = th, y = r, group = panel, color = dens)) +
  geom_path(linewidth = 1) +
  geom_segment(
    data = rug_df,
    aes(x = th, xend = th, y = y, yend = yend, color = dens),
    inherit.aes = FALSE,
    lineend = "butt"
  ) +
  geom_segment(
    data = mu_df,
    aes(x = th, xend = th, y = y0, yend = y1),
    inherit.aes = FALSE,
    color = "red"
  ) +
  # Reference uniform density
  geom_hline(
    yintercept = r_min + 1 / (2 * pi),
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_point(
    data = mu_df,
    aes(x = th, y = y0),
    inherit.aes = FALSE,
    color = "red",
    size = 2.2
  ) +
  geom_point(
    data = mu_df,
    aes(x = th, y = y1),
    inherit.aes = FALSE,
    color = "red",
    size = 2.2
  ) +
  # Inner ring of the donut
  geom_hline(
    yintercept = r_min,
    linewidth = 0.4
  ) +
  # Ring height labels
  geom_text(
    data = ring_lab_df,
    aes(x = th, y = r, label = lab),
    inherit.aes = FALSE,
    angle = 0,
    hjust = 0.5,
    vjust = 0.5,
    size = 3
  ) +
  # Angle labels
  geom_text(
    data = angle_df,
    aes(x = th, y = r, label = lab),
    inherit.aes = FALSE,
    parse = TRUE,
    size = 3.4
  ) +
  facet_wrap(~panel, nrow = 2, ncol = 3, labeller = label_parsed) +
  # Rotate so that angle 0 is at the east and angles increase counterclockwise
  # clip = "off" allows drawing the angle labels outside the panel
  coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "off") +
  # Use x_breaks to create the angular grid lines; hide the x tick labels
  scale_x_continuous(limits = c(0, 2 * pi), breaks = x_breaks,
                     labels = rep("", length(x_breaks))) +
  scale_color_viridis_c(name = "Density", limits = c(0, r_max)) +
  # Keep the outer ring fixed by NOT extending the y limits beyond outer_r
  # oob_keep lets the angle labels be drawn outside without rescaling the rings
  scale_y_continuous(
    limits = c(0, outer_r),
    breaks = y_breaks,
    labels = rep("", length(y_breaks)),
    oob = scales::oob_keep
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(10, 22, 10, 22)
  )

# Increase colorbar size
plot_d1 <- plot_d1 +
  guides(
    color = guide_colorbar(
      title = "Density",
      barwidth = unit(10, "cm"),
      barheight = unit(0.5, "cm")
    )
  )

# Save the plot
ggsave("figs/dens_d1.pdf", plot = plot_d1, width = 12, height = 9)
ggsave("figs/dens_d1.png", plot = plot_d1, width = 12, height = 9)

}

## Plot for d = 2
{

# Number of samples
n <- 2000

# "Camera" parameters
azim <- 0
elev <- 0

# Rotate 3D points X (n x 3) by azimuth (around z) and elevation (around x),
# then project orthographically to (x, y)
rotate_xyz <- function(X, azim, elev) {
  Rz <- matrix(c(cos(azim), -sin(azim), 0,
                 sin(azim), cos(azim), 0,
                 0, 0, 1), 3, 3, byrow = TRUE)
  Rx <- matrix(c(1, 0, 0,
                 0, cos(elev), -sin(elev),
                 0, sin(elev), cos(elev)), 3, 3, byrow = TRUE)
  R <- Rx %*% Rz
  X %*% t(R)
}

# Location mu
th_mu <- pi / 2
ph_mu <- pi / 2
mu <- c(cos(th_mu) * sin(ph_mu), sin(th_mu) * sin(ph_mu), cos(ph_mu))

# Panels to display
panels <- data.frame(
  rho = c(1, 1, -1, 1, 1, -1),
  k = c(1, 2, 2, 3, 4, 4)
)

# Explicit facet order to force the desired 2x3 layout:
# (k=1, rho=1) (k=2, rho=1) (k=2, rho=-1)
# (k=3, rho=1) (k=4, rho=1) (k=4, rho=-1)
panels$panel_raw <- paste0("k == ", panels$k, "*','~~rho == ", panels$rho)
panel_levels <- c(
  "k == 1*','~~rho == 1",
  "k == 2*','~~rho == 1",
  "k == 2*','~~rho == -1",
  "k == 3*','~~rho == 1",
  "k == 4*','~~rho == 1",
  "k == 4*','~~rho == -1"
)

# Simulate points on the sphere for each panel, compute density for coloring,
# then rotate points to the chosen camera orientation
set.seed(1)
df_pts <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  xs <- r_sph_car(n = n, mu = mu, rho = panels$rho[i], k = panels$k[i])
  dens <- c(d_sph_car(xs, mu = mu, rho = panels$rho[i], k = panels$k[i]))
  xs_cam <- rotate_xyz(xs, azim = azim, elev = elev)
  data.frame(
    x = xs_cam[, 1],
    y = xs_cam[, 2],
    z = xs_cam[, 3],
    dens = dens,
    panel_raw = panels$panel_raw[i]
  )
}))

# Convert to an ordered factor for facetting
df_pts$panel <- factor(df_pts$panel_raw, levels = panel_levels)

# Common color scale across panels
dens_rng <- range(df_pts$dens)

# Depth cue: fade the back hemisphere by mapping rotated z to alpha
# (z < 0 gets lower alpha, z > 0 gets higher alpha)
df_pts$alpha <- ifelse(df_pts$z < 0, 0.15, 1)

# Red point for mu, rotated to the same camera orientation; repeated for each facet
mu_cam <- rotate_xyz(matrix(mu, nrow = 1), azim = azim, elev = elev)[1, ]
mu_df <- data.frame(
  x = mu_cam[1],
  y = mu_cam[2],
  panel = factor(panel_levels, levels = panel_levels)
)

# Outline of the orthographic projection of the sphere: the unit circle in (x,y)
t <- seq(0, 2 * pi, length.out = 361)
outline <- expand.grid(panel = factor(panel_levels, levels = panel_levels), t = t)
outline$x <- cos(outline$t)
outline$y <- sin(outline$t)

# Create the plot
plot_d2 <- ggplot(df_pts, aes(x, y)) +
  # Sphere boundary in the projection
  geom_path(
    data = outline,
    aes(x, y),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  # Samples, colored by density and alpha-faded by depth (rotated z)
  geom_point(
    aes(color = dens, alpha = alpha),
    size = 1.5
  ) +
  # Mark mu in red
  geom_point(
    data = mu_df,
    aes(x, y),
    inherit.aes = FALSE,
    color = "red",
    size = 3
  ) +
  # Facet layout with parsed plotmath labels
  facet_wrap(~panel, nrow = 2, ncol = 3, labeller = label_parsed) +
  # Keep the projection circular
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = FALSE) +
  scale_color_viridis_c(name = "Density", limits = dens_rng) +
  # Hide alpha legend
  scale_alpha(range = c(0.15, 1), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Increase colorbar size
plot_d2 <- plot_d2 +
  guides(
    color = guide_colorbar(
      title = "Density",
      barwidth = unit(10, "cm"),
      barheight = unit(0.5, "cm")
    )
  )

# Save the plot
ggsave("figs/dens_d2.pdf", plot = plot_d2, width = 12, height = 9)
ggsave("figs/dens_d2.png", plot = plot_d2, width = 12, height = 9)

}

## Plots AREs
{

# AREs of MM vs. MLE for different p's and k's

# Data frame of AREs
rhos <- seq(0.001, 0.999, by = 0.001)
ps <- 2:11
ks <- 1:2
df_a <- vector("list", length(ks) * length(ps))
idx <- 1
for (kk in ks) {
  for (pp in ps) {

    df_a[[idx]] <- tibble(
      rho = rhos,
      d = pp - 1,
      k = kk,
      Mu = avar_mle_mu(rho = rhos, k = kk, p = pp) /
        avar_mm_mu(rho = rhos, k = kk, p = pp),
      Rho = avar_mle_rho(rho = rhos, k = kk, p = pp) /
        avar_mm_rho(rho = rhos, k = kk, p = pp)
    ) |>
      pivot_longer(
        cols = c(Mu, Rho),
        names_to = "parameter",
        values_to = "ARE"
      )

    idx <- idx + 1

  }
}
df_a <- bind_rows(df_a) |>
  mutate(
    d = factor(d),
    k = factor(k, levels = ks, labels = paste0("k = ", ks)),
    parameter = factor(parameter, levels = c("Mu", "Rho"))
  )

# Create the plot
g_mm_mle <- ggplot(df_a, aes(x = rho, y = ARE, color = d, group = d)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(0, 1), color = "gray", linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(
    k ~ parameter,
    labeller = labeller(
      parameter = as_labeller(c(Mu = "mu", Rho = "rho"), label_parsed),
      k = label_value
    )
  ) +
  #scale_color_viridis_d(end = 0.95) +
  labs(
    x = expression(rho),
    y = "ARE",
    color = expression(d)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Save the plot
ggsave("figs/are_mm_mle.pdf", plot = g_mm_mle, width = 9, height = 9)
ggsave("figs/are_mm_mle.png", plot = g_mm_mle, width = 9, height = 9)

# AREs of MLE vs. MM-mu: fix p = 3, vary k = 3:10

# Data frame of AREs
rhos <- seq(0.001, 0.999, by = 0.001)
p_fix <- 3
ks <- c(1, 2, 3, 4, 5, 6)
df_k <- lapply(ks, function(kk) {
  tibble(
    rho = rhos,
    k = kk,
    ARE = avar_mle_rho(rho = rhos, k = kk, p = p_fix) /
      avar_gm_rho(rho = rhos, k = kk, p = p_fix)
  )
}) |>
  bind_rows() |>
  mutate(k = factor(k))

# Create the plot
g_k <- ggplot(df_k, aes(x = rho, y = ARE, color = k, group = k)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(0, 1), color = "gray", linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  #scale_color_viridis_d(end = 0.95) +
  labs(
    x = expression(rho),
    y = "ARE",
    color = "k"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )
g_k

# Save the plot
ggsave("figs/are_mmmu_mle_k.pdf", plot = g_k, width = 7, height = 7)
ggsave("figs/are_mmmu_mle_k.png", plot = g_k, width = 7, height = 7)

# AREs of MLE vs. MM-mu: fix k = 3, vary p = 3:10

# Data frame of AREs
rhos <- seq(0.001, 0.999, by = 0.001)
k_fix <- 3
ps <- 2:11
df_d <- lapply(ps, function(pp) {
  tibble(
    rho = rhos,
    d = pp - 1,
    ARE = avar_mle_rho(rho = rhos, k = k_fix, p = pp) /
      avar_gm_rho(rho = rhos, k = k_fix, p = pp)
  )
}) |>
  bind_rows() |>
  mutate(d = factor(d))

# Create the plot
g_d <- ggplot(df_d, aes(x = rho, y = ARE, color = d, group = d)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(0, 1), color = "gray", linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  #scale_color_viridis_d(end = 0.95) +
  labs(
    x = expression(rho),
    y = "ARE",
    color = expression(d)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )
g_d

# Save the plot
ggsave("figs/are_mmmu_mle_d.pdf", plot = g_d, width = 7, height = 7)
ggsave("figs/are_mmmu_mle_d.png", plot = g_d, width = 7, height = 7)

}

## Plots asymptotic distributions of MM and MLE estimators
{

# Parameters
p <- 3
d <- p - 1
k <- 2
mu <- c(rep(0, p - 1), 1)
rho <- 0.5
n <- 1e3
M <- 1e4

# Simulation
set.seed(42)
mm_rho <- mm_rho_mu <- mle_rho <- numeric(M)
mm_mu <- mle_mu <- matrix(0, nrow = M, ncol = p)
pb <- txtProgressBar(min = 0, max = M, style = 3)
for (m in 1:M) {

  # Sample
  setTxtProgressBar(pb, m)
  X <- r_sph_car(n = n, mu = mu, rho = rho, k = k)

  ## MM mu
  mm_mu_est <- mm_sph_car(X = X, k = k, mu = mu, rho_trunc = FALSE)
  mm_rho_mu[m] <- mm_mu_est$rho

  ## MM
  mm <- mm_sph_car(X = X, k = k, rho_trunc = FALSE, rho_sign = 1)
  mm_rho[m] <- mm$rho
  mm_mu[m, ] <- mm$mu

  ## MLE
  mle <- mle_sph_car(X = X, k = k, rho0 = rho, mu0 = mu)
  mle_rho[m] <- mle$rho
  mle_mu[m, ] <- mle$mu

}

# CLT stats
stats_mm_rho <- sqrt(n) * (mm_rho - rho)
stats_mm_mu <- sqrt(n) * (mm_mu - matrix(rep(mu, M), nrow = M, byrow = TRUE))
stats_mle_rho <- sqrt(n) * (mle_rho - rho)
stats_mle_mu <- sqrt(n) * (mle_mu - matrix(rep(mu, M), nrow = M, byrow = TRUE))
stats_mm_rho_mu <- sqrt(n) * (mm_rho_mu - rho)

# Asymptotic variances
asymp_var_rho_mm <- avar_mm_rho(rho = rho, k = k, p = p)
asymp_var_mu_mm <- (1 - mu[2]^2) * avar_mm_mu(rho = rho, k = k, p = p)
asymp_var_rho_mle <- avar_mle_rho(rho = rho, k = k, p = p)
asymp_var_mu_mle <- (1 - mu[2]^2) * avar_mle_mu(rho = rho, k = k, p = p)
asymp_var_rho_mu_mm <- avar_gm_rho(rho = rho, k = k, p = p)

# Plot histogram + N(0, sd) density
hist_vs_norm <- function(x, sd, bins = 50) {

  x <- as.numeric(x)
  ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = bins,
      fill = "gray95",
      color = "gray60",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = mean(x),
      color = "gray60",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_vline(
      xintercept = 0,
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = 0, sd = sd),
      color = "red",
      linewidth = 1
    ) +
    theme_minimal() +
    labs(x = "Statistics", y = "Density") +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_continuous(limits = c(0, 0.3)) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92")
    ) +
    guides(fill = "none", color = "none", linetype = "none")

}
p_mm_rho <- hist_vs_norm(x = stats_mm_rho, sd = sqrt(asymp_var_rho_mm))
p_mm_mu <- hist_vs_norm(x = stats_mm_mu[, 1], sd = sqrt(asymp_var_mu_mm))
p_mle_rho <- hist_vs_norm(x = stats_mle_rho, sd = sqrt(asymp_var_rho_mle))
p_mle_mu <- hist_vs_norm(x = stats_mle_mu[, 1], sd = sqrt(asymp_var_mu_mle))

# Save plots
ggsave(paste0("figs/stats_mm_rho_k", k, ".pdf"), plot = p_mm_rho, width = 7, height = 7)
ggsave(paste0("figs/stats_mm_mu_k", k, ".pdf"), plot = p_mm_mu, width = 7, height = 7)
ggsave(paste0("figs/stats_mle_rho_k", k, ".pdf"), plot = p_mle_rho, width = 7, height = 7)
ggsave(paste0("figs/stats_mle_mu_k", k, ".pdf"), plot = p_mle_mu, width = 7, height = 7)

# Match standard deviations: mu (<2% error)
abs(sqrt(asymp_var_mu_mle) - sd(stats_mle_mu[, 1])) / sqrt(asymp_var_mu_mle) * 100
abs(sqrt(asymp_var_mu_mm) - sd(stats_mm_mu[, 1])) / sqrt(asymp_var_mu_mm) * 100

# Match standard deviations: rho (<2% error)
abs(sqrt(asymp_var_rho_mle) - sd(stats_mle_rho)) / sqrt(asymp_var_rho_mle) * 100
abs(sqrt(asymp_var_rho_mm) - sd(stats_mm_rho)) / sqrt(asymp_var_rho_mm) * 100

# Match standard deviations: rho (MM-mu) (<1% error)
abs(sqrt(asymp_var_rho_mu_mm) - sd(stats_mm_rho_mu)) / sqrt(asymp_var_rho_mu_mm) * 100

# Check means (should be close to 0)
mean(stats_mle_mu[, 1])
mean(stats_mm_mu[, 1])
mean(stats_mm_rho) # +0.20 for k=1, +0.33 for k=2
mean(stats_mle_rho)# +0.20 for k=1, +0.33 for k=2
mean(stats_mm_rho_mu)

}

}

## Simulation study for GOF test
{

# Warp-speed Monte Carlo function
warp_speed_mc <- function(r_mod, n, k, M, K, B, type_stat, weight, est_type,
                          mu0 = NULL, rho0 = NULL, rho_sign = NULL) {

  # Monte Carlo study
  progressr::with_progress({
    prog <- progressr::progressor(along = seq_len(M))
    tests <- foreach(j = 1:M, .inorder = TRUE, .packages = "sphunif") %dorng% {

      # Show progress
      if (requireNamespace("progressr", quietly = TRUE)) prog()

      # # Fix seed for reproducibility
      # set.seed(j, kind = "Mersenne-Twister")

      # Sample
      samp <- r_mod(n = n)

      # Tests
      sph_car_gof(x = samp, k = k, K = K, B = 1,
                  type_stat = type_stat, weight = weight, est_type = est_type,
                  mu0 = mu0, rho0 = rho0, rho_sign = rho_sign,
                  plot_boot = FALSE, parallel = FALSE)

    }
  })

  # Extract statistics
  stats <- rbind(sapply(tests, function(x) x$statistic))
  stats_boot <- rbind(sapply(tests, function(x) x$statistic_boot))

  # Compute powers
  alphas <- seq(0, 1, l = 101)
  pow <- t(rbind(sapply(alphas, function(alph) {

      q_boot_alph_null_up <- apply(stats_boot, 1, quantile, probs = 1 - alph)
      return(rowMeans(stats > q_boot_alph_null_up))

  })))
  pow <- cbind(alphas, pow)
  colnames(pow) <- c("alpha", names(tests[[1]]$statistic))
  pow <- as.data.frame(pow)

  # Return results
  return(list(stats = stats, stats_boot = stats_boot, pow = pow))

}

# Regular Monte Carlo
regular_mc <- function(r_mod, n, k, M, K, B, type_stat, weight, est_type,
                       mu0 = NULL, rho0 = NULL, rho_sign = NULL, ...) {

  # Monte Carlo study
  progressr::with_progress({
    prog <- progressr::progressor(along = seq_len(M))
    tests <- foreach(j = 1:M, .inorder = TRUE, .packages = "sphunif") %dorng% {

      # Show progress
      if (requireNamespace("progressr", quietly = TRUE)) prog()

      # # Fix seed for reproducibility
      # set.seed(j, kind = "Mersenne-Twister")

      # Sample
      samp <- r_mod(n = n)

      # Tests
      sph_car_gof(x = samp, k = k, K = K, B = B,
                  type_stat = type_stat, weight = weight, est_type = est_type,
                  mu0 = mu0, rho0 = rho0, rho_sign = rho_sign,
                  plot_boot = FALSE, parallel = FALSE)

    }
  })

  # Extract pvalues
  pvals <- rbind(sapply(tests, function(x) x$pvalue))

  # Compute powers
  alphas <- seq(0, 1, l = 101)
  pow <- t(rbind(sapply(alphas, function(alph) {
      return(rowMeans(pvals < alph))
  })))
  pow <- cbind(alphas, pow)
  colnames(pow) <- c("alpha", names(tests[[1]]$statistic))
  pow <- as.data.frame(pow)

  # Return results
  return(list(pvals = pvals, pow = pow))

}

# Run simulations
sim_gof <- function(d, n, k, M, K, B, sim, sim_par, sim_par_others,
                    type_stat, weight, est_type = "mm",
                    mu_start = NULL, rho_start = NULL, rho_sign = 1,
                    mc_type = c("regular_mc", "warp_speed_mc")[1]) {

  # Type of Monte Carlo
  mc <- get(mc_type)
  B <- ifelse(mc_type == "regular_mc", B, 1)

  # Run simulations
  simus <- list()
  for (j in seq_along(sim_par)) {

    cat(mc_type, "for sim_par =", sim_par[j], "\n")
    sim_par_j <- function(n) {

      sim(n = n, sim_par = sim_par[j], sim_par_others = sim_par_others)

    }
    simus[[j]] <- mc(r_mod = sim_par_j, k = k,
                     n = n, M = M, K = K, B = B,
                     type_stat = type_stat, weight = weight,
                     # est_type = "mle",
                     # mu0 = sim_par_others$mu,
                     # rho0 = sim_par[j],
                     est_type = est_type,
                     mu0 = mu_start,
                     rho0 = rho_start,
                     rho_sign = rho_sign)
    simus[[j]]$sim_par <- sim_par[j]

  }

  # Extract powers
  alphas <- c(0.01, 0.05, 0.10)
  pow <- lapply(simus, function(simus_j) {
    simus_j$pow |>
      filter(alpha %in% alphas) |>
      mutate(sim_par = simus_j$sim_par) |>
      select(sim_par, everything())
  })
  pow <- do.call(rbind, pow)
  return(as.data.frame(pow))

}

# Progress bar
handlers(handler_progress(
  format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
  clear = FALSE))

# Monte Carlo
# doRNG::registerDoRNG()
doFuture::registerDoFuture()
future::plan(future::multisession(), workers = 14)

# Common parameters
M <- 1000
B <- 1000
K <- 100
type_stat <- c("unif_cvm", "unif_mc", "Pn", "Pn", "C_k", "C_k")
weight <- c("CvM", "AD", "CvM", "AD", "CvM", "AD")
# type_stat <- c("unif_cvm", "Pn", "C_k")
# weight <- c("CvM", "CvM", "CvM")
# type_stat <- c("unif_cvm")
# weight <- c("CvM")

# Binomial prediction intervals
predint_binom <- function(n, p, level) {
  a <- 1 - level
  c(lower = qbinom(a / 2, n, p), upper = qbinom(1 - a / 2, n, p))
}

## Null table

n <- 100
rhos_0 <- c(0.25, 0.5, 0.75)
r_mod_0 <- function(n, mu, rho, k) {
  r_sph_car(n = n, mu = mu, rho = rho, k = k)
}
sim_0 <- function(n, sim_par, sim_par_others) {
  r_sph_car(n = n, rho = sim_par, mu = sim_par_others$mu, k = sim_par_others$k)
}

# k = 1, d = 1, 2
set.seed(2026)
res0_d1_k1_mm <- sim_gof(d = 1, n = n, k = 1, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_0,
                         sim_par_others = list(k = 1, d = 1, mu = c(1, 0)),
                         type_stat = type_stat, weight = weight)
res0_d2_k1_mm <- sim_gof(d = 2, n = n, k = 1, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_0,
                         sim_par_others = list(k = 1, d = 2, mu = c(1, 0, 0)),
                         type_stat = type_stat, weight = weight)

# k = 2, d = 1, 2
res0_d1_k2_mm <- sim_gof(d = 1, n = n, k = 2, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_0,
                         sim_par_others = list(k = 2, d = 1, mu = c(1, 0)),
                         type_stat = type_stat, weight = weight)
res0_d2_k2_mm <- sim_gof(d = 2, n = n, k = 2, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_0,
                         sim_par_others = list(k = 2, d = 2, mu = c(1, 0, 0)),
                         type_stat = type_stat, weight = weight)

# Table
alpha_table <- 0.05
res0_table <- bind_rows(
  mutate(k = 1, d = 1, res0_d1_k1_mm),
  mutate(k = 1, d = 2, res0_d2_k1_mm),
  mutate(k = 2, d = 1, res0_d1_k2_mm),
  mutate(k = 2, d = 2, res0_d2_k2_mm)
) |>
  filter(alpha == alpha_table) |>
  mutate(across(.cols = -c(d, k, sim_par), .fns = ~ . * 100)) |>
  select(k, d, sim_par, everything(), -alpha)
res0_table

# Check if counts are within 95%-prediction intervals
pred_set <- predint_binom(n = M, p = alpha_table, level = 0.95)
check_counts_res0 <- res0_table |>
  select(-d, -k, -sim_par) |>
  mutate(across(everything(), ~ . / 100 * M)) |>
  mutate(across(everything(), ~
                  pred_set["lower"] <= . & . <= pred_set["upper"]))
check_counts_res0

# Mark counts outside prediction intervals
res0_table_checked <- res0_table
for (col in 4:ncol(res0_table_checked)) {
  res0_table_checked[[col]] <- ifelse(
    check_counts_res0[[col - 3]],
    res0_table_checked[[col]],
    paste0("\\underline{", res0_table_checked[[col]], "}")
  )
}
res0_table_checked

# LaTeX table
print(xtable::xtable(res0_table_checked, digits = 2),
      include.rownames = FALSE, sanitize.text.function = identity,
      booktabs = TRUE)

## Alternative table

# k = 1 (null) vs. k = 2 (reality), d = 1, 2
n <- 100
rhos_1 <- c(0.1, 0.25, 0.50)
set.seed(2026)
res1_d1_k1_mm <- sim_gof(d = 1, n = n, k = 1, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_1,
                         sim_par_others = list(k = 2, d = 1, mu = c(1, 0)),
                         type_stat = type_stat, weight = weight)
res1_d2_k1_mm <- sim_gof(d = 2, n = n, k = 1, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_1,
                         sim_par_others = list(k = 2, d = 1, mu = c(1, 0)),
                         type_stat = type_stat, weight = weight)

# k = 2 (null) vs. k = 1 (reality), d = 1, 2
res1_d1_k2_mm <- sim_gof(d = 1, n = n, k = 2, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_1,
                         sim_par_others = list(k = 1, d = 1, mu = c(1, 0)),
                         type_stat = type_stat, weight = weight)
res1_d2_k2_mm <- sim_gof(d = 2, n = n, k = 2, M = M, K = K, B = B,
                         sim = sim_0, sim_par = rhos_1,
                         sim_par_others = list(k = 1, d = 2, mu = c(1, 0, 0)),
                         type_stat = type_stat, weight = weight)

# Table
alpha_table <- 0.05
res1_table <- bind_rows(
  mutate(k = 1, d = 1, res1_d1_k1_mm),
  mutate(k = 1, d = 2, res1_d2_k1_mm),
  mutate(k = 2, d = 1, res1_d1_k2_mm),
  mutate(k = 2, d = 2, res1_d2_k2_mm)
) |>
  filter(alpha == alpha_table) |>
  mutate(across(.cols = -c(d, k, sim_par), .fns = ~ . * 100)) |>
  select(k, d, sim_par, everything(), -alpha)
res1_table

# LaTeX table
print(xtable::xtable(res1_table, digits = 1), booktabs = TRUE,
      include.rownames = FALSE)

# Save results
save(list = ls(pattern = "res0_*|res1_*"), file = "gof_simus.RData")

}

## Fisher B6 set data application
{

th_fisher <- circular::fisherB6$set1 / 360 * 2 * pi
x_fisher <- cbind(sin(th_fisher), cos(th_fisher))
sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "unif_cvm", weight = "CvM")

sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "CvM")
sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "AD")
sph_car_gof(x = x_fisher, k = 2, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "CvM")
sph_car_gof(x = x_fisher, k = 2, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "AD")
sph_car_gof(x = x_fisher, k = 3, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "CvM")
sph_car_gof(x = x_fisher, k = 3, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "AD")

with_progress(
sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "unif_mc", weight = "CvM")
)#51s

with_progress(
sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "unif_cvm", weight = "CvM")
)#1.58s

with_progress(
sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
            type_stat = "Pn", weight = "CvM")
)#2.56s

with_progress(
  sph_car_gof(x = x_fisher, k = 1, K = 1e3, B = 1e3, plot_boot = TRUE,
             type_stat = c("Pn", "Pn"), weight = c("CvM", "AD"))
)

}

## Distribution of F_gamma(gamma'X_i)
{

p <- 3
rho <- 0.5
k <- 1
n <- 1e4
gamma <- r_unif_sphere(n = n, p = p)
mu <- c(1, rep(0, p - 1))
X <- c(-sqrt(2) / 2, sqrt(2) / 2, rep(0, p - 2))
F_gamma <- apply(gamma, 1, function(g)
  drop(p_proj_car_gamma(x = drop(g %*% X), rho = rho, k = k,
                        p = p, mu = mu, gamma = g))
)
hist(F_gamma)
}

## Several uniform projections
{

g_tau <- function(uv, p, tau) {

  stopifnot(is.matrix(uv))
  stopifnot(ncol(uv) == 2)
  stopifnot(abs(tau) < 1)
  u <- uv[, 1]
  v <- uv[, 2]
  dens <- 1 - tau^2 - u^2 - v^2 + 2 * tau * u * v
  rotasym::w_p(p = p - 2) / rotasym::w_p(p = p) *
    (1 - tau^2)^(-(p - 3) / 2) *
    ifelse(dens > 0, dens^((p - 4) / 2), 0)

}

g_tau_cond <- function(u, v, p, tau) {

  stopifnot(length(v) == 1)
  g_tau(uv = cbind(u, v), p = p, tau = tau) / drop(d_proj_unif(x = v, p = p))

}

p <- 5
tau <- -0.5
u <- seq(-1, 1, l = 100)
sdetorus::plotSurface2D(x = u, y = u,
                        f = function(x) g_tau(uv = x, p = p, tau = tau),
                        fVect = TRUE)
x <- r_unif_sph(n = 1e4, p = p)[, , 1]
points(x %*% c(1, rep(0, p - 1)),
       x %*% c(tau, rep(sqrt(1 - tau^2) / sqrt(p), p - 1)))
plot(u, g_tau_cond(u = u, v = 0.5, p = p, tau = tau), type = "l")

}

## Sample + visualize densities
{

rho <- 1
th_mu <- pi / 2
mu <- c(cos(th_mu), sin(th_mu))
th <- seq(0, 2 * pi, length.out = 100)
x <- cbind(cos(th), sin(th))
par(mfcol = c(2, 3))
for (k in c(1, 3)) {
  dens <- c(d_sph_car(x, mu = mu, rho = rho, k = k))
  plot(th, dens, type = "l", ylim = c(0, max(dens)))
  abline(v = th_mu, col = 2, pch = 19, cex = 2)
}
for (rho in c(1, -1)) {
  for (k in c(2, 4)) {
    dens <- c(d_sph_car(x, mu = mu, rho = rho, k = k))
    plot(th, dens, type = "l", ylim = c(0, max(dens)))
    abline(v = th_mu, col = 2, pch = 19, cex = 2)
  }
}

n <- 1e3
rho <- 1
mu <- c(0, 0, 1)
mfrow3d(2, 3, byrow = TRUE, sharedMouse = TRUE)
for (k in 1:6) {
  samp <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
  dens <- d_sph_car(samp, mu = mu, rho = rho, k = k)
  plot3d(samp, col = viridis(nrow(samp))[rank(dens)], size = 5)
  points3d(mu, col = 2, size = 15)
}

n <- 1e3
rho <- 1
mu <- c(0, 0, 1)
mfrow3d(2, 3, byrow = TRUE, sharedMouse = TRUE)
for (rho in c(-1, 1)) {
  for (k in c(2, 4, 6)) {
    samp <- r_sph_car(n = n, mu = mu, rho = rho, k = k)
    dens <- d_sph_car(samp, mu = mu, rho = rho, k = k)
    plot3d(samp, col = viridis(nrow(samp))[rank(dens)], size = 5)
    points3d(mu, col = 2, size = 15)
  }
}

}

## Projection for gamma
{

p <- 4
n <- 1e4
rho <- 1
mu <- c(rep(0, p - 1), 1)
k <- 3
samp <- r_sph_car(n = n, mu = mu, rho = rho, k = k)

# gamma <- mu
gamma <- c(sqrt(2) / 2, -sqrt(2) / 2, rep(0, p - 2))
proj_samp <- samp %*% gamma
ks.test(proj_samp, p_proj_car_gamma,
        rho = rho, k = k, p = p, mu = mu, gamma = gamma)

hist(proj_samp, probability = TRUE, breaks = 50, ylim = c(0, 1))
lines(density(proj_samp), col = 1, lwd = 2)
# curve(d_proj_car(x, rho = rho, k = k, p = p),
#       add = TRUE, col = 3, lwd = 2)
curve(d_proj_car_gamma(x, rho = rho, k = k, p = p, mu = mu, gamma = gamma),
      add = TRUE, col = 2, lwd = 2)

plot(ecdf(proj_samp))
# curve(p_proj_car(x, rho = rho, k = k, p = p),
#       from = -1, to = 1, col = 3, lwd = 2, add = TRUE)
curve(p_proj_car_gamma(x, rho = rho, k = k, p = p, mu = mu, gamma = gamma),
      from = -1, to = 1, col = 2, lwd = 2, add = TRUE)

}

## Null moments
{

p <- 3
kron_pow <- function(x, r) {

  stopifnot(is.vector(x))
  if (r == 1) {

    return(x)

  } else if (r == 2) {

    return(x %x% x)

  } else {

    return(x %x% kron_pow(x, r - 1))

  }

}
x_pow <- function(x, r) t(apply(x, 1, kron_pow, r = r))

x1_pow <- colMeans(x_pow(samp, r = 1))
x2_pow <- colMeans(x_pow(samp, r = 2))
x3_pow <- colMeans(x_pow(samp, r = 3))
x4_pow <- colMeans(x_pow(samp, r = 4))
x5_pow <- colMeans(x_pow(samp, r = 5))
x6_pow <- colMeans(x_pow(samp, r = 6))
x7_pow <- colMeans(x_pow(samp, r = 7))
x8_pow <- colMeans(x_pow(samp, r = 8))
x9_pow <- colMeans(x_pow(samp, r = 9))
x2_array <- array(x2_pow, rep(p, 2))
x3_array <- array(x3_pow, rep(p, 3))
x4_array <- array(x4_pow, rep(p, 4))
x5_array <- array(x5_pow, rep(p, 5))
x6_array <- array(x6_pow, rep(p, 6))
x7_array <- array(x7_pow, rep(p, 7))
x8_array <- array(x8_pow, rep(p, 8))
x9_array <- array(x9_pow, rep(p, 9))

sapply(1:p, function(i) x2_array[rbind(rep(i, 2))])
sapply(1:p, function(i) x3_array[rbind(rep(i, 3))])
sapply(1:p, function(i) x4_array[rbind(rep(i, 4))])
sapply(1:p, function(i) x5_array[rbind(rep(i, 5))])
sapply(1:p, function(i) x6_array[rbind(rep(i, 6))])
sapply(1:p, function(i) x7_array[rbind(rep(i, 7))])
sapply(1:p, function(i) x8_array[rbind(rep(i, 8))])
sapply(1:p, function(i) x9_array[rbind(rep(i, 9))])

sum(x3_pow^2)
sum(x5_pow^2)
sum(x7_pow^2)
sum(x9_pow^2)

sum(x2_pow^2); 1/p
sum(x4_pow^2); 1/(p+2)
sum(x6_pow^2); 1/(p+4)
sum(x8_pow^2); 1/(p+6)

# Cardioid distributions
# Legendre polynomial distribution: in Ebner et al

# AR

f_plus_2 <- function(x, mu, rho) {
  d <- length(mu) - 1
  s <- colSums(t(x^2) * (mu^2))
  (2^(d + 1) + 2^d * (d - 1) * ((d + 1) * s - 1)) / w_p(p = d + 1)
}

mu <- c(1, 2, 3)
x <- abs(r_unif_sphere(n = 1e4, p = 3))
plot3d(x, col = viridis(nrow(x))[
  rank(f_plus_2(x, mu = mu/sqrt(sum(mu^2)), rho = 0.5))],
  size = 5)

mean(f_plus_2(x, mu = mu/sqrt(sum(mu^2)), rho = 1)) * w_p(p = 3) / 8

}

## Detect non-uniformity of cardioid
{

library(sphunif)
r_H1 <- function(n, p, M, rho, k) {
  samp <- array(dim = c(n, p, M))
  for (j in 1:M) {
    samp[, , j] <- r_sph_car(n, mu = c(1, rep(0, p - 1)), rho = rho, k = k,
                             rejection = FALSE, odd_trick = TRUE, cubic = TRUE)
  }
  return(samp)
}


p <- 3
n <- 500
stats_0 <- unif_stat_MC(n = n, p = p, M = 1000,
                        type = c("PAD", "Gine_Fn", "Pycke", "Riesz", "Poisson"))

k <- 11
unif_stat_MC(n = n, p = p, M = 1000, r_H1 = r_H1,
             type = c("PAD", "Gine_Fn", "Pycke", "Riesz", "Poisson"),
             rho = 1, k = k, crit_val = stats_0$crit_val_MC)$power_MC

}
