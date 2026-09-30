# Cardioid routines by Eduardo García-Portugués, adapted with permission.
# Retained only where used by the distance-profile GOF implementation.

d_sph_car <- function(x, mu, rho, k, log = FALSE) {

  # Checks
  p <- length(mu)
  stopifnot(length(mu) == ncol(rbind(x)))
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)

  C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
  log_dens <-
    log1p(rho * Gegen_polyn(theta = acos(x %*% mu), k = k, p = p) / C1) -
    w_p(p = p, log = TRUE)
  if (log) {

    return(log_dens)

  } else {

    return(exp(log_dens))

  }

}
r_sph_car <- function(n, mu, rho, k, odd_trick = TRUE, cubic = TRUE) {

  p <- length(mu)
  if (odd_trick && k %% 2 == 1) {

    r_U <- function(n) r_unif_sphere(n = n, p = p - 1)
    r_V <- function(n) {

      # Draw from the symmetric base density
      x <- r_proj_unif(n = n, p = p)
      r <- abs(x)

      # Sign-selection probability
      q  <- drop(sphunif::Gegen_polyn(theta = acos(r), k = k, p = p)) /
        drop(sphunif::Gegen_polyn(theta = 0, k = k, p = p))
      p_plus <- (1 + rho * q) / 2
      p_plus <- pmin(pmax(p_plus, 0), 1)
      s <- ifelse(runif(n) <= p_plus, 1, -1)
      return(s * r)

    }
    samp <- r_tang_norm(n = n, theta = mu, r_U = r_U, r_V = r_V)

  } else if (cubic && p == 3 && k == 2 && rho > 0) {

    r_U <- function(n) r_unif_sphere(n = n, p = p - 1)
    r_V <- function(n) {

      u <- runif(n)
      q3 <- 2 * (1 - 2 * u) / rho
      p3 <- (2 - rho) / rho
      delta <- (q3 / 2)^2 + (p3 / 3)^3
      sqrt_1 <- -q3 / 2 + sqrt(delta)
      sqrt_2 <- -q3 / 2 - sqrt(delta)
      return(sign(sqrt_1) * abs(sqrt_1)^(1 / 3) +
               sign(sqrt_2) * abs(sqrt_2)^(1 / 3))

    }
    samp <- r_tang_norm(n = n, theta = mu, r_U = r_U, r_V = r_V)

  } else {

    M <- 1 + abs(rho)
    X <- rotasym::r_unif_sphere(n = ceiling(1.25 * M * n + 10), p = p)
    U <- runif(ceiling(1.25 * M * n + 10))
    C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
    dens <-
      (1 + rho * drop(Gegen_polyn(theta = acos(X %*% mu), k = k, p = p)) / C1)
    accept <- U <= dens / M
    samp <- X[accept, , drop = FALSE]
    samp <- samp[seq_len(min(n, nrow(samp))), , drop = FALSE]
    if (nrow(samp) < n) {

      message(paste("Rejection sampling: only",
                    nrow(samp), "samples generated out of", n,
                    "required. Generating extra samples."))
      samp_extra <- r_sph_car(n = n - nrow(samp), mu = mu, rho = rho,
                              k = k, odd_trick = odd_trick, cubic = cubic)
      samp <- rbind(samp, samp_extra)

    }

  }
  return(samp)

}
p_proj_car_gamma <- function(x, rho, k, p, mu, gamma) {

  # Checks
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)
  stopifnot(p >= 2)
  stopifnot(length(mu) == length(gamma))

  gamma_mu <- max(min(sum(gamma * mu), 1), -1)
  C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
  cdf_1 <- drop(p_proj_unif(x = x, p = p))
  lambda_k <- rho * w_p(p = p - 1) / w_p(p = p) *
    drop(Gegen_polyn(theta = acos(gamma_mu), k = k, p = p)) / C1^2
  if (p == 2) {

    G_k <- sin(k * acos(x)) / k

  } else {

    G_k <- (p - 2) / (k * (k + p - 2)) *
      drop(Gegen_polyn(theta = acos(x), k = k - 1, p = p + 2)) *
      (1 - x^2)^((p - 1) / 2)

  }
  pmax(pmin(cdf_1 - lambda_k * G_k, 1), 0)

}
