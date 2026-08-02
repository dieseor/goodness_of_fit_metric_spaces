
## Functions

tau_k_d <- function(k, d) {

  if (d == 1) {

    return(2 - (k == 0))

  } else {

    return(1 + 2 * k / (d - 1))

  }

}

vec_kron <- function(x, y) {

  stopifnot(is.matrix(x) && is.matrix(y))
  n <- nrow(x)
  stopifnot(n == nrow(y))
  vec <- matrix(NA, nrow = n, ncol = ncol(x) * ncol(y))
  for (i in 1:n) {

    vec[i, ] <- x[i, ] %x% y[i, ]

  }
  return(vec)

}

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

d_proj_car <- function(x, rho, k, p, log = FALSE) {

  # Checks
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)
  stopifnot(p >= 2)

  C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
  log_dens <- w_p(p = p - 1, log = TRUE) - w_p(p = p, log = TRUE) +
    log1p(rho * drop(Gegen_polyn(theta = acos(x), k = k, p = p)) / C1) +
    ((p - 3) / 2) * log(1 - x^2)
  if (log) {

    return(log_dens)

  } else {

    return(exp(log_dens))

  }

}

p_proj_car <- function(x, rho, k, p) {

  # Checks
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)
  stopifnot(p >= 2)

  C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
  cdf_1 <- drop(p_proj_unif(x = x, p = p))
  lambda_k <- rho * w_p(p = p - 1) / w_p(p = p) / C1
  if (p == 2) {

    G_k <- sin(k * acos(x)) / k

  } else {

    G_k <- (p - 2) / (k * (k + p - 2)) *
      drop(Gegen_polyn(theta = acos(x), k = k - 1, p = p + 2)) *
      (1 - x^2)^((p - 1) / 2)

  }
  pmax(pmin(cdf_1 - lambda_k * G_k, 1), 0)

}

q_proj_car <- function(u, rho, k, p, exact = TRUE) {

  # Invert the cdf
  stopifnot(all(u >= 0 & u <= 1))
  sapply(u, function(uu) {
    uniroot(function(x) p_proj_car(x, rho = rho, k = k, p = p) - uu,
            interval = c(-1, 1))$root
  })

}

r_sph_car <- function(n, mu, rho, k, rejection = TRUE, odd_trick = TRUE,
                      cubic = TRUE) {

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

  } else if (rejection) {

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
                              k = k, rejection = TRUE)
      samp <- rbind(samp, samp_extra)

    }

  } else {

    r_U <- function(n) r_unif_sphere(n = n, p = p - 1)
    r_V <- function(n) q_proj_car(u = runif(n), rho = rho, k = k, p = p)
    samp <- r_tang_norm(n = n, theta = mu, r_U = r_U, r_V = r_V)

  }
  return(samp)

}

d_proj_car_gamma <- function(x, rho, k, p, mu, gamma, log = FALSE) {

  # Checks
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)
  stopifnot(p >= 2)
  stopifnot(length(mu) == length(gamma))

  gamma_mu <- max(min(sum(gamma * mu), 1), -1)
  C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
  log_dens <- w_p(p = p - 1, log = TRUE) - w_p(p = p, log = TRUE) +
    log1p(rho *
            drop(Gegen_polyn(theta = gamma_mu, k = k, p = p)) *
            drop(Gegen_polyn(theta = acos(x), k = k, p = p)) / C1^2
    ) + ((p - 3) / 2) * log(1 - x^2)
  if (log) {

    return(log_dens)

  } else {

    return(exp(log_dens))

  }

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

mom_sph_car <- function(mu, rho, k, m, particular = TRUE) {

  p <- length(mu)
  d <- p - 1
  if (particular && k == m) {

    if (k == 1) {

      moment <- rho / (d + 1) * mu

    } else if (k == 2) {

      vec_Id <- c(diag(rep(1, d + 1)))
      mu <- rbind(mu)
      mu_2 <- vec_kron(mu, mu)
      moment <- vec_Id / (d + 1) +
        ((2 * rho) / (d * (d + 3))) * (vec_kron(mu, mu) - vec_Id / (d + 1))

    } else if (k == 3) {

      vec_Id <- rbind(c(diag(rep(1, d + 1))))
      mu <- rbind(mu)
      mu_2 <- vec_kron(mu, mu)
      mu_3 <- vec_kron(mu_2, mu)
      S3 <- ks:::Sdr(d = d + 1, r = 3)
      moment <- (6 * rho) / (d * (d + 1) * (d + 5)) *
        S3 %*% drop(mu_3 - (3 / (d + 3)) * vec_kron(vec_Id, mu))

    } else if (k == 4) {

      vec_Id <- rbind(c(diag(rep(1, d + 1))))
      vec_Id_2 <- vec_kron(vec_Id, vec_Id)
      mu <- rbind(mu)
      mu_2 <- vec_kron(mu, mu)
      mu_3 <- vec_kron(mu_2, mu)
      mu_4 <- vec_kron(mu_3, mu)
      S4 <- ks:::Sdr(d = d + 1, r = 4)
      moment <- 3 / ((d + 1) * (d + 3)) * S4 %*% (drop(vec_Id_2) +
        (24 * rho) / (d * (d + 1) * (d + 2) * (d + 7)) * drop(
          mu_4 - 6 / (d + 5) * vec_kron(vec_Id, mu_2) +
          3 / ((d + 3) * (d + 5)) * vec_Id_2
        ))

    } else {

      stop("Particular moments only implemented for k = m = 1, 2, 3, 4.")

    }

  } else {

    stop("General moments not implemented yet.")

  }

  return(drop(moment))

}

mom_emp_sph_car <- function(mu, rho, k, m, M = 1e5) {

  # Large Monte Carlo sample
  p <- length(mu)
  X <- r_sph_car(n = M, mu = mu, rho = rho, k = k)

  # Compute the empirical moments
  Xm <- X
  if (m > 1) {

    for (i in seq_len(m - 1)) {

      Xm <- vec_kron(x = Xm, y = X)

    }

  }
  moment <- colMeans(Xm)
  return(moment)

}

mm_sph_car <- function(X, k, mu = NULL, rho_sign = NULL, rho_trunc = TRUE) {

  # Method of moments estimator
  p <- ncol(X)
  d <- p - 1
  eig <- NULL
  if (k == 1 && is.null(mu)) {

    x_bar <- colMeans(X)
    rho_hat <- sqrt(sum(x_bar^2))
    mu_hat <- x_bar / rho_hat
    rho_hat <- (d + 1) * rho_hat

  } else if (k == 2 && is.null(mu)) {

    S <- crossprod(X) / nrow(X)
    eig <- eigen(S, symmetric = TRUE)

    # If rho_sign not provided, see which eigenvalue is more different from the
    # rest, either the first or the last
    if (is.null(rho_sign)) {

      diff_first <- abs(eig$values[1] - mean(eig$values[-1]))
      diff_last <- abs(eig$values[length(eig$values)] -
                        mean(eig$values[-length(eig$values)]))
      rho_sign <- ifelse(diff_first >= diff_last, 1, -1)

    }
    mu_hat <- eig$vectors[, ifelse(rho_sign == 1, 1, length(eig$values))]
    lambda_1 <- eig$values[ifelse(rho_sign == 1, 1, length(eig$values))]
    rho_hat <- (d + 3) / 2 * ((d + 1) * lambda_1 - 1)

  } else if (!is.null(mu)) {

    mu_hat <- mu
    mu_hat <- mu_hat / sqrt(sum(mu_hat^2))
    rho_hat <- tau_k_d(k = k, d = d) *
      mean(Gegen_polyn(theta = acos(X %*% mu_hat), k = k, p = p))

  } else {

    stop("Method of moments only implemented for k = 1, 2 or known mu.")

  }

  # Truncate rho_hat if needed
  if (rho_trunc) {

    rho_hat <- max(min(rho_hat, 1), ifelse(k %% 2 == 0, -1, 0))

  }
  return(list("mu" = mu_hat, "rho" = rho_hat, "k" = k, "eig" = eig))

}

avar_mm_rho <- function(rho, k, p) {

  stopifnot(k == 1 || k == 2)
  d <- p - 1
  switch(k,
         "1" = d + 1 - rho^2,
         "2" = d * (d + 3) / 2 + 2 * (d - 1) * (d + 3) / (d + 5) * rho - rho^2)

}

avar_mm_mu <- function(rho, k, p) {

  stopifnot(k == 1 || k == 2)
  d <- p - 1
  switch(k,
         "1" = (d + 1) / rho^2,
         "2" = (d * (d + 3) * (d * (d + 5) + 2 * rho * (d - 1))) /
                  (4 * rho^2 * (d + 1) * (d + 5)))

}

avar_mle_rho <- function(rho, k, p, exact = TRUE) {

  d <- p - 1
  if (d == 1) {

    if (exact) {

      A_k <- (1 - sqrt(1 - rho^2)) / (rho^2 * sqrt(1 - rho^2))

    } else {

      A_k <- (1 / pi) * sapply(rho, function(r) {
        integrate(function(th) cos(k * th)^2 / (1 + r * cos(k * th)),
                  lower = 0, upper = pi)$value
      })

    }

  } else {

    if (exact && d == 2 && k == 1) {

      A_k <- (atanh(rho) - rho) / rho^3

    } else if (exact && d == 2 && k == 2) {

      A_k <- (2 * atan(sqrt(3 * rho / (2 - rho))) /
         sqrt(3 * rho * (2 - rho)) - 1) / rho^2

    } else {

      C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
      A_k <- rotasym::w_p(p - 1) / rotasym::w_p(p) / C1 *
        sapply(rho, function(r) integrate(function(t)
          Gegen_polyn(theta = acos(t), k = k, p = p)^2 *
          (1 - t^2)^(d / 2 - 1) /
          (C1 + r * Gegen_polyn(theta = acos(t), k = k, p = p)),
          lower = -1, upper = 1)$value)

    }

  }
  return(1 / A_k)

}

avar_mle_mu <- function(rho, k, p, exact = TRUE) {

  d <- p - 1
  if (d == 1) {

    if (exact) {

      B_k <- k^2 * (1 - sqrt(1 - rho^2)) / rho^2

    } else {

      B_k <- (k^2 / pi) * sapply(rho, function(r) {
        integrate(function(th) sin(k * th)^2 / (1 + r * cos(k * th)),
                  lower = 0, upper = pi)$value
      })

    }

  } else {

    if (exact && d == 2 && k == 1) {

      B_k <- (rho - (1 - rho^2) * atanh(rho)) / (2 * rho^3)

    } else if (exact && d == 2 && k == 2) {

      B_k <- (2 + rho - 2 * (2 + rho - rho^2) *
         atan(sqrt(3 * rho / (2 - rho))) / sqrt(3 * rho * (2 - rho)))
      B_k <- B_k / rho^2

    } else {

      C1 <- drop(Gegen_polyn(theta = 0, k = k, p = p))
      B_k <- rotasym::w_p(p - 1) / rotasym::w_p(p) / C1 * (d - 1)^2 / d *
        sapply(rho, function(r) integrate(function(t)
          Gegen_polyn(theta = acos(t), k = k - 1, p = p + 2)^2 *
          (1 - t^2)^(d / 2) /
          (C1 + r * Gegen_polyn(theta = acos(t), k = k, p = p)),
          lower = -1, upper = 1)$value)

    }

  }
  return(1 / (rho^2 * B_k))

}

avar_gm_rho <- function(rho, k, p) {

  d <- p - 1
  if (d == 1) {

    eta <- 0

  } else {

    eta <- ifelse(k %% 2 == 0, (2 * k + d - 1)^2 / ((3 * k + d - 1) * (d - 1)) *
                    factorial(k) / (factorial(k / 2)^3) *
                    gamma((d + k - 1) / 2)^3 * gamma(d + 3 * k / 2 - 1) /
                    (gamma(d + k - 1) * gamma((d - 1) / 2)^2 *
                      gamma((d + 3 * k - 1) / 2)), 0)

  }
  return(sphunif::d_p_k(p = p, k = k) + rho * eta - rho^2)

}

mle_sph_car <- function(X, k, mu0 = NULL, rho0 = NULL, ...) {

  ll <- function(par) {
    -sum(d_sph_car(x = X, mu = par[-1] / sqrt(sum(par[-1]^2)),
                   # rho = max(min(par[1], 1), ifelse(k %% 2 == 0, -1, 0)),
                   rho = max(min(par[1], 1), 0),
                   k = k, log = TRUE))
  }
  if (is.null(mu0)) {
    mu0 <- colMeans(X)
    mu0 <- mu0 / sqrt(sum(mu0^2))
  }
  if (is.null(rho0)) {
    rho0 <- sqrt(sum(mu0^2))
    rho0 <- max(min(rho0, 1), 0)
  }
  opt <- optim(par = c(rho0, mu0), fn = ll, method = "L-BFGS-B",
               # lower = c(ifelse(k %% 2 == 0, -1, 0), rep(-Inf, length(mu0))),
               # upper = c(1, rep(Inf, length(mu0))),
               lower = c(0, rep(-Inf, length(mu0))),
               upper = c(1, rep(Inf, length(mu0))),
               ...)
  opt$par[-1] <- opt$par[-1] / sqrt(sum(opt$par[-1]^2))
  # opt$par[1] <- max(min(opt$par[1], 1), ifelse(k %% 2 == 0, -1, 0))
  opt$par[1] <- max(min(opt$par[1], 1), 0)
  return(list("mu" = opt$par[-1], "rho" = opt$par[1], "k" = k,
              "ll" = -opt$value, "opt" = opt))

}

cvm_stat <- function(u) {

  # Computational form
  u <- sort(u)
  n <- length(u)
  i <- seq_len(n)
  stat <- sum((u - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  return(stat)

}

ad_stat <- function(u) {

  # Trim degenerate cases before computational form
  u <- sort(u)
  n <- length(u)
  if (u[n] == 1) {
    u <- u[-n]
    n <- length(u)
  }
  if (u[1] == 0) {
    u <- u[-1]
    n <- length(u)
  }

  # Computational form
  i <- seq_len(n)
  logs <- (2 * i - 1) * log(u) + (2 * (n - i) + 1) * log(1 - u)
  stat <- - n - mean(logs)
  return(stat)

}

sph_car_stat_gof <- function(X, theta, weight = c("CvM", "AD")[1],
                             type_stat = c("Pn", "unif_mc", "unif_w",
                                           "unif_cvm", "C_k")[1],
                             gammas = NULL) {

  # Match weight and type_stat
  weight <- match.arg(weight, choices = c("CvM", "AD"))
  type_stat <- match.arg(type_stat, choices = c("Pn", "unif_mc", "unif_w",
                                                "unif_cvm", "C_k"))

  # Dimension checks
  stopifnot(is.matrix(X))
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(theta) == p + 2)

  # Parameters
  k <- theta[1]
  rho <- theta[2]
  mu <- theta[-c(1, 2)]

  # Check gammas or generate them for C_k(\hat\mu, \hat\rho) weight
  if (type_stat %in% c("unif_mc", "unif_w") && is.null(gammas)) {

    stop("gammas are needed for 'type_stat = \"unif_mc\", \"unif_w\"'.")

  } else if (type_stat == "C_k") {

    gammas <- r_sph_car(n = n, mu = mu, rho = rho, k = k,
                        rejection = TRUE, odd_trick = TRUE, cubic = TRUE)

  }

  # Statistic computation
  if (type_stat %in% c("unif_mc", "C_k")) {

    # Loop on gammas
    stats <- numeric(nrow(gammas))
    for (i in seq_len(nrow(gammas))) {

      X_gamma <- X %*% gammas[i, ]
      X_gamma <- pmax(pmin(X_gamma, 1), -1)
      P <- p_proj_car_gamma(x = X_gamma, mu = mu, rho = rho, k = k,
                            p = p, gamma = gammas[i, ])
      stats[i] <- switch(weight,
                         "CvM" = cvm_stat(P),
                         "AD" = ad_stat(P))

    }
    return(mean(stats))

  } else if (type_stat == "unif_w") {

    # CvM weights
    if (weight != "CvM") {

      stop("Only 'CvM' weight implemented for 'unif_w' statistic.")

    }
    W <- function(x) x
    W_1 <- function(x) x^2 / 2
    W_2 <- function(x) x^3 / 3
    bias <- W(1) + W_2(1) - W_2(0) - 2 * W_1(1)

    # First term (single sum)
    term_1 <- numeric(n)
    for (i in seq_len(n)) {

      gamma_Xi <- gammas %*% X[i, ]
      gamma_Xi <- pmax(pmin(gamma_Xi, 1), -1)
      Ws <- sapply(seq_len(nrow(gammas)), function(l) {
        F_gamma_i <- p_proj_car_gamma(x = gamma_Xi[l, ], mu = mu, rho = rho,
                                      k = k, p = p, gamma = gammas[l, ])
        W_1(F_gamma_i)
      })
      term_1[i] <- mean(Ws)

    }

    # Second term (double sum)
    term_2 <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      for (j in 1:i) {

        gamma_Xi <- gammas %*% X[i, ]
        gamma_Xj <- gammas %*% X[j, ]
        gamma_Xi <- pmax(pmin(gamma_Xi, 1), -1)
        gamma_Xj <- pmax(pmin(gamma_Xj, 1), -1)
        Ws <- sapply(seq_len(nrow(gammas)), function(l) {
          F_gamma_max <- p_proj_car_gamma(x = max(gamma_Xi[l], gamma_Xj[l]),
                                          mu = mu, rho = rho, k = k, p = p,
                                          gamma = gammas[l, ])
          W(F_gamma_max)
        })
        term_2[i, j] <- mean(Ws)

      }

    }
    term_2 <- term_2 + t(term_2)
    diag(term_2) <- diag(term_2) / 2

    # Statistic
    stat <- bias + 2 * mean(term_1) - mean(term_2)
    return(n * stat)

  } else if (type_stat == "unif_cvm") {

    # Kernels
    kernel_phi <- function(x_mu) {

      if (p == 2) {

        phi <- rho / (2 * pi^2 * k^2) *
          (Gegen_polyn(theta = acos(x_mu), k = k, p = 2) -
            (rho / 4) * (2 - Gegen_polyn(theta = acos(x_mu), k = 2 * k, p = 2)))
        phi <- drop(phi)

      } else if (p == 3) {

        if (k == 1) {

          phi <- rho / 30 * x_mu - rho^2 / 4 * (2 / 35 - 4 * x_mu^2 / 105)

        } else if (k == 2) {

          phi <- rho / 420 * (3 * x_mu^2 - 1) -
            rho^2 / 4 * (1 / 330 + 3 * x_mu^2 / 385 - x_mu^4 / 110)

        } else {

          stop("Only k = 1 and k = 2 implemented for p = 3 in 'unif_cvm' statistic.")

        }

      } else {

        stop("Only p = 2 and p = 3 implemented for 'unif_cvm' statistic.")

      }
      return(phi)

    }
    kernel_psi <- function(x_x, x_mu) {

      x1_mu <- matrix(c(x_mu), nrow = nrow(x_x), ncol = ncol(x_x), byrow = TRUE)
      x2_mu <- matrix(c(x_mu), nrow = nrow(x_x), ncol = ncol(x_x), byrow = FALSE)
      psi_cvm <- sphunif::psi_Pn(theta = acos(x_x), q = p - 1, type = "PCvM")
      if (p == 2) {

        t_12 <- (x1_mu + x2_mu) / sqrt(2 * (1 + x_x))
        t_12 <- pmax(pmin(t_12, 1), -1)
        psi <- rho * (
          (pi - acos(x_x)) / (2 * pi^2 * k) *
            drop(Gegen_polyn(theta = acos(t_12), k = k, p = 2)) *
            sin(k * acos(x_x) / 2)

        )

      } else if (p == 3) {

        if (k == 1) {

          psi <- rho / 32 * sqrt((1 - x_x) / 2) * (x1_mu + x2_mu)

        } else if (k == 2) {

          psi <- rho / 128 *
            sqrt((1 - x_x) / 2) * (
              (1 + x_x) / 2 +
              (3 * (3 * x_x - 1)) / (4 * (1 - x_x)) *
              (x1_mu^2 + x2_mu^2) +
              (3 * (x_x - 3)) / (2 * (1 - x_x)) *
              (x1_mu * x2_mu)
            )

        } else {

          stop("Only k = 1, 2 implemented for p = 3 in 'unif_cvm' statistic.")

        }

      } else {

        stop("Only p = 2, 3 implemented for 'unif_cvm' statistic.")

      }
      psi <- psi_cvm - psi
      return(psi)

    }

    # First term (single sum)
    X_mu <- drop(pmax(pmin(X %*% mu, 1), -1))
    phis <- kernel_phi(x_mu = X_mu)

    # Second term (double sum)
    X_X <- pmax(pmin(tcrossprod(X), 1), -1)
    psis <- kernel_psi(x_x = X_X, x_mu = X_mu)
    psis <- psis[lower.tri(psis)]

    # Statistic
    stat <- (3 - 2 * n) / (6 * n) - mean(phis) + 2 * sum(psis) / n^2
    return(n * stat)

  } else if (type_stat == "Pn") {

    # Loop on sample
    stats <- numeric(n)
    for (i in seq_len(n)) {

      X_gamma <- X %*% X[i, ]
      X_gamma <- pmax(pmin(X_gamma, 1), -1)
      P <- p_proj_car_gamma(x = X_gamma, mu = mu, rho = rho, k = k,
                            p = p, gamma = X[i, ])
      stats[i] <- switch(weight,
                         "CvM" = cvm_stat(P),
                         "AD" = ad_stat(P))

    }
    return(mean(stats))

  }

}

sph_car_gof <- function(x, k, est_type = c("mle", "mm", "simple")[1],
                        mu0 = NULL, rho0 = NULL, rho_sign = NULL,
                        weight = c("CvM", "AD")[1],
                        type_stat = c("Pn", "unif_mc", "unif_w",
                                      "unif_cvm", "C_k")[1],
                        K = 100, B = 5e3, plot_boot = TRUE,
                        parallel = TRUE) {

  # Check weight and type_stat
  n_stats <- length(type_stat)
  stopifnot(weight %in% c("CvM", "AD"))
  stopifnot(type_stat %in% c("Pn", "unif_mc", "unif_w", "unif_cvm", "C_k"))
  stopifnot(n_stats == length(weight))

  # Common random projections
  gammas <- rotasym::r_unif_sphere(n = K, p = ncol(x))

  # Type of theta estimate
  if (est_type == "simple") {

    # No estimation, use mu0 and rho0
    if (is.null(mu0) || is.null(rho0)) {

      stop("For 'simple' estimation, mu0 and rho0 must be provided.")

    }
    theta_hat <- c(k, rho0, mu0)

  } else if (est_type == "mle") {

    # If given, use mu0 and rho0 as initial values
    est <- mle_sph_car(X = x, k = k, mu0 = mu0, rho0 = rho0)
    theta_hat <- c(est$k, est$rho, est$mu)

  } else if (est_type == "mm") {

    # If given, use mu0 as known value and use GM estimator for rho
    est <- mm_sph_car(X = x, k = k, mu = mu0, rho_sign = rho_sign,
                      rho_trunc = TRUE)
    theta_hat <- c(est$k, est$rho, est$mu)

  } else {

    stop("est_type must be either 'simple', 'mle' or 'mm'.")

  }

  # Compute all test statistics
  Tn <- numeric(n_stats)
  for (i in seq_len(n_stats)) {

    Tn[i] <- sph_car_stat_gof(X = x, theta = theta_hat,
                              weight = weight[i], type_stat = type_stat[i],
                              gammas = gammas)

  }
  names(Tn) <- paste0(type_stat, "_", weight)

  # Parallel stuff
  prog <- progressr::progressor(along = seq_len(B))
  if (parallel) {

    `%op%` <- doRNG::`%dorng%`

  } else {

    `%op%` <- foreach::`%do%`

  }

  # Bootstrap resampling
  n <- nrow(x)
  boot <- foreach::foreach(b = 1:B, .combine = rbind) %op% {

    # Progress
    if (requireNamespace("progressr", quietly = TRUE)) prog()

    # Sample cardioid
    x_star <- r_sph_car(n = n, mu = theta_hat[-c(1, 2)],
                        rho = theta_hat[2], k = k)

    # Type of bootstrap theta estimate
    if (est_type == "simple") {

      # No estimation, use mu0 and rho0
      theta_hat_star <- c(k, rho0, mu0)

    } else if (est_type == "mle") {

      # Use mu_hat and rho_hat as initial values
      est_star <- mle_sph_car(X = x_star, k = k, mu0 = theta_hat[-c(1, 2)],
                              rho0 = theta_hat[2])
      theta_hat_star <- c(est_star$k, est_star$rho, est_star$mu)

    } else if (est_type == "mm") {

      # If given, use mu0 as known value and use GM estimator for rho
      est_star <- mm_sph_car(X = x_star, k = k, mu = mu0,
                             rho_sign = rho_sign, rho_trunc = TRUE)
      theta_hat_star <- c(est_star$k, est_star$rho, est_star$mu)

    }

    # Test statistics
    stat_star <- numeric(n_stats)
    for (i in seq_len(n_stats)) {

      stat_star[i] <- tryCatch(
        sph_car_stat_gof(X = x_star, theta = theta_hat_star,
                         weight = weight[i], type_stat = type_stat[i],
                         gammas = gammas),
        error = function(e) NA)

    }

    # Return test statistics and bootstrap estimators
    c(stat_star, theta_hat_star)

  }
  boot <- rbind(boot)
  colnames(boot) <- c(names(Tn), "k", "rho_hat",
                      paste0("mu_hat_", seq_len(ncol(x))))

  # Extract bootstrap statistics and estimators
  Tn_star <- boot[, seq_len(n_stats), drop = FALSE]
  theta_hat_star <- boot[, -seq_len(n_stats), drop = FALSE]

  # Test information
  method <- paste0("Bootstrap-based test for cardioid, k = ", k,
                   " (", paste(paste(weight, type_stat, sep = "-"),
                               collapse = "/"), " variant)")
  alternative <- "any alternative to spherical cardioid"

  # p-values
  pvalues <- (rowSums(t(rbind(Tn_star)) > Tn, na.rm = TRUE) + 1) / (B + 1)

  # Construct an "htest" result
  result <- list(statistic = Tn, pvalue = pvalues,
                 theta_hat = theta_hat, theta_hat_star = unname(theta_hat_star),
                 statistic_boot = Tn_star, B = B,
                 alternative = alternative, method = method,
                 data.name = deparse(substitute(x)))
  class(result) <- "htest"

  # Plot the position of the original statistic with respect to the
  # bootstrap replicates?
  if (plot_boot) {

    old_par <- par(mfrow = c(min(2, n_stats), ceiling(n_stats / 2)))
    for (i in seq_len(n_stats)) {

      hist(result$statistic_boot[, i], probability = TRUE,
           main = paste0(names(result$statistic)[i], " p-value: ",
                         round(result$pvalue[i], 4)),
           xlab = latex2exp::TeX(paste0("$P_n^{", weight[i], ", ",
                                        type_stat[i], ", *}$")))
      rug(result$statistic_boot[, i])
      abline(v = result$statistic[i], col = 2, lwd = 2)
      text(x = result$statistic[i], y = 1.5 * mean(par()$usr[3:4]),
           labels = latex2exp::TeX(paste0("$P_n^{", weight[i], ", ",
                                          type_stat[i], "}$")),
           col = 2, pos = 4)

    }
    par(old_par)

  }

  # Return "htest"
  return(result)

}

sph_car_mgf <- function(t, mu, rho, k) {

  # Checks
  p <- length(mu)
  t <- rbind(t)
  stopifnot(ncol(t) == p)
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)

  # Dimension
  d <- p - 1

  # Norm and normalization of t
  norm_t <- sqrt(rowSums(t^2))
  t <- t / norm_t

  # Bessel functions and e_{k, d}
  bessel_0 <- besselI(x = norm_t, nu = (d - 1) / 2)
  bessel_k <- besselI(x = norm_t, nu = (2 * k + d - 1) / 2)
  e_k_d <- function(k , d) {

    if (d == 1) {

      return(2 - (k == 0))

    } else {

      return((gamma((d - 1) / 2) / 2) * (2 * k + d - 1))

    }

  }

  # Mgf
  mgf <- (2 / norm_t)^((d - 1) / 2) * (
    e_k_d(k = 0, d = d) * bessel_0 +
    rho * e_k_d(k = k, d = d) / sphunif::d_p_k(p = p, k = k) *
    bessel_k * drop(Gegen_polyn(theta = acos(t %*% mu), k = k, p = p))
  )
  return(unname(mgf))

}

sph_car_cf <- function(t, mu, rho, k) {

  # Checks
  p <- length(mu)
  t <- rbind(t)
  stopifnot(ncol(t) == p)
  stopifnot(-1 <= rho && rho <= 1)
  stopifnot(k >= 1)

  # Dimension
  d <- p - 1

  # Norm and normalization of t
  norm_t <- sqrt(rowSums(t^2))
  t <- t / norm_t

  # Bessel functions and e_{k, d}
  bessel_0 <- besselJ(x = norm_t, nu = (d - 1) / 2)
  bessel_k <- besselJ(x = norm_t, nu = (2 * k + d - 1) / 2)
  e_k_d <- function(k , d) {

    if (d == 1) {

      return(2 - (k == 0))

    } else {

      return((gamma((d - 1) / 2) / 2) * (2 * k + d - 1))

    }

  }

  # Cf
  cf <- (2 / norm_t)^((d - 1) / 2) * (
    e_k_d(k = 0, d = d) * bessel_0 +
    rho * (1i)^k * e_k_d(k = k, d = d) / sphunif::d_p_k(p = p, k = k) *
    bessel_k * drop(Gegen_polyn(theta = acos(t %*% mu), k = k, p = p))
  )
  return(unname(cf))

}
