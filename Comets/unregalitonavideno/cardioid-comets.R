
library(sphunif)
library(rotasym)
library(rgl)
library(future)
library(doFuture)
library(foreach)
library(progressr)
library(dplyr)
library(viridis)
library(ggplot2)
library(goftest)

# Load functions
source("cardioid-source.R")

# Progress bar
handlers(handler_progress(
  format = paste("(:spin) [:bar] :percent Iter: :current/:total Rate:",
                 ":tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull"),
  clear = FALSE))

## Load data and preprocess it
{

# Load data
stopifnot(packageVersion("sphunif") >= "1.3.0")
data("comets", package = "sphunif")

# Add normal
comets$normal <- cbind(sin(comets$i) * sin(comets$om),
                       -sin(comets$i) * cos(comets$om),
                       cos(comets$i))

# Oort cloud (long-period comets)
comets_oort <- subset(x = comets,
                      subset = !(class %in% c("HYP", "PAR")) & per_y >= 200)
nrow(comets_oort)

# Comets without fragments
comets_oort <- comets_oort[!comets_oort$frag, ]
nrow(comets_oort)

# Kuiper belt (short-period comets)
comets_kuiper <- subset(x = comets, subset = !(class %in% c("HYP", "PAR")) &
                          per_y < 200)
nrow(comets_kuiper)

# Cleaned data
comets_kuiper <- comets_kuiper[!comets_kuiper$frag, ]
nrow(comets_kuiper)

}

## Goodness-of-fit tests
{

# Setup
ws <- rep(c("CvM", "AD"), 3)
ty2 <- c("unif_cvm", "unif_mc", rep("Pn", 2), rep("C_k", 2))
ty3 <- c("unif_mc", "unif_mc", rep("Pn", 2), rep("C_k", 2))
doRNG::registerDoRNG()
# doFuture::registerDoFuture()
future::plan(future::multisession(), workers = 8)
K <- 1e3
B <- 1e3

# Fit with MLE using mu0 = c(0, 0, 1) as initial value
est <- "mle"
mu0 <- c(0, 0, 1)

# Tests for Oort
set.seed(42)
with_progress({
  gof_oor_0 <- sph_car_gof(comets_oort$normal, k = 1, K = K, B = B,
                           est_type = "simple", mu0 = c(0, 0, 1), rho0 = 0,
                           type_stat = ty2, weight = ws)
  gof_oor_1 <- sph_car_gof(comets_oort$normal, k = 1, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty2, weight = ws)
  gof_oor_2 <- sph_car_gof(comets_oort$normal, k = 2, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty2, weight = ws)
  gof_oor_3 <- sph_car_gof(comets_oort$normal, k = 3, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty3, weight = ws)
  gof_oor_4 <- sph_car_gof(comets_oort$normal, k = 4, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty3, weight = ws)
})

# Tests for Kuiper
set.seed(42)
with_progress({
  gof_kui_0 <- sph_car_gof(comets_kuiper$normal, k = 1, K = K, B = B,
                           est_type = "simple", mu0 = c(0, 0, 1), rho0 = 0,
                           type_stat = ty2, weight = ws)
  gof_kui_1 <- sph_car_gof(comets_kuiper$normal, k = 1, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty2, weight = ws)
  gof_kui_2 <- sph_car_gof(comets_kuiper$normal, k = 2, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty2, weight = ws)
  gof_kui_3 <- sph_car_gof(comets_kuiper$normal, k = 3, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty3, weight = ws)
  gof_kui_4 <- sph_car_gof(comets_kuiper$normal, k = 4, K = K, B = B,
                           est_type = est, mu0 = mu0,
                           type_stat = ty3, weight = ws)
})

# p-value tables
gof_oor <- rbind(
  gof_oor_0$pvalue,
  gof_oor_1$pvalue,
  gof_oor_2$pvalue,
  gof_oor_3$pvalue,
  gof_oor_4$pvalue
)
gof_kui <- rbind(
  gof_kui_0$pvalue,
  gof_kui_1$pvalue,
  gof_kui_2$pvalue,
  gof_kui_3$pvalue,
  gof_kui_4$pvalue
)
gof_oor |>
  knitr::kable(format = "latex", digits = 3, booktabs = TRUE)
gof_kui |>
  knitr::kable(format = "latex", digits = 3, booktabs = TRUE)

# Save results
save(list = c("gof_oor_0", "gof_oor_1", "gof_oor_2", "gof_oor_3", "gof_oor_4",
              "gof_kui_0", "gof_kui_1", "gof_kui_2", "gof_kui_3", "gof_kui_4"),
     file = "comets_gof_BK_1e3.Rdata")

}

## Analyze test results and bootstrap CIs
{

# Load tests
load("comets_gof_BK_1e3.Rdata")

# Estimates and bootstrap CIs for Oort
mu_hat_oort <- gof_oor_2$theta_hat[-c(1, 2)]
rho_hat_oort <- gof_oor_2$theta_hat[2]
t_mu_boot_oort <- drop(gof_oor_2$theta_hat_star[, -c(1, 2)] %*%
                         gof_oor_2$theta_hat[-c(1, 2)])
r_mu_boot_oort <- quantile(abs(t_mu_boot_oort), probs = 0.05)
round(mu_hat_oort, 4)
round(rho_hat_oort, 4)
round(r_mu_boot_oort, 4)
round(quantile(gof_oor_2$theta_hat_star[, 2], probs = c(0.025, 0.975)), 4)

# Estimates and bootstrap CIs for Kuiper
mu_hat_kuiper <- gof_kui_1$theta_hat[-c(1, 2)]
rho_hat_kuiper <- gof_kui_1$theta_hat[2]
t_mu_boot_kuiper <- drop(gof_kui_1$theta_hat_star[, -c(1, 2)] %*%
                           gof_kui_1$theta_hat[-c(1, 2)])
r_mu_boot_kuiper <- quantile(t_mu_boot_kuiper, probs = 0.05)
round(mu_hat_kuiper, 4)
round(rho_hat_kuiper, 4)
round(r_mu_boot_kuiper, 4)
round(quantile(gof_kui_1$theta_hat_star[, 2], probs = c(0.025, 0.975)))

# "Best" cardioid fit
rbind(gof_oor_0$statistic, gof_oor_1$statistic, gof_oor_2$statistic,
      gof_oor_3$statistic, gof_oor_4$statistic) |>
  apply(2, which.min)
rbind(gof_kui_0$statistic, gof_kui_1$statistic, gof_kui_2$statistic,
      gof_kui_3$statistic, gof_kui_4$statistic) |>
  apply(2, which.min)

}

## Plotting functions
{

# Plotting function for sphere
add_sphere <- function(x = 0, y = 0, z = 0, radius = 1, n = 101,
                       plane = FALSE, axes = FALSE, ...) {

  # Sphere surface
  sph <- function(s, t) {
    cbind(radius * cos(t) * cos(s) + x,
          radius * sin(s) + y,
          radius * sin(t) * cos(s) + z)
  }
  rgl::persp3d(sph, slim = c(-pi, pi) / 2, tlim = c(0, 2 * pi), n = n,
               add = TRUE, ...)

  # Add equator plane?
  N <- 201
  if (plane) {
    disk <- function(s, t) {
      cbind(t * radius * cos(s) + x,
            t * radius * sin(s) + y,
            z)
    }
    rgl::persp3d(disk, slim = c(-pi, pi), tlim = c(0, 1), n = N,
                 add = TRUE, lit = FALSE, alpha = 0.1, col = "black")
  }

  # Equator curve
  th <- seq(0, 2 * pi, l = N)
  rgl::lines3d(radius * cbind(cos(th), sin(th), 0), lwd = 2)

  # Add dashed axes inside the sphere?
  if (axes) {

    add_dashed <- function(P, lwd = 2, col = "black") {

      idx <- seq(1, nrow(P) - 1, by = 2)
      sidx <- as.vector(rbind(idx, idx + 1))  # (1,2),(3,4),...
      rgl::segments3d(P[sidx, , drop = FALSE], lwd = lwd, col = col)

    }
    u <- seq(-1, 1, length.out = N) * radius
    # ax_x <- cbind(u, 0, 0) + c(x, y, z)
    # add_dashed(ax_x)
    # ax_y <- cbind(0, u, 0) + c(x, y, z)
    # add_dashed(ax_y)
    ax_z <- cbind(0, 0, u) + c(x, y, z)
    add_dashed(ax_z)

    # Add equator-like meridians?
    # th <- seq(0, 2 * pi, l = N)
    # rgl::lines3d(radius * cbind(0, cos(th), sin(th)), lwd = 2)
    # rgl::lines3d(radius * cbind(cos(th), 0, sin(th)), lwd = 2)

  }

}

# Plotting function for spherical cap at mu with angle theta
add_sphere_cap <- function(x = 0, y = 0, z = 0, radius = 1, n = 101,
                           mu = c(0, 0, 1), theta = pi / 6, border = TRUE,
                           ...) {

  # Orthonormal frame (u, v, mu)
  mu <- mu / sqrt(sum(mu^2))
  a <- if (abs(mu[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)

  u <- a - sum(a * mu) * mu
  u <- u / sqrt(sum(u^2))

  v <- c(mu[2] * u[3] - mu[3] * u[2],
         mu[3] * u[1] - mu[1] * u[3],
         mu[1] * u[2] - mu[2] * u[1])

  # Parametrization of the cap (vectorized!)
  cap <- function(alpha, phi) {

    ca <- cos(alpha)
    sa <- sin(alpha)
    cp <- cos(phi)
    sp <- sin(phi)

    # n x 3 matrices
    Mu <- ca %o% mu
    U  <- (sa * cp) %o% u
    V  <- (sa * sp) %o% v

    P <- Mu + U + V

    cbind(radius * P[, 1] + x,
          radius * P[, 2] + y,
          radius * P[, 3] + z)

  }

  rgl::persp3d(cap,
               slim = c(0, theta),
               tlim = c(0, 2 * pi),
               n = n,
               add = TRUE, ...)

  if (border) {
    N <- 201
    phi <- seq(0, 2 * pi, length.out = N)
    rim <- cap(rep(theta, N), phi)
    rgl::lines3d(rim, lwd = 2)
  }

  invisible(NULL)

}

# Plotting function for ellipses
add_comet_ellipse <- function(a, e, w, normal, L = 1001, log_axis =TRUE,
                              lwd = 3, col = "blue", polygon = TRUE, ...) {

  # The parametric equation of an ellipse with semi-major axis a and semi-minor
  # axis b is (a * cos(t), b * sin(t)), which has foci at (-c, 0) and (+c, 0)
  # for c = sqrt(a^2 - b^2). Therefore, the ellipse with foci (0, 0) and (2c, 0)
  # is parametrized by (c + a * cos(t), b * sin(t)). The eccentricity e is
  # defined as e = c / a = sqrt(1 - b^2 / a^2). Therefore, given a and e,
  # b = a * sqrt(1 - e^2).

  # Log10-transform a?
  if (log_axis) a <- log10(a)

  # Compute semi-minor axis
  b <- a * sqrt(1 - e^2)

  # Compute c
  c <- sqrt(a^2 - b^2)

  # Produce ellipse with one foci at (0, 0) and (0, 2 * c) in the plane
  t <- 2 * acos(seq(1, -1, l = L))
  ellip <- cbind(c + a * cos(t), b * sin(t))

  # Rotate in the plane according to the argument of perihelion
  R <- rbind(c(cos(w), -sin(w)),
             c(sin(w), cos(w)))
  ellip <- t(R %*% t(ellip))

  # Move ellipse to the YZ plane with normal vector (1, 0, 0)'
  ellip <- cbind(0, ellip)

  # Rotate the plane according to the desired normal vector
  R <- cbind(normal, rotasym::Gamma_theta(theta = normal, eig = FALSE))
  ellip <- t(R %*% t(ellip))

  # Plot the ellipse
  if (polygon) {
    rgl::polygon3d(ellip[seq(1, L, by = 10), ], lit = FALSE, alpha = 0.2,
                   col = col, random = FALSE)
  }
  rgl::lines3d(ellip, lwd = lwd, col = col)

}

}

## Plots orbits
{

# Sort by eccentricity for better visualization
comets_oort <- comets_oort[order(comets_oort$e), ]
soft_ecc_oort <- rev(which(comets_oort$e < 0.98))
soft_ecc_oort <- soft_ecc_oort[c(4, 17, 18, 22, 24, 32, 48, 40, 60, 85)]
col <- viridis(length(soft_ecc_oort))

# Common view
radius <- 4.75
# P <- matrix(c(-0.8320085,  -0.5583869,  0.03889489, 0,
#               0.1786495,  -0.2705550, 0.96728635, 0,
#               -0.5462106, 0.7483342, 0.32069334, 0,
#               0.0000000,  0.0000000,  0.00000000, 1),
#             nrow = 4, ncol = 4, byrow = TRUE)

# Plot Oort orbits and normals
close3d()
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = radius, col = "lightblue",
           lit = FALSE, alpha = 0.15, plane = TRUE)
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = 0.15, col = "orange", lit = FALSE)
par3d(windowRect = c(50, 50, 1050, 1050), zoom = 0.6)
for (i in seq_along(soft_ecc_oort)) {
  j <- soft_ecc_oort[i]
  add_comet_ellipse(a = comets_oort$a[j], e = comets_oort$e[j],
                    w = comets_oort$w[j],
                    normal = comets_oort$normal[j, ],
                    lwd = 2, col = col[i])
  arrow3d(p0 = c(0, 0, 0), p1 = radius * comets_oort$normal[j, ],
          type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
          lwd = 3, col = col[i])
}
snapshot3d("figs/fig_1a_orbits.png", webshot = FALSE,
           width = 1000, height = 1000, top = FALSE)
close3d()

# # Plot Oort normals to the orbits
# close3d()
# add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = radius, col = "lightblue",
#            lit = FALSE, alpha = 0.15, plane = TRUE)
# par3d(windowRect = c(50, 50, 1050, 1050), zoom = 0.6)
# for (i in seq_along(soft_ecc_oort)) {
#   j <- soft_ecc_oort[i]
#   arrow3d(p0 = c(0, 0, 0), p1 = radius * comets_oort_clean$normal[j, ],
#           type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
#           lwd = 3, col = col[i])
# }
# add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = 0.15, col = "orange", lit = FALSE)
# snapshot3d("figs/Oort-normal.png", webshot = FALSE,
#            width = 1000, height = 1000, top = FALSE)
# close3d()

# Sort by eccentricity for better visualization
comets_kuiper <- comets_kuiper[order(comets_kuiper$e), ]
soft_ecc_kuiper <- seq(1, nrow(comets_kuiper), l = 10)
col <- viridis(length(soft_ecc_kuiper))

# Plot Kuiper orbits and normals
close3d()
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = radius, col = "lightblue",
           lit = FALSE, alpha = 0.15, plane = TRUE)
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = 0.15, col = "orange", lit = FALSE)
par3d(windowRect = c(50, 50, 1050, 1050), zoom = 0.6)
for (i in seq_along(soft_ecc_kuiper)) {
  j <- soft_ecc_kuiper[i]
  add_comet_ellipse(a = comets_kuiper$a[j], e = comets_kuiper$e[j],
                    w = comets_kuiper$w[j],
                    normal = comets_kuiper$normal[j, ],
                    lwd = 2, col = col[i])
  arrow3d(p0 = c(0, 0, 0), p1 = radius * comets_kuiper$normal[j, ],
          type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
          lwd = 3, col = col[i])
}
snapshot3d("figs/fig_1b_orbits.png", webshot = FALSE,
           width = 1000, height = 1000, top = FALSE)
close3d()

}

## Plots points
{

# Plot Oort points
close3d()
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = radius, col = "lightblue",
           lit = FALSE, alpha = 0.15, plane = FALSE, axes = TRUE)
par3d(windowRect = c(50, 50, 1050, 1050), zoom = 0.6)
for (i in seq_len(nrow(comets_oort$normal))) {
  add_sphere(x = radius * comets_oort$normal[i, 1],
             y = radius * comets_oort$normal[i, 2],
             z = radius * comets_oort$normal[i, 3], radius = 0.03,
             col = "black", lit = FALSE)
}
arrow3d(p0 = c(0, 0, 0), p1 = radius * mu_hat_oort,
        type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
        lwd = 5, col = 2)
arrow3d(p0 = c(0, 0, 0), p1 = -radius * mu_hat_oort,
        type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
        lwd = 5, col = 2)
add_sphere_cap(x = 0, y = 0, z = 0, radius = 1.01 * radius,
               mu = mu_hat_oort, theta = acos(r_mu_boot_oort),
               col = "red", lit = FALSE, alpha = 0.3, border = TRUE)
add_sphere_cap(x = 0, y = 0, z = 0, radius = 0.99 * radius,
               mu = -mu_hat_oort, theta = acos(r_mu_boot_oort),
               col = "red", lit = FALSE, alpha = 0.3, border = TRUE)
snapshot3d("figs/fig_1a_points.png", webshot = FALSE,
           width = 1000, height = 1000, top = FALSE)
close3d()

# Plot Kuiper points
close3d()
add_sphere(x0 = 0, y0 = 0, z0 = 0, radius = radius, col = "lightblue",
           lit = FALSE, alpha = 0.15, plane = FALSE, axes = TRUE)
par3d(windowRect = c(50, 50, 1050, 1050), zoom = 0.6)
for (i in seq_len(nrow(comets_kuiper$normal))) {
  add_sphere(x = radius * comets_kuiper$normal[i, 1],
             y = radius * comets_kuiper$normal[i, 2],
             z = radius * comets_kuiper$normal[i, 3], radius = 0.03,
             col = "black", lit = FALSE)
}
arrow3d(p0 = c(0, 0, 0), p1 = radius * mu_hat_kuiper,
        type = "lines", s = 0.1, width = 2, theta = 0.15, lit = FALSE,
        lwd = 5, col = 2)
add_sphere_cap(x = 0, y = 0, z = 0, radius = 1.01 * radius,
               mu = mu_hat_kuiper, theta = acos(r_mu_boot_kuiper),
               col = "red", lit = FALSE, alpha = 0.3, border = TRUE)
snapshot3d("figs/fig_1b_points.png", webshot = FALSE,
           width = 1000, height = 1000, top = FALSE)
close3d()

}

## Plots ecdf
{

# Plot projected ecdf and cdf
plot_ecdf_fit <- function(X, k, mu_hat, rho_hat) {

  # Projections
  proj <- drop(X %*% mu_hat)

  # Ecdf
  Fn <- ecdf(proj)
  x_jump <- sort(unique(proj))
  y_jump <- Fn(x_jump)

  ggplot() +
    geom_step(
      data = data.frame(x = x_jump, y = y_jump),
      aes(x = x, y = y),
      linewidth = 0.9
    ) +
    geom_point(
      data = data.frame(x = x_jump, y = y_jump),
      aes(x = x, y = y),
      size = 1.4
    ) +
    geom_rug(
      data = data.frame(proj = proj),
      aes(x = proj),
      sides = "b",
      alpha = 0.4
    ) +
    stat_function(
      fun = function(x) {
        p_proj_car_gamma(
          x,
          mu = mu_hat,
          rho = rho_hat,
          k = k,
          gamma = mu_hat,
          p = 3
        )
      },
      linewidth = 1.1,
      color = "red"
    ) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1)) +
    labs(
      x = expression(bold(X)[i]^T * hat(bold(mu))[ML]),
      y = "Probability"
    ) +
    theme_minimal()
}

# Oort
mu_hat_oor <- gof_oor_2$theta_hat[-c(1, 2)]
rho_hat_oor <- gof_oor_2$theta_hat[2]
g_oor <- plot_ecdf_fit(
  X = comets_oort$normal,
  k = 2,
  mu_hat = mu_hat_oor,
  rho_hat = rho_hat_oor
)
g_oor

# Kuiper
mu_hat_kui <- gof_kui_1$theta_hat[-c(1, 2)]
rho_hat_kui <- gof_kui_1$theta_hat[2]
g_kui <- plot_ecdf_fit(
  X = comets_kuiper$normal,
  k = 1,
  mu_hat = mu_hat_kui,
  rho_hat = rho_hat_kui
)
g_kui

# Tests
goftest::cvm.test(comets_oort$normal %*% mu_hat_oor,
                  p_proj_car_gamma, mu = mu_hat_oor, rho = rho_hat_oor, k = 2,
                  gamma = mu_hat_oor, p = 3)
goftest::cvm.test(comets_kuiper$normal %*% mu_hat_kui,
                  p_proj_car_gamma, mu = mu_hat_kui, rho = rho_hat_kui, k = 2,
                  gamma = mu_hat_kui, p = 3)

# Save the plots
ggsave("figs/ecdf_oort.pdf", plot = g_oor, width = 7, height = 7)
ggsave("figs/ecdf_oort.png", plot = g_oor, width = 7, height = 7)
ggsave("figs/ecdf_kui.pdf", plot = g_kui, width = 7, height = 7)
ggsave("figs/ecdf_kui.png", plot = g_kui, width = 7, height = 7)

}
