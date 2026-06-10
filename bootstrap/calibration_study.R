# Empirical calibration study for multiplier bootstrap GOF tests

# Avoid repeated renv autoloader warnings in parallel workers for long-running
# simulation studies. This does not modify the project library; it only stops
# worker sessions from re-checking renv state on startup.
Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_calibration_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

multiplier_bootstrap_path_calibration <- resolve_calibration_path(
  "bootstrap",
  "multiplier_bootstrap.R"
)
required_bootstrap_functions <- c(
  "multiplier_bootstrap_normal",
  "multiplier_bootstrap_vmf",
  "multiplier_bootstrap_jp",
  "multiplier_bootstrap_hvmf",
  "multiplier_bootstrap_logistic_gaussian",
  "multiplier_bootstrap_beta_mixture2",
  "multiplier_bootstrap_logitnormal_mixture2",
  "multiplier_bootstrap_cardioid",
  "multiplier_bootstrap_spherical_cauchy",
  "multiplier_bootstrap_small_circle"
)
if (any(!vapply(required_bootstrap_functions, exists, logical(1), mode = "function"))) {
  source(multiplier_bootstrap_path_calibration)
}

default_calibration_output_dir <- function(mode = c(
                                            "smoke",
                                            "full",
                                            "smoke_composite",
                                            "full_composite"
                                          )) {
  mode <- match.arg(mode)
  file.path("output", "calibration", "bootstrap", mode)
}

format_calibration_number_tag <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x)) {
    stop("`x` must be a finite scalar.")
  }

  tag <- format(x, scientific = FALSE, trim = TRUE)
  tag <- gsub("-", "neg", tag, fixed = TRUE)
  gsub("\\.", "p", tag)
}

make_normal_simple_calibration_scenario <- function() {
  omega_max <- stats::qnorm(0.995)
  list(
    id = "normal_simple_n01",
    model = "normal",
    label = "Normal simple: N(0,1)",
    null = list(type = "simple", theta = list(mu = 0, sigma = 1)),
    sample_params = list(mu = 0, sigma = 1),
    ks_grid = list(
      omega_grid = seq(-omega_max, omega_max, length.out = 10),
      t_grid = seq(0, omega_max, length.out = 10)
    )
  )
}

make_normal_composite_calibration_scenario <- function() {
  omega_max <- stats::qnorm(0.995)
  list(
    id = "normal_composite_both_n01",
    model = "normal",
    label = "Normal composite (mu, sigma): N(0,1)",
    null = list(type = "composite"),
    sample_params = list(mu = 0, sigma = 1),
    unknown_param = "both",
    ks_grid = list(
      omega_grid = seq(-omega_max, omega_max, length.out = 10),
      t_grid = seq(0, omega_max, length.out = 10)
    )
  )
}

default_logistic_gaussian_omega_grid <- function(ambient_dim,
                                                 dominant_mass = 0.6) {
  ambient_dim <- as.integer(ambient_dim)
  if (length(ambient_dim) != 1L || !is.finite(ambient_dim) || ambient_dim < 2L) {
    stop("`ambient_dim` must be an integer greater than or equal to 2.")
  }
  dominant_mass <- as.numeric(dominant_mass)
  if (length(dominant_mass) != 1L || !is.finite(dominant_mass)) {
    stop("`dominant_mass` must be a finite scalar.")
  }
  if (dominant_mass <= 1 / ambient_dim || dominant_mass >= 1) {
    stop("`dominant_mass` must belong to (1 / ambient_dim, 1).")
  }

  barycenter <- rep.int(1 / ambient_dim, ambient_dim)
  off_mass <- (1 - dominant_mass) / (ambient_dim - 1L)
  corners <- t(vapply(seq_len(ambient_dim), function(i) {
    point <- rep.int(off_mass, ambient_dim)
    point[[i]] <- dominant_mass
    point
  }, numeric(ambient_dim)))

  rbind(barycenter, corners)
}

make_logistic_gaussian_ks_grid <- function(mu_ilr,
                                           Sigma_ilr,
                                           ambient_dim = length(mu_ilr) + 1L,
                                           dominant_mass = 0.6,
                                           n_t = 10L) {
  theta <- normalize_logistic_gaussian_theta(
    list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    ambient_dim = ambient_dim
  )
  omega_grid <- default_logistic_gaussian_omega_grid(
    ambient_dim = theta$ambient_dim,
    dominant_mass = dominant_mass
  )
  omega_ilr <- logistic_gaussian_ilr_matrix(omega_grid)
  shift_matrix <- matrix(
    rep(theta$mu_ilr, each = nrow(omega_ilr)),
    nrow = nrow(omega_ilr),
    ncol = length(theta$mu_ilr)
  ) - omega_ilr
  shift_norm_max <- max(sqrt(rowSums(shift_matrix^2)))
  lambda_max <- max(theta$eigenvalues_full)
  t_max <- shift_norm_max + 5 * sqrt(max(lambda_max, 0) * theta$ilr_dim)
  if (!is.finite(t_max) || t_max <= 0) {
    t_max <- 1
  }

  list(
    omega_grid = omega_grid,
    t_grid = seq(0, t_max, length.out = as.integer(n_t))
  )
}

make_logistic_gaussian_simple_calibration_scenario <- function(id,
                                                               label,
                                                               mu_ilr,
                                                               Sigma_ilr,
                                                               dominant_mass = 0.6,
                                                               n_t = 10L) {
  ambient_dim <- length(mu_ilr) + 1L
  list(
    id = id,
    model = "logistic_gaussian",
    label = label,
    null = list(type = "simple", theta = list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr)),
    sample_params = list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    distance_type = "aitchison",
    unknown_param = "both",
    ks_grid = make_logistic_gaussian_ks_grid(
      mu_ilr = mu_ilr,
      Sigma_ilr = Sigma_ilr,
      ambient_dim = ambient_dim,
      dominant_mass = dominant_mass,
      n_t = n_t
    )
  )
}

make_logistic_gaussian_composite_calibration_scenario <- function(id,
                                                                  label,
                                                                  mu_ilr,
                                                                  Sigma_ilr,
                                                                  dominant_mass = 0.6,
                                                                  n_t = 10L) {
  ambient_dim <- length(mu_ilr) + 1L
  list(
    id = id,
    model = "logistic_gaussian",
    label = label,
    null = list(type = "composite"),
    sample_params = list(mu_ilr = mu_ilr, Sigma_ilr = Sigma_ilr),
    distance_type = "aitchison",
    unknown_param = "both",
    ks_grid = make_logistic_gaussian_ks_grid(
      mu_ilr = mu_ilr,
      Sigma_ilr = Sigma_ilr,
      ambient_dim = ambient_dim,
      dominant_mass = dominant_mass,
      n_t = n_t
    )
  )
}

default_logistic_gaussian_simple_calibration_scenarios <- function() {
  list(
    make_logistic_gaussian_simple_calibration_scenario(
      id = "lg_simple_delta2_balanced",
      label = "Logistic Gaussian simple Delta^2 Aitchison: balanced",
      mu_ilr = c(0, 0),
      Sigma_ilr = matrix(c(0.35, 0.12, 0.12, 0.25), nrow = 2L)
    ),
    make_logistic_gaussian_simple_calibration_scenario(
      id = "lg_simple_delta2_shifted",
      label = "Logistic Gaussian simple Delta^2 Aitchison: shifted",
      mu_ilr = c(0.9, -0.5),
      Sigma_ilr = matrix(c(0.60, -0.18, -0.18, 0.45), nrow = 2L)
    )
  )
}

default_logistic_gaussian_composite_calibration_scenarios <- function() {
  list(
    make_logistic_gaussian_composite_calibration_scenario(
      id = "lg_composite_delta2_balanced",
      label = "Logistic Gaussian composite Delta^2 Aitchison: balanced",
      mu_ilr = c(0, 0),
      Sigma_ilr = matrix(c(0.35, 0.12, 0.12, 0.25), nrow = 2L)
    ),
    make_logistic_gaussian_composite_calibration_scenario(
      id = "lg_composite_delta2_shifted",
      label = "Logistic Gaussian composite Delta^2 Aitchison: shifted",
      mu_ilr = c(0.9, -0.5),
      Sigma_ilr = matrix(c(0.60, -0.18, -0.18, 0.45), nrow = 2L)
    )
  )
}

make_vmf_simple_calibration_scenario <- function(kappa) {
  mu <- c(1, 0, 0)
  list(
    id = sprintf("vmf_simple_s2_geodesic_kappa_%0.1f", kappa),
    model = "vmf",
    label = sprintf("vMF simple S^2 geodesic: kappa=%0.1f", kappa),
    null = list(type = "simple", theta = list(mu = mu, kappa = kappa)),
    sample_params = list(mu = mu, kappa = kappa),
    distance_type = "geodesic",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_vmf_composite_calibration_scenario <- function(kappa) {
  mu <- c(1, 0, 0)
  list(
    id = sprintf("vmf_composite_s2_geodesic_kappa_%0.1f", kappa),
    model = "vmf",
    label = sprintf("vMF composite S^2 geodesic: kappa=%0.1f", kappa),
    null = list(type = "composite"),
    sample_params = list(mu = mu, kappa = kappa),
    distance_type = "geodesic",
    unknown_param = "xi",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_jp_simple_calibration_scenario <- function(kappa, psi) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "jp_simple_s2_geodesic_kappa_%s_psi_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(psi)
    ),
    model = "jp",
    label = sprintf("JP simple S^2 geodesic: kappa=%s, psi=%s", kappa, psi),
    null = list(type = "simple", theta = list(mu = mu, kappa = kappa, psi = psi)),
    sample_params = list(mu = mu, kappa = kappa, psi = psi),
    distance_type = "geodesic",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_jp_composite_calibration_scenario <- function(kappa, psi) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "jp_composite_s2_geodesic_kappa_%s_psi_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(psi)
    ),
    model = "jp",
    label = sprintf("JP composite S^2 geodesic: kappa=%s, psi=%s", kappa, psi),
    null = list(type = "composite"),
    sample_params = list(mu = mu, kappa = kappa, psi = psi),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_cardioid_simple_calibration_scenario <- function(k,
                                                      rho = 0.5) {
  mu <- c(0, 0, 1)
  k <- as.integer(k)

  list(
    id = sprintf(
      "cardioid_simple_s2_geodesic_k%d_rho_%s",
      k,
      format_calibration_number_tag(rho)
    ),
    model = "cardioid",
    label = sprintf("Cardioide simple S^2 geodesic: k=%d, rho=%s", k, rho),
    null = list(type = "simple", theta = list(mu = mu, rho = rho, k = k)),
    sample_params = list(mu = mu, rho = rho, k = k),
    distance_type = "geodesic",
    unknown_param = "both",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_cardioid_composite_calibration_scenario <- function(k,
                                                          rho = 0.5) {
  mu <- c(0, 0, 1)
  k <- as.integer(k)

  list(
    id = sprintf(
      "cardioid_composite_s2_geodesic_k%d_rho_%s",
      k,
      format_calibration_number_tag(rho)
    ),
    model = "cardioid",
    label = sprintf("Cardioide compuesta S^2 geodesic: k=%d, rho=%s", k, rho),
    null = list(type = "composite"),
    sample_params = list(mu = mu, rho = rho, k = k),
    distance_type = "geodesic",
    unknown_param = "both",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_spherical_cauchy_simple_calibration_scenario <- function(rho = 0.5) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "spherical_cauchy_simple_s2_geodesic_rho_%s",
      format_calibration_number_tag(rho)
    ),
    model = "spherical_cauchy",
    label = sprintf("Spherical Cauchy simple S^2 geodesic: rho=%s", rho),
    null = list(type = "simple", theta = list(mu = mu, rho = rho)),
    sample_params = list(mu = mu, rho = rho),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_spherical_cauchy_composite_calibration_scenario <- function(rho = 0.5) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "spherical_cauchy_composite_s2_geodesic_rho_%s",
      format_calibration_number_tag(rho)
    ),
    model = "spherical_cauchy",
    label = sprintf("Spherical Cauchy composite S^2 geodesic: rho=%s", rho),
    null = list(type = "composite"),
    sample_params = list(mu = mu, rho = rho),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_beta_mixture2_simple_calibration_scenario <- function(weight1 = 0.4,
                                                                      alpha1 = 2,
                                                                      beta1 = 8,
                                                                      alpha2 = 8,
                                                                      beta2 = 2) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "beta_mixture2_simple_s2_geodesic_w%s_a1_%s_b1_%s_a2_%s_b2_%s",
      format_calibration_number_tag(weight1),
      format_calibration_number_tag(alpha1),
      format_calibration_number_tag(beta1),
      format_calibration_number_tag(alpha2),
      format_calibration_number_tag(beta2)
    ),
    model = "beta_mixture2",
    label = sprintf(
      "Rotational beta-mixture simple S^2 geodesic: w=%s, (a1,b1)=(%s,%s), (a2,b2)=(%s,%s)",
      weight1,
      alpha1,
      beta1,
      alpha2,
      beta2
    ),
    null = list(
      type = "simple",
      theta = list(
        mu = mu,
        weight1 = weight1,
        alpha1 = alpha1,
        beta1 = beta1,
        alpha2 = alpha2,
        beta2 = beta2
      )
    ),
    sample_params = list(
      mu = mu,
      weight1 = weight1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha2 = alpha2,
      beta2 = beta2
    ),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      beta_mixture2_L_max = 150L,
      beta_mixture2_quad_n = 400L,
      beta_mixture2_tol = 1e-10,
      beta_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_beta_mixture2_composite_calibration_scenario <- function(weight1 = 0.4,
                                                                         alpha1 = 2,
                                                                         beta1 = 8,
                                                                         alpha2 = 8,
                                                                         beta2 = 2) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "beta_mixture2_composite_s2_geodesic_w%s_a1_%s_b1_%s_a2_%s_b2_%s",
      format_calibration_number_tag(weight1),
      format_calibration_number_tag(alpha1),
      format_calibration_number_tag(beta1),
      format_calibration_number_tag(alpha2),
      format_calibration_number_tag(beta2)
    ),
    model = "beta_mixture2",
    label = sprintf(
      "Rotational beta-mixture composite S^2 geodesic: w=%s, (a1,b1)=(%s,%s), (a2,b2)=(%s,%s)",
      weight1,
      alpha1,
      beta1,
      alpha2,
      beta2
    ),
    null = list(type = "composite"),
    sample_params = list(
      mu = mu,
      weight1 = weight1,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha2 = alpha2,
      beta2 = beta2
    ),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      beta_mixture2_L_max = 150L,
      beta_mixture2_quad_n = 400L,
      beta_mixture2_tol = 1e-10,
      beta_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_logitnormal_mixture2_simple_calibration_scenario <- function(weight1 = 0.45,
                                                                             mean1 = -1.2,
                                                                             sd1 = 0.5,
                                                                             mean2 = 1.0,
                                                                             sd2 = 0.6) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "logitnormal_mixture2_simple_s2_geodesic_w%s_m1_%s_s1_%s_m2_%s_s2_%s",
      format_calibration_number_tag(weight1),
      format_calibration_number_tag(mean1),
      format_calibration_number_tag(sd1),
      format_calibration_number_tag(mean2),
      format_calibration_number_tag(sd2)
    ),
    model = "logitnormal_mixture2",
    label = sprintf(
      "Rotational logit-normal-mixture simple S^2 geodesic: w=%s, (m1,s1)=(%s,%s), (m2,s2)=(%s,%s)",
      weight1,
      mean1,
      sd1,
      mean2,
      sd2
    ),
    null = list(
      type = "simple",
      theta = list(
        mu = mu,
        weight1 = weight1,
        mean1 = mean1,
        sd1 = sd1,
        mean2 = mean2,
        sd2 = sd2
      )
    ),
    sample_params = list(
      mu = mu,
      weight1 = weight1,
      mean1 = mean1,
      sd1 = sd1,
      mean2 = mean2,
      sd2 = sd2
    ),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      logitnormal_mixture2_L_max = 150L,
      logitnormal_mixture2_quad_n = 400L,
      logitnormal_mixture2_tol = 1e-10,
      logitnormal_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_logitnormal_mixture2_composite_calibration_scenario <- function(weight1 = 0.45,
                                                                                mean1 = -1.2,
                                                                                sd1 = 0.5,
                                                                                mean2 = 1.0,
                                                                                sd2 = 0.6) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "logitnormal_mixture2_composite_s2_geodesic_w%s_m1_%s_s1_%s_m2_%s_s2_%s",
      format_calibration_number_tag(weight1),
      format_calibration_number_tag(mean1),
      format_calibration_number_tag(sd1),
      format_calibration_number_tag(mean2),
      format_calibration_number_tag(sd2)
    ),
    model = "logitnormal_mixture2",
    label = sprintf(
      "Rotational logit-normal-mixture composite S^2 geodesic: w=%s, (m1,s1)=(%s,%s), (m2,s2)=(%s,%s)",
      weight1,
      mean1,
      sd1,
      mean2,
      sd2
    ),
    null = list(type = "composite"),
    sample_params = list(
      mu = mu,
      weight1 = weight1,
      mean1 = mean1,
      sd1 = sd1,
      mean2 = mean2,
      sd2 = sd2
    ),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      logitnormal_mixture2_L_max = 150L,
      logitnormal_mixture2_quad_n = 400L,
      logitnormal_mixture2_tol = 1e-10,
      logitnormal_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_small_circle_symmetric_mixture2_simple_calibration_scenario <- function(kappa = 12,
                                                                             nu = 0.6) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "small_circle_symmetric_mixture2_simple_s2_geodesic_w50_50_kappa_%s_nu_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(nu)
    ),
    model = "small_circle_symmetric_mixture2",
    label = sprintf(
      "Symmetric two-small-circles simple S^2 geodesic: w=(0.5,0.5), kappa=%s, nu=%s",
      kappa,
      nu
    ),
    null = list(type = "simple", theta = list(mu = mu, kappa = kappa, nu = nu)),
    sample_params = list(mu = mu, kappa = kappa, nu = nu),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      small_circle_symmetric_mixture2_L_max = 150L,
      small_circle_symmetric_mixture2_quad_n = 400L,
      small_circle_symmetric_mixture2_tol = 1e-10,
      small_circle_symmetric_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_small_circle_symmetric_mixture2_composite_calibration_scenario <- function(kappa = 12,
                                                                                nu = 0.6) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "small_circle_symmetric_mixture2_composite_s2_geodesic_w50_50_kappa_%s_nu_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(nu)
    ),
    model = "small_circle_symmetric_mixture2",
    label = sprintf(
      "Symmetric two-small-circles composite S^2 geodesic: w=(0.5,0.5), kappa=%s, nu=%s",
      kappa,
      nu
    ),
    null = list(type = "composite"),
    sample_params = list(mu = mu, kappa = kappa, nu = nu),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      small_circle_symmetric_mixture2_L_max = 150L,
      small_circle_symmetric_mixture2_quad_n = 400L,
      small_circle_symmetric_mixture2_tol = 1e-10,
      small_circle_symmetric_mixture2_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

make_small_circle_simple_calibration_scenario <- function(kappa,
                                                          nu) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "small_circle_simple_s2_geodesic_kappa_%s_nu_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(nu)
    ),
    model = "small_circle",
    label = sprintf("Small Circle simple S^2 geodesic: kappa=%s, nu=%s", kappa, nu),
    null = list(type = "simple", theta = list(mu = mu, kappa = kappa, nu = nu)),
    sample_params = list(mu = mu, kappa = kappa, nu = nu),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    )
  )
}

make_small_circle_composite_calibration_scenario <- function(kappa,
                                                             nu) {
  mu <- c(0, 0, 1)
  list(
    id = sprintf(
      "small_circle_composite_s2_geodesic_kappa_%s_nu_%s",
      format_calibration_number_tag(kappa),
      format_calibration_number_tag(nu)
    ),
    model = "small_circle",
    label = sprintf("Small Circle composite S^2 geodesic: kappa=%s, nu=%s", kappa, nu),
    null = list(type = "composite"),
    sample_params = list(mu = mu, kappa = kappa, nu = nu),
    distance_type = "geodesic",
    unknown_param = "theta",
    ks_grid = list(
      omega_grid = generate_canonical_lattice(10, dim = 3),
      t_grid = seq(1e-8, pi - 1e-8, length.out = 10)
    ),
    control = list(
      small_circle_L_max = 200L,
      small_circle_quad_n = 400L,
      small_circle_tol = 1e-10,
      small_circle_optim_control = list(maxit = 300L, reltol = 1e-9)
    )
  )
}

default_cardioid_simple_calibration_scenarios <- function(rho = 0.5) {
  list(
    make_cardioid_simple_calibration_scenario(k = 1L, rho = rho),
    make_cardioid_simple_calibration_scenario(k = 2L, rho = rho)
  )
}

default_cardioid_composite_calibration_scenarios <- function(rho = 0.5) {
  list(
    make_cardioid_composite_calibration_scenario(k = 1L, rho = rho),
    make_cardioid_composite_calibration_scenario(k = 2L, rho = rho)
  )
}

default_spherical_cauchy_simple_calibration_scenarios <- function(rho_values = c(0.3, 0.7)) {
  lapply(as.numeric(rho_values), make_spherical_cauchy_simple_calibration_scenario)
}

default_spherical_cauchy_composite_calibration_scenarios <- function(rho_values = c(0.3, 0.7)) {
  lapply(as.numeric(rho_values), make_spherical_cauchy_composite_calibration_scenario)
}

default_beta_mixture2_simple_calibration_scenarios <- function() {
  list(
    make_beta_mixture2_simple_calibration_scenario(
      weight1 = 0.4,
      alpha1 = 2,
      beta1 = 8,
      alpha2 = 8,
      beta2 = 2
    ),
    make_beta_mixture2_simple_calibration_scenario(
      weight1 = 0.55,
      alpha1 = 4,
      beta1 = 12,
      alpha2 = 10,
      beta2 = 3
    )
  )
}

default_beta_mixture2_composite_calibration_scenarios <- function() {
  list(
    make_beta_mixture2_composite_calibration_scenario(
      weight1 = 0.4,
      alpha1 = 2,
      beta1 = 8,
      alpha2 = 8,
      beta2 = 2
    ),
    make_beta_mixture2_composite_calibration_scenario(
      weight1 = 0.55,
      alpha1 = 4,
      beta1 = 12,
      alpha2 = 10,
      beta2 = 3
    )
  )
}

default_logitnormal_mixture2_simple_calibration_scenarios <- function() {
  list(
    make_logitnormal_mixture2_simple_calibration_scenario(
      weight1 = 0.45,
      mean1 = -1.2,
      sd1 = 0.5,
      mean2 = 1.0,
      sd2 = 0.6
    ),
    make_logitnormal_mixture2_simple_calibration_scenario(
      weight1 = 0.35,
      mean1 = -0.7,
      sd1 = 0.35,
      mean2 = 1.4,
      sd2 = 0.45
    )
  )
}

default_logitnormal_mixture2_composite_calibration_scenarios <- function() {
  list(
    make_logitnormal_mixture2_composite_calibration_scenario(
      weight1 = 0.45,
      mean1 = -1.2,
      sd1 = 0.5,
      mean2 = 1.0,
      sd2 = 0.6
    ),
    make_logitnormal_mixture2_composite_calibration_scenario(
      weight1 = 0.35,
      mean1 = -0.7,
      sd1 = 0.35,
      mean2 = 1.4,
      sd2 = 0.45
    )
  )
}

default_small_circle_symmetric_mixture2_simple_calibration_scenarios <- function(kappa_values = c(12, 20),
                                                                                 nu_values = c(0.35, 0.6)) {
  grid <- expand.grid(
    kappa = as.numeric(kappa_values),
    nu = as.numeric(nu_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  lapply(seq_len(nrow(grid)), function(i) {
    make_small_circle_symmetric_mixture2_simple_calibration_scenario(
      kappa = grid$kappa[[i]],
      nu = grid$nu[[i]]
    )
  })
}

default_small_circle_symmetric_mixture2_composite_calibration_scenarios <- function(kappa_values = c(12, 20),
                                                                                    nu_values = c(0.35, 0.6)) {
  grid <- expand.grid(
    kappa = as.numeric(kappa_values),
    nu = as.numeric(nu_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  lapply(seq_len(nrow(grid)), function(i) {
    make_small_circle_symmetric_mixture2_composite_calibration_scenario(
      kappa = grid$kappa[[i]],
      nu = grid$nu[[i]]
    )
  })
}

default_small_circle_simple_calibration_scenarios <- function(kappa_values = c(5, 20),
                                                              nu_values = c(0, 0.3, 0.7)) {
  grid <- expand.grid(
    kappa = as.numeric(kappa_values),
    nu = as.numeric(nu_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  lapply(seq_len(nrow(grid)), function(i) {
    make_small_circle_simple_calibration_scenario(
      kappa = grid$kappa[[i]],
      nu = grid$nu[[i]]
    )
  })
}

default_small_circle_composite_calibration_scenarios <- function(kappa_values = c(5, 20),
                                                                 nu_values = c(0, 0.3, 0.7)) {
  grid <- expand.grid(
    kappa = as.numeric(kappa_values),
    nu = as.numeric(nu_values),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  lapply(seq_len(nrow(grid)), function(i) {
    make_small_circle_composite_calibration_scenario(
      kappa = grid$kappa[[i]],
      nu = grid$nu[[i]]
    )
  })
}

default_jp_simple_calibration_scenarios <- function() {
  list(
    make_jp_simple_calibration_scenario(1.0, 0.5),
    make_jp_simple_calibration_scenario(2.0, -0.5)
  )
}

default_jp_composite_calibration_scenarios <- function() {
  list(
    make_jp_composite_calibration_scenario(1.0, 0.5),
    make_jp_composite_calibration_scenario(2.0, -0.5)
  )
}

default_hvmf_calibration_data_dir <- function() {
  candidates <- c(
    file.path("data", "hvmf_typeiv_calibration"),
    file.path("sim_data_hvmf", "hvmf_typeiv_calibration"),
    file.path("..", "data", "hvmf_typeiv_calibration"),
    file.path("..", "sim_data_hvmf", "hvmf_typeiv_calibration"),
    file.path("..", "..", "data", "hvmf_typeiv_calibration"),
    file.path("..", "..", "sim_data_hvmf", "hvmf_typeiv_calibration")
  )

  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(candidate)
    }
  }

  stop("Could not resolve HvMF calibration data directory in `data/` or `sim_data_hvmf/`.")
}

hvmf_typeiv_fixed_mu <- function() {
  t0 <- stats::qnorm(0.25, mean = 0, sd = 1 / 4)
  c(
    cosh(abs(t0)),
    sinh(abs(t0)) * sign(t0) / sqrt(2),
    sinh(abs(t0)) * sign(t0) / sqrt(2)
  )
}

make_hvmf_ks_grid <- function(data,
                              mu = hvmf_typeiv_fixed_mu(),
                              n_omega = 10L,
                              n_t = 10L,
                              chi_margin = 0.25,
                              t_margin = 0.25) {
  x <- normalize_hvmf_h2_data(data)
  mu <- normalize_hvmf_theta(list(mu = mu, kappa = 1), control = list())$mu
  n_omega <- as.integer(n_omega)
  n_t <- as.integer(n_t)

  if (!is.finite(n_omega) || n_omega < 2L) {
    stop("`n_omega` must be an integer >= 2.")
  }
  if (!is.finite(n_t) || n_t < 2L) {
    stop("`n_t` must be an integer >= 2.")
  }

  chi_values <- acosh(pmax(as.numeric(x[, 1L]), 1))
  chi_max <- max(chi_values, na.rm = TRUE) + as.numeric(chi_margin)
  if (!is.finite(chi_max) || chi_max <= 0) {
    chi_max <- max(acosh(mu[[1L]]), 0.25)
  }

  theta0 <- atan2(mu[[3L]], mu[[2L]])
  n_radial <- ceiling(n_omega / 2)
  radial_grid <- seq(0, chi_max, length.out = n_radial)
  angular_grid <- theta0 + c(0, pi / 2)

  omega_grid <- do.call(
    rbind,
    lapply(radial_grid, function(chi) {
      do.call(
        rbind,
        lapply(angular_grid, function(theta) {
          c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta))
        })
      )
    })
  )
  omega_grid <- omega_grid[seq_len(min(nrow(omega_grid), n_omega)), , drop = FALSE]

  distance_matrix <- hvmf_distance_matrix(omega_grid, x)
  t_max <- max(distance_matrix, na.rm = TRUE) + as.numeric(t_margin)
  if (!is.finite(t_max) || t_max <= 0) {
    t_max <- chi_max + as.numeric(t_margin)
  }

  list(
    omega_grid = omega_grid,
    t_grid = seq(1e-8, t_max, length.out = n_t)
  )
}

make_hvmf_composite_calibration_scenario <- function(kappa,
                                                     data_dir = default_hvmf_calibration_data_dir()) {
  mu <- hvmf_typeiv_fixed_mu()
  list(
    id = sprintf("hvmf_composite_h2_geodesic_kappa_%d", as.integer(kappa)),
    model = "hvmf",
    label = sprintf("HvMF composite H^2 geodesic: kappa=%d", as.integer(kappa)),
    null = list(type = "composite"),
    sample_params = list(mu = mu, kappa = as.numeric(kappa)),
    distance_type = "geodesic",
    unknown_param = "both",
    ks_grid = NULL,
    data_dir = data_dir
  )
}

make_hvmf_simple_calibration_scenario <- function(kappa,
                                                  data_dir = default_hvmf_calibration_data_dir()) {
  mu <- hvmf_typeiv_fixed_mu()
  list(
    id = sprintf("hvmf_simple_h2_geodesic_kappa_%d", as.integer(kappa)),
    model = "hvmf",
    label = sprintf("HvMF simple H^2 geodesic: kappa=%d", as.integer(kappa)),
    null = list(type = "simple", theta = list(mu = mu, kappa = as.numeric(kappa))),
    sample_params = list(mu = mu, kappa = as.numeric(kappa)),
    distance_type = "geodesic",
    unknown_param = "both",
    ks_grid = NULL,
    data_dir = data_dir
  )
}

default_hvmf_simple_calibration_scenarios <- function(data_dir = default_hvmf_calibration_data_dir()) {
  list(
    make_hvmf_simple_calibration_scenario(50, data_dir = data_dir),
    make_hvmf_simple_calibration_scenario(200, data_dir = data_dir)
  )
}

default_hvmf_composite_calibration_scenarios <- function(data_dir = default_hvmf_calibration_data_dir()) {
  list(
    make_hvmf_composite_calibration_scenario(50, data_dir = data_dir),
    make_hvmf_composite_calibration_scenario(200, data_dir = data_dir)
  )
}

default_bootstrap_calibration_scenarios <- function() {
  list(
    make_normal_simple_calibration_scenario(),
    make_vmf_simple_calibration_scenario(0.5),
    make_vmf_simple_calibration_scenario(2.0),
    default_beta_mixture2_simple_calibration_scenarios()[[1L]],
    default_logitnormal_mixture2_simple_calibration_scenarios()[[1L]]
  )
}

default_bootstrap_composite_calibration_scenarios <- function() {
  list(
    make_normal_composite_calibration_scenario(),
    make_vmf_composite_calibration_scenario(0.5),
    make_vmf_composite_calibration_scenario(2.0),
    default_beta_mixture2_composite_calibration_scenarios()[[1L]],
    default_logitnormal_mixture2_composite_calibration_scenarios()[[1L]]
  )
}

list_hvmf_calibration_files <- function(scenario, n) {
  if (is.null(scenario$data_dir)) {
    stop("HvMF calibration scenario requires `data_dir`.")
  }

  folder <- file.path(
    scenario$data_dir,
    sprintf("kappa%d", as.integer(scenario$sample_params$kappa)),
    sprintf("n%d", as.integer(n))
  )
  if (!dir.exists(folder)) {
    stop(sprintf("HvMF calibration data folder not found: %s", folder))
  }

  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) {
    stop(sprintf("No HvMF calibration CSV files found in: %s", folder))
  }

  sample_ids <- suppressWarnings(as.integer(sub(".*samp_([0-9]+).*", "\\1", basename(files))))
  if (any(is.na(sample_ids))) {
    stop(sprintf("Could not parse HvMF sample ids from filenames in: %s", folder))
  }

  files[order(sample_ids)]
}

load_hvmf_calibration_sample <- function(scenario, n, replicate_id) {
  files <- list_hvmf_calibration_files(scenario, n)
  replicate_id <- as.integer(replicate_id)

  if (!is.finite(replicate_id) || replicate_id < 1L || replicate_id > length(files)) {
    stop(sprintf(
      "HvMF replicate_id=%d is out of bounds for n=%d; available range is 1..%d.",
      replicate_id,
      as.integer(n),
      length(files)
    ))
  }

  raw_data <- utils::read.csv(files[[replicate_id]])
  required_columns <- c("V1", "V2", "V3")
  if (!all(required_columns %in% names(raw_data))) {
    stop(sprintf(
      "HvMF calibration file %s does not contain columns %s.",
      files[[replicate_id]],
      paste(required_columns, collapse = ", ")
    ))
  }

  data_matrix <- as.matrix(raw_data[, required_columns, drop = FALSE])
  if (nrow(data_matrix) != as.integer(n)) {
    stop(sprintf(
      "HvMF calibration file %s has %d rows but n=%d was requested.",
      files[[replicate_id]],
      nrow(data_matrix),
      as.integer(n)
    ))
  }

  normalize_hvmf_h2_data(data_matrix)
}

simulate_h0_sample <- function(scenario, n, replicate_id = NULL) {
  if (identical(scenario$model, "normal")) {
    return(stats::rnorm(
      n = n,
      mean = scenario$sample_params$mu,
      sd = scenario$sample_params$sigma
    ))
  }

  if (identical(scenario$model, "vmf")) {
    if (!requireNamespace("rotasym", quietly = TRUE)) {
      stop("The `rotasym` package is required for vMF calibration studies.")
    }
    return(rotasym::r_vMF(
      n = n,
      mu = scenario$sample_params$mu,
      kappa = scenario$sample_params$kappa
    ))
  }

  if (identical(scenario$model, "jp")) {
    return(r_sph_jp(
      n = n,
      mu = scenario$sample_params$mu,
      kappa = scenario$sample_params$kappa,
      psi = scenario$sample_params$psi
    ))
  }

  if (identical(scenario$model, "hvmf")) {
    if (is.null(replicate_id)) {
      stop("HvMF calibration sampling from disk requires `replicate_id`.")
    }

    return(load_hvmf_calibration_sample(
      scenario = scenario,
      n = n,
      replicate_id = replicate_id
    ))
  }

  if (identical(scenario$model, "logistic_gaussian")) {
    return(rlogistic_gaussian_simplex(
      n = n,
      mu_ilr = scenario$sample_params$mu_ilr,
      Sigma_ilr = scenario$sample_params$Sigma_ilr
    ))
  }

  if (identical(scenario$model, "cardioid")) {
    return(r_sph_car(
      n = n,
      mu = scenario$sample_params$mu,
      rho = scenario$sample_params$rho,
      k = as.integer(scenario$sample_params$k)
    ))
  }

  if (identical(scenario$model, "spherical_cauchy")) {
    return(r_sph_spherical_cauchy(
      n = n,
      mu = scenario$sample_params$mu,
      rho = scenario$sample_params$rho
    ))
  }

  if (identical(scenario$model, "beta_mixture2")) {
    return(r_sph_beta_mixture2(
      n = n,
      mu = scenario$sample_params$mu,
      weight1 = scenario$sample_params$weight1,
      alpha1 = scenario$sample_params$alpha1,
      beta1 = scenario$sample_params$beta1,
      alpha2 = scenario$sample_params$alpha2,
      beta2 = scenario$sample_params$beta2
    ))
  }

  if (identical(scenario$model, "logitnormal_mixture2")) {
    return(r_sph_logitnormal_mixture2(
      n = n,
      mu = scenario$sample_params$mu,
      weight1 = scenario$sample_params$weight1,
      mean1 = scenario$sample_params$mean1,
      sd1 = scenario$sample_params$sd1,
      mean2 = scenario$sample_params$mean2,
      sd2 = scenario$sample_params$sd2
    ))
  }

  if (identical(scenario$model, "small_circle")) {
    return(r_sph_small_circle(
      n = n,
      mu = scenario$sample_params$mu,
      kappa = scenario$sample_params$kappa,
      nu = scenario$sample_params$nu
    ))
  }

  if (identical(scenario$model, "small_circle_symmetric_mixture2")) {
    return(r_sph_small_circle_symmetric_mixture2(
      n = n,
      mu = scenario$sample_params$mu,
      kappa = scenario$sample_params$kappa,
      nu = scenario$sample_params$nu
    ))
  }

  stop(sprintf("Unsupported scenario model: %s", scenario$model))
}

run_bootstrap_for_scenario <- function(data,
                                       scenario,
                                       B,
                                       alpha_nominal,
                                       seed,
                                       statistics = c("ks", "cvm")) {
  if (identical(scenario$model, "normal")) {
    return(multiplier_bootstrap_normal(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      unknown_param = scenario$unknown_param %||% NULL
    ))
  }

  if (identical(scenario$model, "vmf")) {
    return(multiplier_bootstrap_vmf(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      distance_type = scenario$distance_type,
      unknown_param = scenario$unknown_param %||% "xi"
    ))
  }

  if (identical(scenario$model, "jp")) {
    return(multiplier_bootstrap_jp(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      distance_type = scenario$distance_type,
      control = scenario$control %||% list()
    ))
  }

  if (identical(scenario$model, "hvmf")) {
    ks_grid <- scenario$ks_grid %||% if ("ks" %in% tolower(statistics)) {
      make_hvmf_ks_grid(
        data = data,
        mu = scenario$sample_params$mu
      )
    } else {
      NULL
    }

    return(multiplier_bootstrap_hvmf(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      unknown_param = scenario$unknown_param %||% "both"
    ))
  }

  if (identical(scenario$model, "logistic_gaussian")) {
    return(multiplier_bootstrap_logistic_gaussian(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      control = scenario$control %||% list(),
      unknown_param = scenario$unknown_param %||% "both"
    ))
  }

  if (identical(scenario$model, "beta_mixture2")) {
    return(multiplier_bootstrap_beta_mixture2(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic"
    ))
  }

  if (identical(scenario$model, "logitnormal_mixture2")) {
    return(multiplier_bootstrap_logitnormal_mixture2(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic"
    ))
  }

  if (identical(scenario$model, "cardioid")) {
    return(multiplier_bootstrap_cardioid(
      data = data,
      null = scenario$null,
      k = as.integer(scenario$sample_params$k),
      statistics = statistics,
      ks_grid = scenario$ks_grid %||% NULL,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic",
      unknown_param = scenario$unknown_param %||% "both"
    ))
  }

  if (identical(scenario$model, "spherical_cauchy")) {
    return(multiplier_bootstrap_spherical_cauchy(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic"
    ))
  }

  if (identical(scenario$model, "small_circle")) {
    return(multiplier_bootstrap_small_circle(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic"
    ))
  }

  if (identical(scenario$model, "small_circle_symmetric_mixture2")) {
    return(multiplier_bootstrap_small_circle_symmetric_mixture2(
      data = data,
      null = scenario$null,
      statistics = statistics,
      ks_grid = scenario$ks_grid,
      B = B,
      alpha = alpha_nominal,
      seed = seed,
      n_cores = 1,
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = TRUE
      ),
      control = scenario$control %||% list(),
      distance_type = scenario$distance_type %||% "geodesic"
    ))
  }

  stop(sprintf("Unsupported scenario model: %s", scenario$model))
}

alpha_to_suffix <- function(alpha_value) {
  gsub("\\.", "_", format(alpha_value, nsmall = 2))
}

build_pvalue_rejection_columns <- function(p_value, alphas) {
  output <- vector("list", length(alphas))
  names(output) <- sprintf("reject_p_%0.2f", alphas)

  for (i in seq_along(alphas)) {
    alpha_value <- alphas[[i]]
    column_name <- sprintf("reject_p_%s", alpha_to_suffix(alpha_value))
    output[[i]] <- as.integer(p_value <= alpha_value)
    names(output)[[i]] <- column_name
  }

  output
}

calibration_result_rows_from_bootstrap <- function(bootstrap_result,
                                                   scenario,
                                                   n,
                                                   replicate_id,
                                                   alpha_nominal,
                                                   alphas) {
  statistics <- names(bootstrap_result$inference)

  rows <- lapply(statistics, function(stat_name) {
    stat_inference <- bootstrap_result$inference[[stat_name]]
    pvalue_rejections <- build_pvalue_rejection_columns(stat_inference$p_value, alphas)

    row_df <- data.frame(
      model = scenario$model,
      scenario = scenario$id,
      scenario_label = scenario$label,
      null_type = scenario$null$type,
      unknown_param = scenario$unknown_param %||% NA_character_,
      n = n,
      statistic = stat_name,
      replicate_id = replicate_id,
      p_value = stat_inference$p_value,
      observed_statistic = stat_inference$observed,
      stringsAsFactors = FALSE
    )

    for (column_name in names(pvalue_rejections)) {
      row_df[[column_name]] <- pvalue_rejections[[column_name]]
    }

    # Backward-compatible aliases used by previous exploratory summaries.
    row_df$reject_0_01 <- row_df$reject_p_0_01 %||% NA_integer_
    row_df$reject_0_05 <- row_df$reject_p_0_05 %||% NA_integer_
    row_df$reject_0_10 <- row_df$reject_p_0_10 %||% NA_integer_

    row_df
  })

  do.call(rbind, rows)
}

build_calibration_tasks <- function(n_values,
                                    M_outer,
                                    seed = NULL) {
  task_grid <- expand.grid(
    n = as.integer(n_values),
    replicate_id = seq_len(M_outer),
    KEEP.OUT.ATTRS = FALSE
  )
  task_grid <- task_grid[order(task_grid$n, task_grid$replicate_id), , drop = FALSE]

  if (!is.null(seed)) {
    set.seed(seed)
  }
  task_grid$sample_seed <- sample.int(.Machine$integer.max, nrow(task_grid))
  task_grid$bootstrap_seed <- sample.int(.Machine$integer.max, nrow(task_grid))

  task_grid
}

build_calibration_tasks_single_n <- function(n_value,
                                             M_outer,
                                             seed = NULL) {
  task_grid <- data.frame(
    n = rep.int(as.integer(n_value), M_outer),
    replicate_id = seq_len(M_outer),
    stringsAsFactors = FALSE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }
  task_grid$sample_seed <- sample.int(.Machine$integer.max, nrow(task_grid))
  task_grid$bootstrap_seed <- sample.int(.Machine$integer.max, nrow(task_grid))

  task_grid
}

run_calibration_task <- function(task_row,
                                 scenario,
                                 B,
                                 alpha_nominal,
                                 alphas,
                                 statistics) {
  set.seed(task_row$sample_seed)
  data <- simulate_h0_sample(
    scenario = scenario,
    n = task_row$n,
    replicate_id = task_row$replicate_id
  )
  bootstrap_result <- run_bootstrap_for_scenario(
    data = data,
    scenario = scenario,
    B = B,
    alpha_nominal = alpha_nominal,
    seed = task_row$bootstrap_seed,
    statistics = statistics
  )

  calibration_result_rows_from_bootstrap(
    bootstrap_result = bootstrap_result,
    scenario = scenario,
    n = task_row$n,
    replicate_id = task_row$replicate_id,
    alpha_nominal = alpha_nominal,
    alphas = alphas
  )
}

run_calibration_task_chunk <- function(task_chunk,
                                       scenario,
                                       B,
                                       alpha_nominal,
                                       alphas,
                                       statistics) {
  chunk_rows <- lapply(seq_len(nrow(task_chunk)), function(i) {
    run_calibration_task(
      task_row = task_chunk[i, , drop = FALSE],
      scenario = scenario,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics
    )
  })

  do.call(rbind, chunk_rows)
}

split_task_grid <- function(task_grid, n_chunks) {
  indices <- split(seq_len(nrow(task_grid)), rep(seq_len(n_chunks), length.out = nrow(task_grid)))
  lapply(indices, function(idx) task_grid[idx, , drop = FALSE])
}

initialize_calibration_cluster <- function(n_workers) {
  n_workers <- max(1L, as.integer(n_workers))
  calibration_path_worker <- normalizePath(
    resolve_calibration_path("bootstrap", "calibration_study.R"),
    winslash = "/",
    mustWork = TRUE
  )

  cl <- parallel::makeCluster(n_workers)

  parallel::clusterExport(
    cl,
    c("calibration_path_worker"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    source(calibration_path_worker)
    NULL
  })

  worker_symbols <- c(
    "run_calibration_task_chunk",
    "run_calibration_task",
    "simulate_h0_sample",
    "run_bootstrap_for_scenario",
    "calibration_result_rows_from_bootstrap",
    "build_pvalue_rejection_columns",
    "alpha_to_suffix"
  )
  parallel::clusterExport(cl, worker_symbols, envir = environment())

  cl
}

run_calibration_task_grid <- function(task_grid,
                                      scenario,
                                      B,
                                      alpha_nominal,
                                      alphas,
                                      statistics,
                                      n_cores_outer = 1,
                                      cluster = NULL) {
  n_cores_outer <- max(1L, as.integer(n_cores_outer))

  if (n_cores_outer == 1L) {
    return(run_calibration_task_chunk(
      task_chunk = task_grid,
      scenario = scenario,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics
    ))
  }

  task_chunks <- split_task_grid(task_grid, min(n_cores_outer, nrow(task_grid)))
  cl <- cluster
  created_here <- FALSE
  if (is.null(cl)) {
    cl <- initialize_calibration_cluster(length(task_chunks))
    created_here <- TRUE
  }
  if (created_here) {
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  parallel::clusterExport(
    cl,
    c("scenario", "B", "alpha_nominal", "alphas", "statistics"),
    envir = environment()
  )

  chunk_results <- parallel::parLapply(cl, task_chunks, function(chunk) {
    run_calibration_task_chunk(
      task_chunk = chunk,
      scenario = scenario,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics
    )
  })

  do.call(rbind, chunk_results)
}

run_calibration_scenario <- function(scenario,
                                     n_values,
                                     M_outer,
                                     B,
                                     alpha_nominal = 0.05,
                                     alphas = c(0.01, 0.05, 0.10),
                                     statistics = c("ks", "cvm"),
                                     n_cores_outer = 1,
                                     seed = NULL,
                                     cluster = NULL,
                                     progress_callback = NULL,
                                     verbose = TRUE) {
  n_values <- as.integer(n_values)
  M_outer <- as.integer(M_outer)
  B <- as.integer(B)
  statistics <- normalize_requested_statistics(statistics)
  n_cores_outer <- max(1L, as.integer(n_cores_outer))
  scenario_results <- vector("list", length(n_values))

  for (i in seq_along(n_values)) {
    n_value <- as.integer(n_values[[i]])
    n_seed <- if (is.null(seed)) NULL else seed + 100000L * i

    if (identical(scenario$model, "hvmf")) {
      n_available <- length(list_hvmf_calibration_files(scenario, n_value))
      if (M_outer > n_available) {
        stop(sprintf(
          "HvMF calibration requested M_outer=%d for n=%d, but only %d disk replicates are available.",
          M_outer,
          n_value,
          n_available
        ))
      }
    }

    if (isTRUE(verbose)) {
      message(sprintf(
        "  Scenario '%s': starting n = %d (%d/%d), M_outer = %d, B = %d",
        scenario$id,
        n_value,
        i,
        length(n_values),
        M_outer,
        B
      ))
    }

    task_grid <- build_calibration_tasks_single_n(
      n_value = n_value,
      M_outer = M_outer,
      seed = n_seed
    )

    start_time_n <- Sys.time()
    scenario_results[[i]] <- run_calibration_task_grid(
      task_grid = task_grid,
      scenario = scenario,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics,
      n_cores_outer = n_cores_outer,
      cluster = cluster
    )
    elapsed_n <- as.numeric(difftime(Sys.time(), start_time_n, units = "secs"))
    if (isTRUE(verbose)) {
      message(sprintf(
        "  Scenario '%s': completed n = %d in %.1f seconds",
        scenario$id,
        n_value,
        elapsed_n
      ))
    }
    if (is.function(progress_callback)) {
      progress_callback(
        scenario = scenario,
        n_value = n_value,
        index_n = i,
        n_values = n_values,
        elapsed_seconds = elapsed_n
      )
    }
  }

  do.call(rbind, scenario_results)
}

compute_uniform_ks_distance <- function(p_values) {
  p_values <- sort(as.numeric(p_values))
  n <- length(p_values)
  if (n == 0) {
    return(NA_real_)
  }

  d_plus <- max((seq_len(n) / n) - p_values)
  d_minus <- max(p_values - ((seq_len(n) - 1) / n))
  max(d_plus, d_minus)
}

summarize_calibration_results <- function(raw_results,
                                          alphas = c(0.01, 0.05, 0.10)) {
  summary_groups <- split(
    raw_results,
    list(
      raw_results$model,
      raw_results$scenario,
      raw_results$scenario_label,
      raw_results$n,
      raw_results$statistic
    ),
    drop = TRUE
  )

  summary_rows <- list()

  for (group_df in summary_groups) {
    p_values <- group_df$p_value
    n_outer <- length(p_values)

    base_summary <- list(
      model = group_df$model[[1]],
      scenario = group_df$scenario[[1]],
      scenario_label = group_df$scenario_label[[1]],
      null_type = group_df$null_type[[1]],
      unknown_param = group_df$unknown_param[[1]],
      n = group_df$n[[1]],
      statistic = group_df$statistic[[1]],
      M_outer = n_outer,
      p_value_mean = mean(p_values),
      p_value_median = stats::median(p_values),
      p_value_q05 = as.numeric(stats::quantile(p_values, probs = 0.05, names = FALSE, type = 8)),
      p_value_q25 = as.numeric(stats::quantile(p_values, probs = 0.25, names = FALSE, type = 8)),
      p_value_q75 = as.numeric(stats::quantile(p_values, probs = 0.75, names = FALSE, type = 8)),
      p_value_q95 = as.numeric(stats::quantile(p_values, probs = 0.95, names = FALSE, type = 8)),
      ks_uniform_distance = compute_uniform_ks_distance(p_values)
    )

    for (alpha_value in alphas) {
      suffix <- alpha_to_suffix(alpha_value)
      p_col <- sprintf("reject_p_%s", suffix)
      rejection_rate_p <- mean(group_df[[p_col]])
      mc_se <- sqrt(alpha_value * (1 - alpha_value) / n_outer)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        model = base_summary$model,
        scenario = base_summary$scenario,
        scenario_label = base_summary$scenario_label,
        null_type = base_summary$null_type,
        unknown_param = base_summary$unknown_param,
        n = base_summary$n,
        statistic = base_summary$statistic,
        alpha = alpha_value,
        rejection_rate_p = rejection_rate_p,
        rejection_gap_p = rejection_rate_p - alpha_value,
        mc_se = mc_se,
        M_outer = base_summary$M_outer,
        p_value_mean = base_summary$p_value_mean,
        p_value_median = base_summary$p_value_median,
        p_value_q05 = base_summary$p_value_q05,
        p_value_q25 = base_summary$p_value_q25,
        p_value_q75 = base_summary$p_value_q75,
        p_value_q95 = base_summary$p_value_q95,
        ks_uniform_distance = base_summary$ks_uniform_distance,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, summary_rows)
}

plot_calibration_pvalue_ecdf <- function(raw_results_scenario) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The `ggplot2` package is required to generate calibration plots.")
  }

  raw_results_scenario <- raw_results_scenario[, , drop = FALSE]
  raw_results_scenario$n_label <- factor(
    paste0("n=", raw_results_scenario$n),
    levels = paste0("n=", sort(unique(raw_results_scenario$n)))
  )

  ggplot2::ggplot(raw_results_scenario, ggplot2::aes(x = p_value)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.9, color = "#1f5aa6") +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "#333333"
    ) +
    ggplot2::facet_grid(statistic ~ n_label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = expression(alpha),
      y = "Empirical probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(size = 11)
    )
}

plot_calibration_rejection_rates <- function(summary_results_scenario) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The `ggplot2` package is required to generate calibration plots.")
  }

  n_levels <- paste0("n=", sort(unique(summary_results_scenario$n)))

  long_results <- data.frame(
    scenario_label = summary_results_scenario$scenario_label,
    n_label = factor(paste0("n=", summary_results_scenario$n), levels = n_levels),
    statistic = summary_results_scenario$statistic,
    alpha = summary_results_scenario$alpha,
    rejection_rate = summary_results_scenario$rejection_rate_p,
    method = "p_value",
    stringsAsFactors = FALSE
  )

  group_counts <- stats::aggregate(
    alpha ~ n_label + statistic,
    data = long_results,
    FUN = length
  )
  names(group_counts)[names(group_counts) == "alpha"] <- "n_group_points"
  long_results <- merge(
    long_results,
    group_counts,
    by = c("n_label", "statistic"),
    all.x = TRUE,
    sort = FALSE
  )

  ggplot2::ggplot(
    long_results,
    ggplot2::aes(
      x = alpha,
      y = rejection_rate,
      color = statistic,
      group = statistic
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "#333333"
    ) +
    ggplot2::geom_line(
      data = long_results[long_results$n_group_points > 1L, , drop = FALSE],
      linewidth = 0.9
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ n_label) +
    ggplot2::coord_cartesian(xlim = c(0, max(long_results$alpha) * 1.05), ylim = c(0, max(long_results$alpha, long_results$rejection_rate) * 1.15)) +
    ggplot2::labs(
      x = expression(alpha),
      y = "Empirical rejection probability",
      color = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(size = 11)
    )
}

save_calibration_plots <- function(raw_results,
                                   summary_results,
                                   output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The `ggplot2` package is required to save calibration plots.")
  }

  scenarios <- unique(raw_results$scenario)
  saved_files <- list()

  for (scenario_id in scenarios) {
    raw_scenario <- raw_results[raw_results$scenario == scenario_id, , drop = FALSE]
    summary_scenario <- summary_results[summary_results$scenario == scenario_id, , drop = FALSE]

    ecdf_plot <- plot_calibration_pvalue_ecdf(raw_scenario)
    rejection_plot <- plot_calibration_rejection_rates(summary_scenario)

    ecdf_png <- file.path(output_dir, sprintf("%s_pvalue_ecdf.png", scenario_id))
    ecdf_pdf <- file.path(output_dir, sprintf("%s_pvalue_ecdf.pdf", scenario_id))
    rejection_png <- file.path(output_dir, sprintf("%s_rejection_rates.png", scenario_id))
    rejection_pdf <- file.path(output_dir, sprintf("%s_rejection_rates.pdf", scenario_id))

    ggplot2::ggsave(ecdf_png, plot = ecdf_plot, width = 10, height = 6, dpi = 300)
    ggplot2::ggsave(ecdf_pdf, plot = ecdf_plot, width = 10, height = 6)
    ggplot2::ggsave(rejection_png, plot = rejection_plot, width = 10, height = 6, dpi = 300)
    ggplot2::ggsave(rejection_pdf, plot = rejection_plot, width = 10, height = 6)

    saved_files[[scenario_id]] <- list(
      ecdf_png = ecdf_png,
      ecdf_pdf = ecdf_pdf,
      rejection_png = rejection_png,
      rejection_pdf = rejection_pdf
    )
  }

  saved_files
}

persist_calibration_outputs <- function(raw_results,
                                        summary_results,
                                        output_dir,
                                        completed_scenarios = NULL,
                                        progress_state = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  raw_csv <- file.path(output_dir, "bootstrap_calibration_raw.csv")
  summary_csv <- file.path(output_dir, "bootstrap_calibration_summary.csv")
  utils::write.csv(raw_results, raw_csv, row.names = FALSE)
  utils::write.csv(summary_results, summary_csv, row.names = FALSE)

  progress_path <- file.path(output_dir, "progress_status.txt")
  progress_lines <- c(
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("completed_scenarios: %d", length(completed_scenarios %||% character(0)))
  )
  if (!is.null(completed_scenarios) && length(completed_scenarios) > 0) {
    progress_lines <- c(
      progress_lines,
      sprintf("scenario_ids: %s", paste(completed_scenarios, collapse = ", "))
    )
  }
  if (!is.null(progress_state) && length(progress_state) > 0) {
    progress_lines <- c(progress_lines, progress_state)
  }
  writeLines(progress_lines, con = progress_path)

  list(
    raw_csv = raw_csv,
    summary_csv = summary_csv,
    progress_path = progress_path
  )
}

is_pid_alive <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || pid <= 0L) {
    return(FALSE)
  }

  probe <- tryCatch(
    suppressWarnings(system2("ps", args = c("-p", as.character(pid), "-o", "pid="), stdout = TRUE, stderr = TRUE)),
    error = function(...) character(0)
  )

  length(trimws(probe)) > 0L
}

acquire_calibration_run_lock <- function(output_dir) {
  lock_path <- file.path(output_dir, ".calibration_run.lock")
  pid <- Sys.getpid()

  if (file.exists(lock_path)) {
    lock_contents <- tryCatch(readLines(lock_path, warn = FALSE), error = function(...) character(0))
    lock_pid_line <- lock_contents[grepl("^pid=", lock_contents)]
    lock_pid <- if (length(lock_pid_line) > 0L) sub("^pid=", "", lock_pid_line[[1]]) else NA_character_

    if (is_pid_alive(lock_pid)) {
      stop(sprintf(
        paste0(
          "A calibration run is already active for this output directory. ",
          "output_dir='%s', active_pid=%s, lock_path='%s'. ",
          "Stop the existing run or choose a different output_dir."
        ),
        output_dir,
        lock_pid,
        lock_path
      ))
    }

    unlink(lock_path, force = TRUE)
  }

  created <- file.create(lock_path)
  if (!isTRUE(created)) {
    stop(sprintf("Unable to create calibration lock file at '%s'.", lock_path))
  }

  writeLines(
    c(
      sprintf("pid=%d", pid),
      sprintf("started_at=%s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      sprintf("output_dir=%s", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
    ),
    con = lock_path
  )

  lock_path
}

run_bootstrap_calibration_study <- function(scenarios = default_bootstrap_calibration_scenarios(),
                                            n_values = c(50, 100, 200),
                                            M_outer = 1000,
                                            B = 1000,
                                            alpha_nominal = 0.05,
                                            alphas = c(0.01, 0.05, 0.10),
                                            statistics = c("ks", "cvm"),
                                            n_cores_outer = 1,
                                            seed = 123,
                                            output_dir = NULL,
                                            show_progress = FALSE,
                                            verbose = TRUE) {
  if (is.null(output_dir)) {
    output_dir <- default_calibration_output_dir("full")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  lock_path <- acquire_calibration_run_lock(output_dir)
  on.exit(unlink(lock_path, force = TRUE), add = TRUE)

  shared_cluster <- NULL
  if (as.integer(n_cores_outer) > 1L) {
    shared_cluster <- initialize_calibration_cluster(as.integer(n_cores_outer))
    on.exit(parallel::stopCluster(shared_cluster), add = TRUE)
  }

  total_progress_steps <- length(scenarios) * length(as.integer(n_values))
  progress_counter <- 0L
  progress_bar <- NULL
  if (isTRUE(show_progress)) {
    progress_bar <- utils::txtProgressBar(min = 0, max = total_progress_steps, style = 3)
    on.exit(close(progress_bar), add = TRUE)
  }

  progress_state <- character(0)
  progress_callback <- if (isTRUE(show_progress)) {
    function(scenario, n_value, index_n, n_values, elapsed_seconds) {
      progress_counter <<- progress_counter + 1L
      utils::setTxtProgressBar(progress_bar, progress_counter)
      progress_state <<- c(
        sprintf("last_scenario: %s", scenario$id),
        sprintf("last_n: %d", as.integer(n_value)),
        sprintf("completed_blocks: %d/%d", progress_counter, total_progress_steps),
        sprintf("last_block_elapsed_seconds: %.3f", elapsed_seconds)
      )
    }
  } else {
    NULL
  }

  scenario_results <- vector("list", length(scenarios))
  plot_files <- list()

  for (i in seq_along(scenarios)) {
    scenario_seed <- if (is.null(seed)) NULL else seed + 1000L * i
    if (isTRUE(verbose)) {
      message(sprintf(
        "Running calibration scenario %d/%d: %s",
        i,
        length(scenarios),
        scenarios[[i]]$label
      ))
    }
    scenario_results[[i]] <- run_calibration_scenario(
      scenario = scenarios[[i]],
      n_values = n_values,
      M_outer = M_outer,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics,
      n_cores_outer = n_cores_outer,
      seed = scenario_seed,
      cluster = shared_cluster,
      progress_callback = progress_callback,
      verbose = verbose
    )
    raw_results_partial <- do.call(rbind, scenario_results[seq_len(i)])
    summary_results_partial <- summarize_calibration_results(raw_results_partial, alphas = alphas)
    persist_info <- persist_calibration_outputs(
      raw_results = raw_results_partial,
      summary_results = summary_results_partial,
      output_dir = output_dir,
      completed_scenarios = vapply(scenarios[seq_len(i)], `[[`, character(1), "id"),
      progress_state = progress_state
    )
    plot_files <- save_calibration_plots(
      raw_results = raw_results_partial,
      summary_results = summary_results_partial,
      output_dir = output_dir
    )
    saveRDS(
      list(
        raw_results = raw_results_partial,
        summary_results = summary_results_partial,
        completed_scenarios = vapply(scenarios[seq_len(i)], `[[`, character(1), "id"),
        config = list(
          scenarios = scenarios,
          n_values = n_values,
          M_outer = M_outer,
          B = B,
          alpha_nominal = alpha_nominal,
          alphas = alphas,
          statistics = statistics,
          n_cores_outer = n_cores_outer,
          seed = seed
        )
      ),
      file = file.path(output_dir, "checkpoint_run_object.rds")
    )
  }

  raw_results <- do.call(rbind, scenario_results)
  summary_results <- summarize_calibration_results(raw_results, alphas = alphas)

  persist_info <- persist_calibration_outputs(
    raw_results = raw_results,
    summary_results = summary_results,
    output_dir = output_dir,
    completed_scenarios = vapply(scenarios, `[[`, character(1), "id"),
    progress_state = c(progress_state, "status: finished")
  )
  plot_files <- save_calibration_plots(
    raw_results = raw_results,
    summary_results = summary_results,
    output_dir = output_dir
  )

  list(
    raw_results = raw_results,
    summary_results = summary_results,
    output_dir = output_dir,
    raw_csv = persist_info$raw_csv,
    summary_csv = persist_info$summary_csv,
    plot_files = plot_files,
    config = list(
      scenarios = scenarios,
      n_values = n_values,
      M_outer = M_outer,
      B = B,
      alpha_nominal = alpha_nominal,
      alphas = alphas,
      statistics = statistics,
      n_cores_outer = n_cores_outer,
      seed = seed
    )
  )
}

run_smoke_bootstrap_calibration_study <- function(output_dir = NULL,
                                                  n_cores_outer = 1,
                                                  seed = 123) {
  run_bootstrap_calibration_study(
    scenarios = default_bootstrap_calibration_scenarios(),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% default_calibration_output_dir("smoke")
  )
}

run_full_bootstrap_calibration_study <- function(output_dir = NULL,
                                                 n_cores_outer = 1,
                                                 seed = 123) {
  run_bootstrap_calibration_study(
    scenarios = default_bootstrap_calibration_scenarios(),
    n_values = c(50, 100, 200),
    M_outer = 1000,
    B = 1000,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% default_calibration_output_dir("full")
  )
}

run_smoke_bootstrap_composite_calibration_study <- function(output_dir = NULL,
                                                            n_cores_outer = 1,
                                                            seed = 123) {
  run_bootstrap_calibration_study(
    scenarios = default_bootstrap_composite_calibration_scenarios(),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% default_calibration_output_dir("smoke_composite")
  )
}

run_full_bootstrap_composite_calibration_study <- function(output_dir = NULL,
                                                           n_cores_outer = 1,
                                                           seed = 123,
                                                           M_outer = 1000,
                                                           B = 1000) {
  run_bootstrap_calibration_study(
    scenarios = default_bootstrap_composite_calibration_scenarios(),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% default_calibration_output_dir("full_composite")
  )
}

run_smoke_bootstrap_jp_composite_calibration_study <- function(output_dir = NULL,
                                                               n_cores_outer = 1,
                                                               seed = 123,
                                                               jp_mle_maxit = 30L,
                                                               jp_mle_reltol = 1e-5) {
  scenarios <- default_jp_composite_calibration_scenarios()
  jp_control <- list(
    jp_mle_maxit = as.integer(jp_mle_maxit),
    jp_mle_reltol = as.numeric(jp_mle_reltol)
  )
  scenarios <- lapply(scenarios, function(scenario) {
    scenario$control <- modifyList(scenario$control %||% list(), jp_control)
    scenario
  })

  run_bootstrap_calibration_study(
    scenarios = scenarios,
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "jp_composite_smoke")
  )
}

run_full_bootstrap_jp_composite_calibration_study <- function(output_dir = NULL,
                                                              n_cores_outer = 1,
                                                              seed = 123,
                                                              M_outer = 1000,
                                                              B = 1000,
                                                              jp_mle_maxit = 30L,
                                                              jp_mle_reltol = 1e-5) {
  scenarios <- default_jp_composite_calibration_scenarios()
  jp_control <- list(
    jp_mle_maxit = as.integer(jp_mle_maxit),
    jp_mle_reltol = as.numeric(jp_mle_reltol)
  )
  scenarios <- lapply(scenarios, function(scenario) {
    scenario$control <- modifyList(scenario$control %||% list(), jp_control)
    scenario
  })

  run_bootstrap_calibration_study(
    scenarios = scenarios,
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "jp_composite_full")
  )
}

run_full_bootstrap_jp_composite_nm_boot_local_calibration_study <- function(output_dir = NULL,
                                                                             n_cores_outer = 1,
                                                                             seed = 123,
                                                                             M_outer = 500,
                                                                             B = 500,
                                                                             jp_mle_maxit = 500L,
                                                                             jp_mle_reltol = 1e-10,
                                                                             jp_mle_max_abs_kappa_psi = 6,
                                                                             statistics = c("ks", "cvm"),
                                                                             n_values = c(50, 100, 200)) {
  scenarios <- default_jp_composite_calibration_scenarios()
  jp_control <- list(
    jp_mle_method = "Nelder-Mead",
    jp_mle_maxit = as.integer(jp_mle_maxit),
    jp_mle_reltol = as.numeric(jp_mle_reltol),
    jp_mle_max_abs_kappa_psi = as.numeric(jp_mle_max_abs_kappa_psi),
    jp_profile_method = "tabulated",
    jp_profile_n_u = 1025L,
    jp_profile_n_delta = 257L,
    jp_vmf_switch_abs_kappa_psi = 1e-3
  )

  scenarios <- lapply(scenarios, function(scenario) {
    scenario$control <- modifyList(scenario$control %||% list(), jp_control)
    scenario
  })

  run_bootstrap_calibration_study(
    scenarios = scenarios,
    n_values = as.integer(n_values),
    M_outer = as.integer(M_outer),
    B = as.integer(B),
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = statistics,
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "jp_composite_nm_boot_local_full")
  )
}

run_smoke_hvmf_composite_cvm_calibration_study <- function(output_dir = NULL,
                                                           n_cores_outer = 1,
                                                           seed = 123,
                                                           show_progress = FALSE,
                                                           verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_composite_calibration_scenarios(),
    n_values = 50,
    M_outer = 5,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = "cvm",
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_composite_cvm_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_hvmf_composite_cvm_calibration_study <- function(output_dir = NULL,
                                                          n_cores_outer = 1,
                                                          seed = 123,
                                                          M_outer = 1000,
                                                          B = 1000,
                                                          show_progress = FALSE,
                                                          verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_composite_calibration_scenarios(),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = "cvm",
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_composite_cvm_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_hvmf_simple_cvm_calibration_study <- function(output_dir = NULL,
                                                        n_cores_outer = 1,
                                                        seed = 123,
                                                        show_progress = FALSE,
                                                        verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_simple_calibration_scenarios(),
    n_values = 50,
    M_outer = 5,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = "cvm",
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_simple_cvm_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_hvmf_simple_cvm_calibration_study <- function(output_dir = NULL,
                                                       n_cores_outer = 1,
                                                       seed = 123,
                                                       M_outer = 1000,
                                                       B = 1000,
                                                       show_progress = FALSE,
                                                       verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_simple_calibration_scenarios(),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = "cvm",
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_simple_cvm_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_hvmf_simple_ks_cvm_calibration_study <- function(output_dir = NULL,
                                                           n_cores_outer = 1,
                                                           seed = 123,
                                                           show_progress = FALSE,
                                                           verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_simple_calibration_scenarios(),
    n_values = 50,
    M_outer = 5,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_simple_ks_cvm_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_hvmf_simple_ks_cvm_calibration_study <- function(output_dir = NULL,
                                                          n_cores_outer = 1,
                                                          seed = 123,
                                                          M_outer = 1000,
                                                          B = 1000,
                                                          show_progress = FALSE,
                                                          verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_hvmf_simple_calibration_scenarios(),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "hvmf_simple_ks_cvm_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_logistic_gaussian_calibration_study <- function(output_dir = NULL,
                                                          n_cores_outer = 1,
                                                          seed = 123,
                                                          show_progress = FALSE,
                                                          verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_simple_calibration_scenarios(),
    n_values = 30,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "logistic_gaussian_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_logistic_gaussian_calibration_study <- function(output_dir = NULL,
                                                         n_cores_outer = 1,
                                                         seed = 123,
                                                         M_outer = 500,
                                                         B = 500,
                                                         show_progress = FALSE,
                                                         verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_simple_calibration_scenarios(),
    n_values = c(50, 100),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "logistic_gaussian_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_logistic_gaussian_composite_calibration_study <- function(output_dir = NULL,
                                                                    n_cores_outer = 1,
                                                                    seed = 123,
                                                                    show_progress = FALSE,
                                                                    verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_composite_calibration_scenarios(),
    n_values = 24,
    M_outer = 6,
    B = 11,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "logistic_gaussian_composite_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_logistic_gaussian_composite_calibration_study <- function(output_dir = NULL,
                                                                   n_cores_outer = 1,
                                                                   seed = 123,
                                                                   M_outer = 100,
                                                                   B = 100,
                                                                   show_progress = FALSE,
                                                                   verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_logistic_gaussian_composite_calibration_scenarios(),
    n_values = c(30, 50),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "logistic_gaussian_composite_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_small_circle_simple_calibration_study <- function(output_dir = NULL,
                                                            n_cores_outer = 1,
                                                            seed = 123,
                                                            kappa_values = c(5, 20),
                                                            nu_values = c(0, 0.3, 0.7),
                                                            show_progress = FALSE,
                                                            verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_simple_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_simple_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_small_circle_simple_calibration_study <- function(output_dir = NULL,
                                                           n_cores_outer = 1,
                                                           seed = 123,
                                                           M_outer = 500,
                                                           B = 500,
                                                           kappa_values = c(5, 20),
                                                           nu_values = c(0, 0.3, 0.7),
                                                           show_progress = FALSE,
                                                           verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_simple_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_simple_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_small_circle_composite_calibration_study <- function(output_dir = NULL,
                                                               n_cores_outer = 1,
                                                               seed = 123,
                                                               kappa_values = c(5, 20),
                                                               nu_values = c(0, 0.3, 0.7),
                                                               show_progress = FALSE,
                                                               verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_composite_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_composite_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_small_circle_composite_calibration_study <- function(output_dir = NULL,
                                                              n_cores_outer = 1,
                                                              seed = 123,
                                                              M_outer = 500,
                                                              B = 500,
                                                              kappa_values = c(5, 20),
                                                              nu_values = c(0, 0.3, 0.7),
                                                              show_progress = FALSE,
                                                              verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_composite_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_composite_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_small_circle_symmetric_mixture2_simple_calibration_study <- function(output_dir = NULL,
                                                                                n_cores_outer = 1,
                                                                                seed = 123,
                                                                                kappa_values = c(12, 20),
                                                                                nu_values = c(0.35, 0.6),
                                                                                show_progress = FALSE,
                                                                                verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_simple_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_symmetric_mixture2_simple_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_small_circle_symmetric_mixture2_simple_calibration_study <- function(output_dir = NULL,
                                                                               n_cores_outer = 1,
                                                                               seed = 123,
                                                                               M_outer = 500,
                                                                               B = 500,
                                                                               kappa_values = c(12, 20),
                                                                               nu_values = c(0.35, 0.6),
                                                                               show_progress = FALSE,
                                                                               verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_simple_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_symmetric_mixture2_simple_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_smoke_small_circle_symmetric_mixture2_composite_calibration_study <- function(output_dir = NULL,
                                                                                   n_cores_outer = 1,
                                                                                   seed = 123,
                                                                                   kappa_values = c(12, 20),
                                                                                   nu_values = c(0.35, 0.6),
                                                                                   show_progress = FALSE,
                                                                                   verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_composite_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = 50,
    M_outer = 10,
    B = 19,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_symmetric_mixture2_composite_smoke"),
    show_progress = show_progress,
    verbose = verbose
  )
}

run_full_small_circle_symmetric_mixture2_composite_calibration_study <- function(output_dir = NULL,
                                                                                  n_cores_outer = 1,
                                                                                  seed = 123,
                                                                                  M_outer = 500,
                                                                                  B = 500,
                                                                                  kappa_values = c(12, 20),
                                                                                  nu_values = c(0.35, 0.6),
                                                                                  show_progress = FALSE,
                                                                                  verbose = TRUE) {
  run_bootstrap_calibration_study(
    scenarios = default_small_circle_symmetric_mixture2_composite_calibration_scenarios(
      kappa_values = kappa_values,
      nu_values = nu_values
    ),
    n_values = c(50, 100, 200),
    M_outer = M_outer,
    B = B,
    alpha_nominal = 0.05,
    alphas = c(0.01, 0.05, 0.10),
    statistics = c("ks", "cvm"),
    n_cores_outer = n_cores_outer,
    seed = seed,
    output_dir = output_dir %||% file.path("output", "calibration", "bootstrap", "small_circle_symmetric_mixture2_composite_full"),
    show_progress = show_progress,
    verbose = verbose
  )
}
