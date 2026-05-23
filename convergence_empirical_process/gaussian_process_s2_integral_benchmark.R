# ============================================================================
# S2 EXACT-INTEGRAL BENCHMARKS FOR THE LIMIT GAUSSIAN PROCESS
# These functions generate validation artifacts. When saved to disk, they
# should go under tests/benchmark_outputs/, not under routine output/ folders.
# ============================================================================

utils_path_s2 <- if (file.exists("utils.R")) {
  "utils.R"
} else if (file.exists(file.path("..", "utils.R"))) {
  file.path("..", "utils.R")
} else {
  stop("Could not find utils.R in current directory or parent directory.")
}
source(utils_path_s2)

vmf_gp_path_s2 <- if (file.exists(file.path("convergence_empirical_process", "gaussian_process_vmf.R"))) {
  file.path("convergence_empirical_process", "gaussian_process_vmf.R")
} else if (file.exists("gaussian_process_vmf.R")) {
  "gaussian_process_vmf.R"
} else {
  stop("Could not find gaussian_process_vmf.R in the expected locations.")
}
source(vmf_gp_path_s2)


# ============================================================================
# INTERNAL HELPERS
# ============================================================================

s2_point_from_angles <- function(azimuth, colatitude) {
  c(
    sin(colatitude) * cos(azimuth),
    sin(colatitude) * sin(azimuth),
    cos(colatitude)
  )
}

distance_from_geodesic_radius <- function(radius, distance_type = "geodesic") {
  distance_type <- match.arg(distance_type, choices = c("chordal", "geodesic"))
  if (distance_type == "geodesic") {
    return(radius)
  }
  2 * sin(radius / 2)
}

default_s2_benchmark_cases <- function(distance_type = "geodesic") {
  make_case <- function(name, omega1, omega2, radius1, radius2) {
    list(
      name = name,
      omega1 = omega1 / sqrt(sum(omega1^2)),
      omega2 = omega2 / sqrt(sum(omega2^2)),
      t1 = distance_from_geodesic_radius(radius1, distance_type),
      t2 = distance_from_geodesic_radius(radius2, distance_type),
      geodesic_radius1 = radius1,
      geodesic_radius2 = radius2
    )
  }

  list(
    make_case(
      name = "same_center",
      omega1 = c(0, 0, 1),
      omega2 = c(0, 0, 1),
      radius1 = 0.60,
      radius2 = 0.95
    ),
    make_case(
      name = "moderate_overlap",
      omega1 = c(0, 0, 1),
      omega2 = s2_point_from_angles(azimuth = 0.40, colatitude = 0.90),
      radius1 = 0.85,
      radius2 = 1.05
    ),
    make_case(
      name = "oblique",
      omega1 = s2_point_from_angles(azimuth = 0.25, colatitude = 0.85),
      omega2 = s2_point_from_angles(azimuth = 2.10, colatitude = 1.00),
      radius1 = 0.90,
      radius2 = 0.95
    ),
    make_case(
      name = "near_tangent",
      omega1 = s2_point_from_angles(azimuth = 0.55, colatitude = 0.95),
      omega2 = s2_point_from_angles(azimuth = 1.55, colatitude = 1.00),
      radius1 = 0.48,
      radius2 = 0.50
    )
  )
}

default_s2_t_grid <- function(t_points, distance_type = "geodesic") {
  radii <- seq(0.35, 2.35, length.out = t_points)
  distance_from_geodesic_radius(radii, distance_type = distance_type)
}

default_exact_settings_s2 <- function() {
  list(
    list(
      name = "integral_fast",
      rel.tol_outer = 1e-5,
      abs.tol_outer = 1e-7,
      rel.tol_inner = 1e-6,
      abs.tol_inner = 1e-8,
      subdivisions_outer = 150L,
      subdivisions_inner = 150L
    ),
    list(
      name = "integral_balanced",
      rel.tol_outer = 1e-6,
      abs.tol_outer = 1e-8,
      rel.tol_inner = 1e-7,
      abs.tol_inner = 1e-9,
      subdivisions_outer = 200L,
      subdivisions_inner = 200L
    )
  )
}

default_reference_settings_s2 <- function() {
  list(
    rel.tol_outer = 1e-8,
    abs.tol_outer = 1e-10,
    rel.tol_inner = 1e-8,
    abs.tol_inner = 1e-10,
    subdivisions_outer = 300L,
    subdivisions_inner = 300L
  )
}

vmf_s2_mu_projection_tail_closed_form <- function(a, kappa, tol = 1e-12) {
  a <- pmin(pmax(a, -1), 1)
  if (kappa <= tol) {
    return((1 - a) / 2)
  }
  (exp(kappa) - exp(kappa * a)) / (exp(kappa) - exp(-kappa))
}

vmf_s2_mu_projection_interval_closed_form <- function(lower, upper, kappa, tol = 1e-12) {
  lower <- pmin(pmax(lower, -1), 1)
  upper <- pmin(pmax(upper, -1), 1)
  if (upper < lower) {
    return(0)
  }
  if (kappa <= tol) {
    return((upper - lower) / 2)
  }
  (exp(kappa * upper) - exp(kappa * lower)) / (exp(kappa) - exp(-kappa))
}

default_s2_closed_form_cases <- function(mu = c(0, 0, 1), distance_type = "geodesic") {
  mu <- as.numeric(mu)
  mu <- mu / sqrt(sum(mu^2))
  minus_mu <- -mu
  make_case <- function(name, omega1, omega2, radius1, radius2, truth_type) {
    list(
      name = name,
      omega1 = omega1,
      omega2 = omega2,
      t1 = distance_from_geodesic_radius(radius1, distance_type),
      t2 = distance_from_geodesic_radius(radius2, distance_type),
      geodesic_radius1 = radius1,
      geodesic_radius2 = radius2,
      truth_type = truth_type
    )
  }

  list(
    make_case("mu_mu_nested", mu, mu, 0.60, 0.95, "mu_mu"),
    make_case("mu_mu_reversed", mu, mu, 1.20, 0.70, "mu_mu"),
    make_case("minusmu_minusmu", minus_mu, minus_mu, 0.65, 1.00, "minusmu_minusmu"),
    make_case("mu_minusmu_overlap", mu, minus_mu, 0.85, 0.90, "mu_minusmu"),
    make_case("mu_minusmu_disjoint", mu, minus_mu, 0.40, 0.55, "mu_minusmu")
  )
}

closed_form_joint_probability_vmf_s2_special <- function(case,
                                                          kappa,
                                                          distance_type = "geodesic") {
  a1 <- sphere_distance_to_dot_threshold(case$t1, distance_type = distance_type)
  a2 <- sphere_distance_to_dot_threshold(case$t2, distance_type = distance_type)

  switch(
    case$truth_type,
    mu_mu = {
      vmf_s2_mu_projection_tail_closed_form(max(a1, a2), kappa = kappa)
    },
    minusmu_minusmu = {
      vmf_s2_mu_projection_interval_closed_form(
        lower = -1,
        upper = min(-a1, -a2),
        kappa = kappa
      )
    },
    mu_minusmu = {
      vmf_s2_mu_projection_interval_closed_form(
        lower = a1,
        upper = -a2,
        kappa = kappa
      )
    },
    stop("Unknown closed-form case type.")
  )
}


# ============================================================================
# JOINT-PROBABILITY BENCHMARK
# ============================================================================

#' Benchmark the exact S^2 integral formula against Monte Carlo
#' @param mu Mean direction on S^2
#' @param kappas Concentration parameters to benchmark
#' @param distance_type Either "chordal" or "geodesic"
#' @param cases Optional list of benchmark cases
#' @param exact_settings Optional list of non-reference exact settings
#' @param reference_settings Optional list of high-accuracy reference settings
#' @param mc_sizes Monte Carlo sample sizes
#' @param seed Optional seed
#' @param output_dir Optional output directory for CSV/plots
#' @param save_plots Whether to save plots when output_dir is provided
#' @return List with metrics and optional plots
benchmark_joint_probability_vmf_s2 <- function(mu = c(0, 0, 1),
                                               kappas = c(0.5, 2, 5),
                                               distance_type = c("geodesic", "chordal"),
                                               cases = NULL,
                                               exact_settings = default_exact_settings_s2(),
                                               reference_settings = default_reference_settings_s2(),
                                               mc_sizes = c(1000, 5000, 10000, 50000),
                                               seed = 123,
                                               output_dir = NULL,
                                               save_plots = FALSE) {
  distance_type <- match.arg(distance_type)
  mu <- as.numeric(mu)
  mu <- mu / sqrt(sum(mu^2))
  if (is.null(cases)) {
    cases <- default_s2_benchmark_cases(distance_type = distance_type)
  }

  if (!is.null(seed)) set.seed(seed)
  rows <- list()
  row_idx <- 0

  for (kappa in kappas) {
    for (case in cases) {
      ref_start <- Sys.time()
      ref_prob <- joint_probability_vmf_s2_simple_integral(
        omega1 = case$omega1,
        t1 = case$t1,
        omega2 = case$omega2,
        t2 = case$t2,
        mu = mu,
        kappa = kappa,
        distance_type = distance_type,
        rel.tol_outer = reference_settings$rel.tol_outer,
        abs.tol_outer = reference_settings$abs.tol_outer,
        rel.tol_inner = reference_settings$rel.tol_inner,
        abs.tol_inner = reference_settings$abs.tol_inner,
        subdivisions_outer = reference_settings$subdivisions_outer,
        subdivisions_inner = reference_settings$subdivisions_inner
      )
      ref_time <- as.numeric(difftime(Sys.time(), ref_start, units = "secs"))

      row_idx <- row_idx + 1
      rows[[row_idx]] <- data.frame(
        distance_type = distance_type,
        kappa = kappa,
        case = case$name,
        method = "integral_reference",
        tuning = NA_character_,
        estimate = ref_prob,
        abs_error = 0,
        rel_error = 0,
        time_sec = ref_time,
        stringsAsFactors = FALSE
      )

      for (setting in exact_settings) {
        exact_start <- Sys.time()
        exact_prob <- joint_probability_vmf_s2_simple_integral(
          omega1 = case$omega1,
          t1 = case$t1,
          omega2 = case$omega2,
          t2 = case$t2,
          mu = mu,
          kappa = kappa,
          distance_type = distance_type,
          rel.tol_outer = setting$rel.tol_outer,
          abs.tol_outer = setting$abs.tol_outer,
          rel.tol_inner = setting$rel.tol_inner,
          abs.tol_inner = setting$abs.tol_inner,
          subdivisions_outer = setting$subdivisions_outer,
          subdivisions_inner = setting$subdivisions_inner
        )
        exact_time <- as.numeric(difftime(Sys.time(), exact_start, units = "secs"))

        row_idx <- row_idx + 1
        rows[[row_idx]] <- data.frame(
          distance_type = distance_type,
          kappa = kappa,
          case = case$name,
          method = "integral",
          tuning = setting$name,
          estimate = exact_prob,
          abs_error = abs(exact_prob - ref_prob),
          rel_error = abs(exact_prob - ref_prob) / ifelse(ref_prob > 0, ref_prob, 1),
          time_sec = exact_time,
          stringsAsFactors = FALSE
        )
      }

      for (mc_n in mc_sizes) {
        if (!is.null(seed)) set.seed(seed)
        mc_start <- Sys.time()
        sample_points <- rotasym::r_vMF(n = mc_n, mu = mu, kappa = kappa)
        dots1 <- as.numeric(sample_points %*% case$omega1)
        dots2 <- as.numeric(sample_points %*% case$omega2)
        a1 <- sphere_distance_to_dot_threshold(case$t1, distance_type = distance_type)
        a2 <- sphere_distance_to_dot_threshold(case$t2, distance_type = distance_type)
        mc_prob <- mean(dots1 >= a1 & dots2 >= a2)
        mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))

        row_idx <- row_idx + 1
        rows[[row_idx]] <- data.frame(
          distance_type = distance_type,
          kappa = kappa,
          case = case$name,
          method = "mc",
          tuning = paste0("n_mc=", mc_n),
          estimate = mc_prob,
          abs_error = abs(mc_prob - ref_prob),
          rel_error = abs(mc_prob - ref_prob) / ifelse(ref_prob > 0, ref_prob, 1),
          time_sec = mc_time,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  metrics <- do.call(rbind, rows)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparable <- subset(metrics, method != "integral_reference")
    error_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = abs_error, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Absolute error vs reference integral",
        title = paste("S2 Joint-Probability Benchmark (", distance_type, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    time_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = time_sec, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Time (seconds)",
        title = paste("S2 Joint-Probability Timing (", distance_type, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    csv_name <- paste0("s2_joint_probability_benchmark_", distance_type, ".csv")
    utils::write.csv(metrics, file.path(output_dir, csv_name), row.names = FALSE)
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_joint_probability_error_", distance_type, ".png")),
        error_plot,
        width = 10,
        height = 7,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_joint_probability_time_", distance_type, ".png")),
        time_plot,
        width = 10,
        height = 7,
        dpi = 300
      )
    }
  }

  list(
    metrics = metrics,
    error_plot = error_plot,
    time_plot = time_plot
  )
}

#' Benchmark S^2 special cases with truly closed-form joint probabilities
#' @param mu Mean direction on S^2
#' @param kappas Concentration parameters to benchmark
#' @param distance_type Either "chordal" or "geodesic"
#' @param cases Optional list of special cases with closed-form truth
#' @param exact_settings Optional list of exact integral settings
#' @param mc_sizes Monte Carlo sample sizes
#' @param seed Optional seed
#' @param output_dir Optional output directory for CSV/plots
#' @param save_plots Whether to save plots when output_dir is provided
#' @return List with metrics and optional plots
benchmark_joint_probability_vmf_s2_closed_form <- function(mu = c(0, 0, 1),
                                                           kappas = c(0.5, 2, 5),
                                                           distance_type = c("geodesic", "chordal"),
                                                           cases = NULL,
                                                           exact_settings = default_exact_settings_s2(),
                                                           mc_sizes = c(1000, 5000, 10000, 50000),
                                                           seed = 123,
                                                           output_dir = NULL,
                                                           save_plots = FALSE) {
  distance_type <- match.arg(distance_type)
  mu <- as.numeric(mu)
  mu <- mu / sqrt(sum(mu^2))
  if (is.null(cases)) {
    cases <- default_s2_closed_form_cases(mu = mu, distance_type = distance_type)
  }

  rows <- list()
  row_idx <- 0

  for (kappa in kappas) {
    for (case in cases) {
      truth_prob <- closed_form_joint_probability_vmf_s2_special(
        case = case,
        kappa = kappa,
        distance_type = distance_type
      )

      row_idx <- row_idx + 1
      rows[[row_idx]] <- data.frame(
        distance_type = distance_type,
        kappa = kappa,
        case = case$name,
        method = "closed_form_truth",
        tuning = NA_character_,
        estimate = truth_prob,
        abs_error = 0,
        rel_error = 0,
        time_sec = 0,
        stringsAsFactors = FALSE
      )

      for (setting in exact_settings) {
        exact_start <- Sys.time()
        exact_prob <- joint_probability_vmf_s2_simple_integral(
          omega1 = case$omega1,
          t1 = case$t1,
          omega2 = case$omega2,
          t2 = case$t2,
          mu = mu,
          kappa = kappa,
          distance_type = distance_type,
          rel.tol_outer = setting$rel.tol_outer,
          abs.tol_outer = setting$abs.tol_outer,
          rel.tol_inner = setting$rel.tol_inner,
          abs.tol_inner = setting$abs.tol_inner,
          subdivisions_outer = setting$subdivisions_outer,
          subdivisions_inner = setting$subdivisions_inner
        )
        exact_time <- as.numeric(difftime(Sys.time(), exact_start, units = "secs"))

        row_idx <- row_idx + 1
        rows[[row_idx]] <- data.frame(
          distance_type = distance_type,
          kappa = kappa,
          case = case$name,
          method = "integral",
          tuning = setting$name,
          estimate = exact_prob,
          abs_error = abs(exact_prob - truth_prob),
          rel_error = abs(exact_prob - truth_prob) / ifelse(truth_prob > 0, truth_prob, 1),
          time_sec = exact_time,
          stringsAsFactors = FALSE
        )
      }

      for (mc_n in mc_sizes) {
        if (!is.null(seed)) set.seed(seed)
        mc_start <- Sys.time()
        sample_points <- rotasym::r_vMF(n = mc_n, mu = mu, kappa = kappa)
        dots1 <- as.numeric(sample_points %*% case$omega1)
        dots2 <- as.numeric(sample_points %*% case$omega2)
        a1 <- sphere_distance_to_dot_threshold(case$t1, distance_type = distance_type)
        a2 <- sphere_distance_to_dot_threshold(case$t2, distance_type = distance_type)
        mc_prob <- mean(dots1 >= a1 & dots2 >= a2)
        mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))

        row_idx <- row_idx + 1
        rows[[row_idx]] <- data.frame(
          distance_type = distance_type,
          kappa = kappa,
          case = case$name,
          method = "mc",
          tuning = paste0("n_mc=", mc_n),
          estimate = mc_prob,
          abs_error = abs(mc_prob - truth_prob),
          rel_error = abs(mc_prob - truth_prob) / ifelse(truth_prob > 0, truth_prob, 1),
          time_sec = mc_time,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  metrics <- do.call(rbind, rows)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparable <- subset(metrics, method != "closed_form_truth")
    error_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = abs_error, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Absolute error vs closed-form truth",
        title = paste("S2 Closed-Form Special Cases (", distance_type, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    time_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = time_sec, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_grid(case ~ kappa, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Time (seconds)",
        title = paste("S2 Closed-Form Timing (", distance_type, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    csv_name <- paste0("s2_joint_probability_closed_form_benchmark_", distance_type, ".csv")
    utils::write.csv(metrics, file.path(output_dir, csv_name), row.names = FALSE)
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_joint_probability_closed_form_error_", distance_type, ".png")),
        error_plot,
        width = 10,
        height = 8,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_joint_probability_closed_form_time_", distance_type, ".png")),
        time_plot,
        width = 10,
        height = 8,
        dpi = 300
      )
    }
  }

  list(
    metrics = metrics,
    error_plot = error_plot,
    time_plot = time_plot
  )
}


# ============================================================================
# MATRIX BENCHMARK
# ============================================================================

#' Benchmark the exact S^2 covariance matrix against Monte Carlo
#' @param mu Mean direction on S^2
#' @param kappa Concentration parameter
#' @param distance_type Either "chordal" or "geodesic"
#' @param scenarios Optional list of grid scenarios
#' @param exact_settings Optional list of non-reference exact settings
#' @param reference_settings Optional list of high-accuracy reference settings
#' @param mc_sizes Monte Carlo covariance sample sizes
#' @param n_cores Number of cores for the MC covariance route
#' @param seed Optional seed
#' @param output_dir Optional output directory for CSV/plots
#' @param save_plots Whether to save plots when output_dir is provided
#' @return List with metrics and optional plots
benchmark_covariance_vmf_s2 <- function(mu = c(0, 0, 1),
                                        kappa = 0.5,
                                        distance_type = c("geodesic", "chordal"),
                                        scenarios = list(
                                          list(name = "small", n_omega = 6, t_points = 4),
                                          list(name = "medium", n_omega = 8, t_points = 6)
                                        ),
                                        exact_settings = default_exact_settings_s2(),
                                        reference_settings = default_reference_settings_s2(),
                                        mc_sizes = c(10000, 50000, 100000),
                                        n_cores = 1,
                                        seed = 123,
                                        output_dir = NULL,
                                        save_plots = FALSE) {
  distance_type <- match.arg(distance_type)
  mu <- as.numeric(mu)
  mu <- mu / sqrt(sum(mu^2))

  rows <- list()
  row_idx <- 0

  for (scenario in scenarios) {
    omega_grid <- generate_canonical_lattice(scenario$n_omega, dim = 3)
    t_grid <- default_s2_t_grid(scenario$t_points, distance_type = distance_type)

    ref_start <- Sys.time()
    sigma_ref <- cov_vmf_s2_simple_integral(
      omega_grid = omega_grid,
      t_grid = t_grid,
      mu = mu,
      kappa = kappa,
      distance_type = distance_type,
      rel.tol_outer = reference_settings$rel.tol_outer,
      abs.tol_outer = reference_settings$abs.tol_outer,
      rel.tol_inner = reference_settings$rel.tol_inner,
      abs.tol_inner = reference_settings$abs.tol_inner,
      subdivisions_outer = reference_settings$subdivisions_outer,
      subdivisions_inner = reference_settings$subdivisions_inner
    )
    ref_time <- as.numeric(difftime(Sys.time(), ref_start, units = "secs"))
    ref_diag <- validate_covariance_matrix(
      sigma_ref,
      symmetry_tol = 1e-10,
      psd_tol = 1e-10,
      stop_on_failure = FALSE
    )
    fro_ref <- sqrt(sum(sigma_ref^2))

    row_idx <- row_idx + 1
    rows[[row_idx]] <- data.frame(
      distance_type = distance_type,
      scenario = scenario$name,
      n_omega = scenario$n_omega,
      t_points = scenario$t_points,
      method = "integral_reference",
      tuning = NA_character_,
      time_sec = ref_time,
      max_abs_diff = 0,
      rel_fro_error = 0,
      symmetry_gap = ref_diag$symmetry_gap,
      min_eigenvalue = ref_diag$min_eigenvalue,
      stringsAsFactors = FALSE
    )

    for (setting in exact_settings) {
      exact_start <- Sys.time()
      sigma_exact <- cov_vmf_s2_simple_integral(
        omega_grid = omega_grid,
        t_grid = t_grid,
        mu = mu,
        kappa = kappa,
        distance_type = distance_type,
        rel.tol_outer = setting$rel.tol_outer,
        abs.tol_outer = setting$abs.tol_outer,
        rel.tol_inner = setting$rel.tol_inner,
        abs.tol_inner = setting$abs.tol_inner,
        subdivisions_outer = setting$subdivisions_outer,
        subdivisions_inner = setting$subdivisions_inner
      )
      exact_time <- as.numeric(difftime(Sys.time(), exact_start, units = "secs"))
      exact_diag <- validate_covariance_matrix(
        sigma_exact,
        symmetry_tol = 1e-10,
        psd_tol = 1e-10,
        stop_on_failure = FALSE
      )
      diff_mat <- sigma_exact - sigma_ref

      row_idx <- row_idx + 1
      rows[[row_idx]] <- data.frame(
        distance_type = distance_type,
        scenario = scenario$name,
        n_omega = scenario$n_omega,
        t_points = scenario$t_points,
        method = "integral",
        tuning = setting$name,
        time_sec = exact_time,
        max_abs_diff = max(abs(diff_mat)),
        rel_fro_error = sqrt(sum(diff_mat^2)) / ifelse(fro_ref > 0, fro_ref, 1),
        symmetry_gap = exact_diag$symmetry_gap,
        min_eigenvalue = exact_diag$min_eigenvalue,
        stringsAsFactors = FALSE
      )
    }

    for (mc_n in mc_sizes) {
      mc_start <- Sys.time()
      sigma_mc <- cov_vmf(
        omega_grid = omega_grid,
        t_grid = t_grid,
        mu = mu,
        kappa = kappa,
        distance_type = distance_type,
        n_mc_samples = mc_n,
        n_cores = n_cores,
        seed = seed,
        h0 = "simple",
        cov_method = "mc"
      )
      mc_time <- as.numeric(difftime(Sys.time(), mc_start, units = "secs"))
      mc_diag <- validate_covariance_matrix(
        sigma_mc,
        symmetry_tol = 1e-10,
        psd_tol = 1e-10,
        stop_on_failure = FALSE
      )
      diff_mat <- sigma_mc - sigma_ref

      row_idx <- row_idx + 1
      rows[[row_idx]] <- data.frame(
        distance_type = distance_type,
        scenario = scenario$name,
        n_omega = scenario$n_omega,
        t_points = scenario$t_points,
        method = "mc",
        tuning = paste0("n_mc=", mc_n),
        time_sec = mc_time,
        max_abs_diff = max(abs(diff_mat)),
        rel_fro_error = sqrt(sum(diff_mat^2)) / ifelse(fro_ref > 0, fro_ref, 1),
        symmetry_gap = mc_diag$symmetry_gap,
        min_eigenvalue = mc_diag$min_eigenvalue,
        stringsAsFactors = FALSE
      )
    }
  }

  metrics <- do.call(rbind, rows)
  error_plot <- NULL
  time_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparable <- subset(metrics, method != "integral_reference")
    error_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = rel_fro_error, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~scenario, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Relative Frobenius error vs reference integral",
        title = paste("S2 Covariance Benchmark (", distance_type, ", kappa=", kappa, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    time_plot <- ggplot2::ggplot(
      comparable,
      ggplot2::aes(x = tuning, y = time_sec, color = method, group = method)
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~scenario, scales = "free_x") +
      ggplot2::labs(
        x = "Tuning parameter",
        y = "Time (seconds)",
        title = paste("S2 Covariance Timing (", distance_type, ", kappa=", kappa, ")", sep = "")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    csv_name <- paste0("s2_covariance_benchmark_", distance_type, "_kappa", gsub("\\.", "p", sprintf("%.2f", kappa)), ".csv")
    utils::write.csv(metrics, file.path(output_dir, csv_name), row.names = FALSE)
    if (isTRUE(save_plots) && !is.null(error_plot) && !is.null(time_plot)) {
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_covariance_error_", distance_type, "_kappa", gsub("\\.", "p", sprintf("%.2f", kappa)), ".png")),
        error_plot,
        width = 9,
        height = 5,
        dpi = 300
      )
      ggplot2::ggsave(
        file.path(output_dir, paste0("s2_covariance_time_", distance_type, "_kappa", gsub("\\.", "p", sprintf("%.2f", kappa)), ".png")),
        time_plot,
        width = 9,
        height = 5,
        dpi = 300
      )
    }
  }

  list(
    metrics = metrics,
    error_plot = error_plot,
    time_plot = time_plot
  )
}

cat("S2 exact-integral benchmark functions loaded successfully!\n")
