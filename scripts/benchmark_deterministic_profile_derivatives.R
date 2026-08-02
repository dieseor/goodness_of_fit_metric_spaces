args <- commandArgs(trailingOnly = TRUE)
cores_arg <- grep("^--cores=", args, value = TRUE)
cores <- if (length(cores_arg)) {
  as.integer(sub("^--cores=", "", cores_arg[[1L]]))
} else {
  8L
}
cores <- max(1L, min(8L, cores))

source(file.path("bootstrap", "multiplier_bootstrap.R"))

bench_unit <- function(x) x / sqrt(sum(x^2))
bench_hpoint <- function(rho, direction) {
  direction <- bench_unit(direction)
  c(cosh(rho), sinh(rho) * direction)
}

bench_vmf_F <- function(xi, omega, t) {
  threshold <- cos(t)
  if (threshold <= -1) return(1)
  if (threshold >= 1) return(0)
  integrate(
    function(s) vmf_projected_density_canonical(s, xi, omega)$density,
    threshold,
    1,
    rel.tol = 2e-11,
    abs.tol = 2e-13,
    subdivisions = 3000L
  )$value
}

bench_hvmf_F <- function(xi, omega, t) {
  kappa <- sqrt(-hvmf_minkowski_inner_product(xi, xi))
  a <- hvmf_minkowski_inner_product(xi, omega)
  b_sq <- profile_derivative_nonnegative_square(
    a^2 - kappa^2,
    max(a^2, kappa^2),
    "benchmark HvMF b^2"
  )
  chi <- asinh(sqrt(b_sq) / kappa)
  integrate(
    function(r) hvmf_radial_density(
      r,
      q = length(xi) - 1L,
      kappa = kappa,
      chi = chi
    ),
    0,
    t,
    rel.tol = 5e-10,
    abs.tol = 2e-12,
    subdivisions = 3000L
  )$value
}

bench_fd_matrix <- function(F_function, xi, omega, t_values, relative_step) {
  vapply(seq_along(xi), function(j) {
    step <- relative_step * max(1, abs(xi[[j]]))
    plus <- minus <- xi
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    (
      vapply(t_values, function(t) F_function(plus, omega, t), numeric(1)) -
        vapply(t_values, function(t) F_function(minus, omega, t), numeric(1))
    ) / (2 * step)
  }, numeric(length(t_values)))
}

bench_error_metrics <- function(estimate, reference) {
  error <- estimate - reference
  relative_index <- abs(reference) >= 1e-3
  c(
    max_abs_error = max(abs(error)),
    rmse = sqrt(mean(error^2)),
    max_relative_error = if (any(relative_index)) {
      max(abs(error[relative_index] / reference[relative_index]))
    } else {
      NA_real_
    }
  )
}

bench_vmf_mc <- function(scenario, M, seed) {
  set.seed(seed)
  sample <- rotasym::r_vMF(M, scenario$mu, scenario$kappa)
  score <- sweep(sample, 2L, A_q(scenario$kappa, scenario$q) * scenario$mu, "-")
  distances <- acos(pmin(pmax(drop(sample %*% scenario$omega), -1), 1))
  indicators <- outer(distances, scenario$t, "<=") * 1
  crossprod(indicators, score) / M
}

bench_hvmf_mc <- function(scenario, M, seed) {
  set.seed(seed)
  sample <- if (scenario$q == 2L) {
    rhvmf_h2_polar(M, scenario$mu, scenario$kappa)
  } else {
    rhvmf_polar(M, scenario$mu, scenario$kappa)
  }
  score <- hvmf_canonical_score_matrix(sample, scenario$xi)
  distances <- drop(hvmf_distance_matrix_hq(
    matrix(scenario$omega, nrow = 1L),
    sample
  ))
  indicators <- outer(distances, scenario$t, "<=") * 1
  crossprod(indicators, score) / M
}

scenarios <- list(
  list(
    id = "vmf_paper_s2_k2",
    family = "vmf",
    q = 2L,
    kappa = 2,
    mu = c(1, 0, 0),
    omega = bench_unit(c(0.3, 0.8, -0.2)),
    t = seq(0.25, 2.75, length.out = 9L),
    fd_step = 1e-5
  ),
  list(
    id = "vmf_section6_s9_k10",
    family = "vmf",
    q = 9L,
    kappa = 10,
    mu = c(1, rep.int(0, 9L)),
    omega = bench_unit(c(0.6, 0.8, rep.int(0, 8L))),
    t = seq(0.25, 2.2, length.out = 9L),
    fd_step = 8e-6
  ),
  list(
    id = "hvmf_calibration_h2_k50",
    family = "hvmf",
    q = 2L,
    kappa = 50,
    mu = bench_hpoint(0.7, c(cos(0.4), sin(0.4))),
    omega = bench_hpoint(0.9, c(cos(-0.2), sin(-0.2))),
    t = seq(0.15, 1.2, length.out = 9L),
    fd_step = 3e-6
  ),
  list(
    id = "hvmf_paper_h2_k200",
    family = "hvmf",
    q = 2L,
    kappa = 200,
    mu = c(sqrt(2), 1 / sqrt(2), 1 / sqrt(2)),
    omega = bench_hpoint(0.8, c(-0.8, 0.6)),
    t = seq(0.35, 1.7, length.out = 9L),
    fd_step = 2e-6
  ),
  list(
    id = "hvmf_section6_h10_k10",
    family = "hvmf",
    q = 10L,
    kappa = 10,
    mu = c(sqrt(2), 1, rep.int(0, 9L)),
    omega = c(sqrt(2), 0, 1, rep.int(0, 8L)),
    t = seq(0.45, 2.0, length.out = 9L),
    fd_step = 4e-6
  )
)
scenarios <- lapply(scenarios, function(scenario) {
  scenario$xi <- scenario$kappa * scenario$mu
  scenario
})

precision_one <- function(scenario) {
  reference_coarse <- if (identical(scenario$family, "vmf")) {
    bench_fd_matrix(
      bench_vmf_F,
      scenario$xi,
      scenario$omega,
      scenario$t,
      scenario$fd_step
    )
  } else {
    bench_fd_matrix(
      bench_hvmf_F,
      scenario$xi,
      scenario$omega,
      scenario$t,
      scenario$fd_step
    )
  }
  reference <- reference_coarse
  deterministic <- if (identical(scenario$family, "vmf")) {
    vmf_profile_and_derivative_xi(
      scenario$omega,
      scenario$xi,
      scenario$t,
      "geodesic",
      4097L
    )$derivative
  } else {
    hvmf_profile_and_derivative_xi(
      scenario$omega,
      scenario$xi,
      scenario$t,
      4097L
    )$derivative
  }
  deterministic_metrics <- bench_error_metrics(deterministic, reference)
  rows <- data.frame(
    scenario = scenario$id,
    method = "quadrature",
    M = NA_integer_,
    seeds = 1L,
    max_abs_error = deterministic_metrics[["max_abs_error"]],
    rmse = deterministic_metrics[["rmse"]],
    max_relative_error = deterministic_metrics[["max_relative_error"]],
    max_abs_bias = deterministic_metrics[["max_abs_error"]],
    mean_mc_sd = 0,
    max_mc_sd = 0,
    stringsAsFactors = FALSE
  )

  for (M in c(250L, 1000L, 10000L)) {
    estimates <- lapply(1:5, function(seed_offset) {
      if (identical(scenario$family, "vmf")) {
        bench_vmf_mc(scenario, M, 9000L + 100L * M + seed_offset)
      } else {
        bench_hvmf_mc(scenario, M, 9000L + 100L * M + seed_offset)
      }
    })
    estimate_array <- simplify2array(estimates)
    mean_estimate <- apply(estimate_array, c(1L, 2L), mean)
    sd_estimate <- apply(estimate_array, c(1L, 2L), stats::sd)
    all_errors <- vapply(estimates, function(estimate) {
      as.numeric(estimate - reference)
    }, numeric(length(reference)))
    relative_index <- abs(reference) >= 1e-3
    rows <- rbind(rows, data.frame(
      scenario = scenario$id,
      method = "score_mc",
      M = M,
      seeds = length(estimates),
      max_abs_error = max(abs(all_errors)),
      rmse = sqrt(mean(all_errors^2)),
      max_relative_error = if (any(relative_index)) {
        max(abs(all_errors[relative_index, , drop = FALSE] /
          as.numeric(reference)[relative_index]))
      } else {
        NA_real_
      },
      max_abs_bias = max(abs(mean_estimate - reference)),
      mean_mc_sd = mean(sd_estimate),
      max_mc_sd = max(sd_estimate),
      stringsAsFactors = FALSE
    ))
  }
  rows
}

precision_list <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(scenarios, precision_one, mc.cores = min(cores, length(scenarios)))
} else {
  lapply(scenarios, precision_one)
}
precision_results <- do.call(rbind, precision_list)

bench_time <- function(expression, repetitions = 7L) {
  invisible(force(expression()))
  timings <- replicate(
    repetitions,
    system.time(invisible(expression()))[["elapsed"]]
  )
  c(median = median(timings), min = min(timings))
}

timing_one <- function(scenario,
                       center_count = 8L,
                       threshold_count = 9L,
                       repetitions = 7L) {
  t_values <- seq(
    min(scenario$t),
    max(scenario$t),
    length.out = threshold_count
  )
  if (identical(scenario$family, "vmf")) {
    center_base <- rbind(
      scenario$omega,
      scenario$mu,
      -scenario$mu,
      bench_unit(rev(scenario$omega))
    )
    centers <- center_base[rep(seq_len(nrow(center_base)), length.out = center_count), , drop = FALSE]
    deterministic_prepare <- function() {
      lapply(seq_len(nrow(centers)), function(i) {
        vmf_profile_derivative_table_xi(centers[i, ], scenario$xi, 4097L)
      })
    }
    tables <- deterministic_prepare()
    deterministic_eval <- function() {
      lapply(seq_len(nrow(centers)), function(i) {
        vmf_profile_and_derivative_xi(
          centers[i, ],
          scenario$xi,
          t_values,
          "geodesic",
          table = tables[[i]]
        )$derivative
      })
    }
  } else {
    spatial_directions <- diag(scenario$q)[seq_len(min(scenario$q, center_count)), , drop = FALSE]
    centers <- t(vapply(seq_len(center_count), function(i) {
      direction <- spatial_directions[(i - 1L) %% nrow(spatial_directions) + 1L, ]
      bench_hpoint(0.3 + 0.1 * ((i - 1L) %% 7L), direction)
    }, numeric(scenario$q + 1L)))
    deterministic_prepare <- function() {
      lapply(seq_len(nrow(centers)), function(i) {
        hvmf_profile_derivative_table_xi(
          centers[i, ],
          scenario$xi,
          max(t_values),
          4097L
        )
      })
    }
    tables <- deterministic_prepare()
    deterministic_eval <- function() {
      lapply(seq_len(nrow(centers)), function(i) {
        hvmf_profile_and_derivative_xi(
          centers[i, ],
          scenario$xi,
          t_values,
          table = tables[[i]]
        )$derivative
      })
    }
  }
  det_prepare_time <- bench_time(deterministic_prepare, repetitions)
  det_eval_time <- bench_time(deterministic_eval, repetitions)
  rows <- data.frame(
    scenario = scenario$id,
    method = "quadrature",
    M = NA_integer_,
    centers = nrow(centers),
    thresholds = length(t_values),
    prep_median_seconds = det_prepare_time[["median"]],
    prep_min_seconds = det_prepare_time[["min"]],
    eval_median_seconds = det_eval_time[["median"]],
    eval_min_seconds = det_eval_time[["min"]],
    stringsAsFactors = FALSE
  )

  for (M in c(1000L, 10000L)) {
    mc_prepare <- function() {
      if (identical(scenario$family, "vmf")) {
        sample <- rotasym::r_vMF(M, scenario$mu, scenario$kappa)
        score <- sweep(
          sample,
          2L,
          A_q(scenario$kappa, scenario$q) * scenario$mu,
          "-"
        )
        distances <- acos(pmin(pmax(sample %*% t(centers), -1), 1))
      } else {
        sample <- if (scenario$q == 2L) {
          rhvmf_h2_polar(M, scenario$mu, scenario$kappa)
        } else {
          rhvmf_polar(M, scenario$mu, scenario$kappa)
        }
        score <- hvmf_canonical_score_matrix(sample, scenario$xi)
        distances <- t(hvmf_distance_matrix_hq(centers, sample))
      }
      list(score = score, distances = distances)
    }
    mc_state <- mc_prepare()
    mc_eval <- function() {
      lapply(seq_len(ncol(mc_state$distances)), function(i) {
        indicators <- outer(mc_state$distances[, i], t_values, "<=") * 1
        crossprod(indicators, mc_state$score) / M
      })
    }
    prep_time <- bench_time(mc_prepare, repetitions)
    eval_time <- bench_time(mc_eval, repetitions)
    rows <- rbind(rows, data.frame(
      scenario = scenario$id,
      method = "score_mc",
      M = M,
      centers = nrow(centers),
      thresholds = length(t_values),
      prep_median_seconds = prep_time[["median"]],
      prep_min_seconds = prep_time[["min"]],
      eval_median_seconds = eval_time[["median"]],
      eval_min_seconds = eval_time[["min"]],
      stringsAsFactors = FALSE
    ))
  }
  rows
}

timing_jobs <- unlist(lapply(seq_along(scenarios), function(i) {
  base_jobs <- list(list(
    scenario_index = i,
    centers = 8L,
    thresholds = 9L,
    repetitions = 7L
  ), list(
    scenario_index = i,
    centers = 50L,
    thresholds = 50L,
    repetitions = 3L
  ))
  if (scenarios[[i]]$id %in% c("vmf_paper_s2_k2", "hvmf_paper_h2_k200")) {
    base_jobs <- c(base_jobs, list(list(
      scenario_index = i,
      centers = 200L,
      thresholds = 100L,
      repetitions = 3L
    )))
  }
  base_jobs
}), recursive = FALSE)
timing_results <- do.call(rbind, lapply(timing_jobs, function(job) {
  timing_one(
    scenario = scenarios[[job$scenario_index]],
    center_count = job$centers,
    threshold_count = job$thresholds,
    repetitions = job$repetitions
  )
}))

precision_path <- file.path(
  "benchmarks",
  "deterministic_profile_derivative_precision.csv"
)
timing_path <- file.path(
  "benchmarks",
  "deterministic_profile_derivative_timing.csv"
)
write.csv(precision_results, precision_path, row.names = FALSE)
write.csv(timing_results, timing_path, row.names = FALSE)

cat("\nPRECISION\n")
print(precision_results, row.names = FALSE)
cat("\nTIMING\n")
print(timing_results, row.names = FALSE)
cat(sprintf(
  "\nWrote %s and %s using %d core(s).\n",
  precision_path,
  timing_path,
  cores
))
