resolve_spline_benchmark_path <- function(...) {
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

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

source(resolve_spline_benchmark_path("bootstrap", "multiplier_bootstrap.R"))

benchmark_elapsed_seconds <- function(expr) {
  start <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - start
}

hyperboloid_point_benchmark <- function(chi, theta) {
  c(cosh(chi), sinh(chi) * cos(theta), sinh(chi) * sin(theta))
}

distance_profile_jp_grid_linear_reference <- function(omega_grid,
                                                      mu,
                                                      kappa,
                                                      psi,
                                                      t_grid,
                                                      distance_type = "geodesic",
                                                      n_u = 4097L,
                                                      n_delta = 1025L) {
  params <- jp_validate_parameters(mu = mu, kappa = kappa, psi = psi)
  omega_grid <- jp_normalize_unit_matrix(
    omega_grid,
    arg_name = "`omega_grid`",
    min_ncol = length(params$mu)
  )
  rho <- as.numeric(omega_grid %*% params$mu)
  cdf_object <- build_jp_projected_cdf_matrix(
    rho = rho,
    q = params$q,
    alpha = params$alpha,
    beta = params$beta,
    n_u = n_u,
    n_delta = n_delta
  )
  thresholds <- sphere_distance_to_dot_threshold(t_grid, distance_type = distance_type)
  thresholds <- pmin(pmax(thresholds, -1), 1)

  output <- t(vapply(seq_along(rho), function(i) {
    jp_interpolate_upper_tail_linear(
      threshold = thresholds,
      x_grid = cdf_object$u,
      cdf_grid = cdf_object$cdf[, i]
    )
  }, numeric(length(t_grid))))

  output[, t_grid <= 0] <- 0
  if (identical(distance_type, "geodesic")) {
    output[, t_grid >= pi] <- 1
  } else {
    output[, t_grid >= 2] <- 1
  }
  output
}

run_hvmf_case <- function(replicate_id) {
  angles <- seq(-0.9, 0.9, length.out = 60L)
  data <- t(vapply(seq_along(angles), function(i) {
    hyperboloid_point_benchmark(0.55 + 0.01 * i, angles[[i]])
  }, numeric(3)))
  theta <- hvmf_mle_h2(data)
  theta$kappa <- 80
  distance_matrix <- hvmf_distance_matrix(data, data)

  linear_time <- benchmark_elapsed_seconds(
    linear_matrix <- hvmf_cvm_profile_matrix_tabulated_linear(
      data = data,
      theta = theta,
      grid_size = 4097L,
      distance_matrix = distance_matrix
    )
  )
  spline_time <- benchmark_elapsed_seconds(
    spline_matrix <- hvmf_cvm_profile_matrix_tabulated(
      data = data,
      theta = theta,
      grid_size = 4097L,
      distance_matrix = distance_matrix
    )
  )

  data.frame(
    model = "hvmf",
    benchmark = "cvm_profile_matrix",
    replicate_id = replicate_id,
    linear_seconds = linear_time,
    spline_seconds = spline_time,
    speedup = linear_time / spline_time,
    max_abs_diff = max(abs(linear_matrix - spline_matrix)),
    stringsAsFactors = FALSE
  )
}

run_jp_case <- function(replicate_id) {
  mu <- c(0, 0, 1)
  omega_grid <- generate_canonical_lattice(60L, dim = 3)
  t_grid <- seq(1e-4, pi - 1e-4, length.out = 200L)

  linear_time <- benchmark_elapsed_seconds(
    linear_matrix <- distance_profile_jp_grid_linear_reference(
      omega_grid = omega_grid,
      mu = mu,
      kappa = 2,
      psi = 0.5,
      t_grid = t_grid,
      distance_type = "geodesic",
      n_u = 4097L,
      n_delta = 1025L
    )
  )
  spline_time <- benchmark_elapsed_seconds(
    spline_matrix <- distance_profile_jp_grid(
      omega_grid = omega_grid,
      mu = mu,
      kappa = 2,
      psi = 0.5,
      t_grid = t_grid,
      distance_type = "geodesic",
      n_u = 4097L,
      n_delta = 1025L
    )
  )

  data.frame(
    model = "jp",
    benchmark = "ks_profile_grid",
    replicate_id = replicate_id,
    linear_seconds = linear_time,
    spline_seconds = spline_time,
    speedup = linear_time / spline_time,
    max_abs_diff = max(abs(linear_matrix - spline_matrix)),
    stringsAsFactors = FALSE
  )
}

summarize_spline_benchmark <- function(raw_df) {
  do.call(rbind, lapply(split(raw_df, raw_df$model), function(df) {
    data.frame(
      model = df$model[[1L]],
      benchmark = df$benchmark[[1L]],
      n_replicates = nrow(df),
      median_linear_seconds = stats::median(df$linear_seconds),
      median_spline_seconds = stats::median(df$spline_seconds),
      median_speedup = stats::median(df$speedup),
      max_abs_diff_max = max(df$max_abs_diff),
      max_abs_diff_median = stats::median(df$max_abs_diff),
      stringsAsFactors = FALSE
    )
  }))
}

run_splinefun_tabulation_benchmark <- function(output_root = file.path("output", "tabulation_splinefun_benchmark"),
                                               n_repeats = 6L,
                                               n_cores = 12L) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  cases <- list(run_hvmf_case, run_jp_case)
  jobs <- expand.grid(case_id = seq_along(cases), replicate_id = seq_len(as.integer(n_repeats)))

  rows <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    case_fun <- cases[[jobs$case_id[[i]]]]
    case_fun(jobs$replicate_id[[i]])
  }, mc.cores = as.integer(n_cores))

  raw_df <- do.call(rbind, rows)
  summary_df <- summarize_spline_benchmark(raw_df)

  raw_csv <- file.path(output_root, "benchmark_raw.csv")
  summary_csv <- file.path(output_root, "benchmark_summary.csv")
  utils::write.csv(raw_df, raw_csv, row.names = FALSE)
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  list(raw = raw_df, summary = summary_df, raw_csv = raw_csv, summary_csv = summary_csv)
}

parse_spline_benchmark_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else TRUE
    out[[key]] <- value
  }
  out
}

if (sys.nframe() == 0L) {
  args <- parse_spline_benchmark_args(commandArgs(trailingOnly = TRUE))
  result <- run_splinefun_tabulation_benchmark(
    output_root = args$output_root %||% file.path("output", "tabulation_splinefun_benchmark"),
    n_repeats = as.integer(args$n_repeats %||% 6L),
    n_cores = as.integer(args$n_cores %||% 12L)
  )
  cat("Raw CSV:", result$raw_csv, "\n")
  cat("Summary CSV:", result$summary_csv, "\n")
}
