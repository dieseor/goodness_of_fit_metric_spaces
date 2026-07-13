#!/usr/bin/env Rscript

resolve_sunspots_gof_path <- function(...) {
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

utils_path_sunspots_gof <- resolve_sunspots_gof_path("utils.R")
model_specs_path_sunspots_gof <- resolve_sunspots_gof_path("bootstrap", "model_specs.R")
bootstrap_path_sunspots_gof <- resolve_sunspots_gof_path("bootstrap", "multiplier_bootstrap.R")
spec_path_sunspots_gof <- resolve_sunspots_gof_path("bootstrap", "small_circle_symmetric_mixture2_model_spec.R")
prep_path_sunspots_gof <- resolve_sunspots_gof_path("real_data", "sunspots", "sunspots.R")

source(utils_path_sunspots_gof)
source(model_specs_path_sunspots_gof)
source(spec_path_sunspots_gof)
source(bootstrap_path_sunspots_gof)

write_lines_if_possible_sunspots <- function(lines, path) {
  writeLines(as.character(lines), con = path)
  invisible(path)
}

sunspots_projected_gof <- function(z, cdf_fun, grid_size = 1001L) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  if (length(z) == 0L) {
    stop("`z` must contain at least one finite projected value.")
  }

  ecdf_z <- stats::ecdf(z)
  z_grid <- sort(unique(c(seq(-1, 1, length.out = as.integer(grid_size)), z)))
  fitted_cdf <- cdf_fun(z_grid)
  empirical_cdf <- ecdf_z(z_grid)

  list(
    ks = max(abs(empirical_cdf - fitted_cdf)),
    cvm = mean((empirical_cdf - fitted_cdf)^2),
    z_grid = z_grid,
    empirical_cdf = empirical_cdf,
    fitted_cdf = fitted_cdf
  )
}

plot_sunspots_projected_diagnostics <- function(z,
                                                cdf_grid_df,
                                                output_dir) {
  cdf_path <- file.path(output_dir, "sunspots_cycle24_north_pole_projected_ecdf_cdf.png")

  grDevices::png(cdf_path, width = 1200, height = 900, res = 140)
  plot(
    NA,
    NA,
    xlim = c(-1, 1),
    ylim = c(0, 1),
    xlab = "z = <e3, x> = x3",
    ylab = "CDF",
    main = "Cycle 24 sunspots: projected CDF on the north-pole axis"
  )
  z_emp <- sort(unique(c(-1, z, 1)))
  lines(z_emp, stats::ecdf(z)(z_emp), type = "s", lwd = 2, col = "black")
  lines(cdf_grid_df$z_grid, cdf_grid_df$cdf, lwd = 2, col = "#1f78b4")
  legend(
    "topleft",
    legend = c("ECDF", "Fitted CDF"),
    col = c("black", "#1f78b4"),
    lwd = c(2, 2),
    lty = c(1, 1),
    bty = "n"
  )
  grid(col = "#d9d9d9")
  grDevices::dev.off()

  list(cdf_path = cdf_path)
}

extract_sunspots_stat_inference <- function(result, statistic_name, field) {
  inference <- result$inference[[statistic_name]]
  if (is.null(inference)) {
    return(NA_real_)
  }

  candidate_names <- switch(
    field,
    observed = c("observed", "observed_statistic", "statistic", "Tn"),
    p_value = c("p_value", "p.value", "pvalue", "p"),
    stop("Unknown inference field.")
  )

  for (candidate_name in candidate_names) {
    value <- NULL
    if (is.list(inference) && !is.null(inference[[candidate_name]])) {
      value <- inference[[candidate_name]]
    } else if (is.data.frame(inference) && candidate_name %in% names(inference)) {
      value <- inference[[candidate_name]][[1L]]
    }
    if (!is.null(value) && length(value) >= 1L && is.finite(as.numeric(value[[1L]]))) {
      return(as.numeric(value[[1L]]))
    }
  }

  NA_real_
}

normalize_sunspots_statistics <- function(statistics) {
  if (length(statistics) == 0L || all(is.na(statistics))) {
    stop("At least one statistic must be supplied via --statistics or --statistic.")
  }

  statistics <- unlist(strsplit(paste(statistics, collapse = ","), ",", fixed = TRUE), use.names = FALSE)
  statistics <- trimws(statistics)
  statistics <- statistics[nzchar(statistics)]
  statistics <- tolower(statistics)

  supported <- c("ks", "cvm")
  if (any(!statistics %in% supported)) {
    invalid <- unique(statistics[!statistics %in% supported])
    stop(
      sprintf(
        "Unsupported statistic(s): %s. Supported values are: ks, cvm.",
        paste(invalid, collapse = ", ")
      )
    )
  }

  unique(statistics)
}

parse_sunspots_gof_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  help_text <- paste(
    "Usage:",
    "  Rscript real_data/sunspots/run_sunspots_cycle24_small_circle_symmetric_mixture_gof.R [options]",
    "",
    "Options:",
    "  --help                  Show this help message and exit.",
    "  --input_csv=PATH        Input CSV. Default: real_data/sunspots/output/sunspots_cycle24_s2_all.csv",
    "  --output_dir=PATH       Output directory. Default: real_data/sunspots/output/cycle24_small_circle_symmetric_mixture",
    "  --B=INTEGER             Number of bootstrap replicates.",
    "  --n_cores=INTEGER       Number of worker cores.",
    "  --seed=INTEGER          Random seed.",
    "  --M_value=INTEGER       Canonical lattice resolution for KS.",
    "  --ks_t_points=INTEGER   Number of t-grid points for KS.",
    "  --distance_type=TYPE    Distance type passed to the model spec.",
    "  --statistics=LIST       Comma-separated statistics from {ks,cvm}.",
    "  --statistic=LIST        Alias for --statistics.",
    "  --keep_bootstrap_thetas=TRUE|FALSE",
    "  --verbose_timing=TRUE|FALSE",
    sep = "\n"
  )

  if (length(args) == 0L) {
    return(list())
  }

  if (any(args %in% c("--help", "-h"))) {
    cat(help_text, sep = "\n")
    quit(save = "no", status = 0L)
  }

  parsed <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unrecognized positional argument: %s\n\n%s", arg, help_text))
    }

    key_value <- sub("^--", "", arg)
    parts <- strsplit(key_value, "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else NULL

    if (!nzchar(key)) {
      stop(sprintf("Malformed option: %s\n\n%s", arg, help_text))
    }
    if (is.null(value)) {
      stop(sprintf("Option '--%s' requires the form --%s=value.\n\n%s", key, key, help_text))
    }
    if (key %in% names(parsed)) {
      stop(sprintf("Option '--%s' was supplied more than once.", key))
    }

    parsed[[key]] <- value
  }

  if ("statistics" %in% names(parsed) && "statistic" %in% names(parsed)) {
    stop("Use only one of --statistics or --statistic, not both.")
  }

  unknown <- setdiff(names(parsed), c(
    "input_csv", "output_dir", "B", "n_cores", "seed", "M_value",
    "ks_t_points", "distance_type", "statistics", "statistic",
    "keep_bootstrap_thetas", "verbose_timing"
  ))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown option(s): %s\n\n%s",
        paste(paste0("--", unknown), collapse = ", "),
        help_text
      )
    )
  }

  parse_integer <- function(value, key) {
    if (!grepl("^-?[0-9]+$", value)) {
      stop(sprintf("Option '--%s' must be an integer. Received: %s", key, value))
    }
    as.integer(value)
  }

  parse_logical <- function(value, key) {
    value_lower <- tolower(value)
    if (!value_lower %in% c("true", "false")) {
      stop(sprintf("Option '--%s' must be TRUE or FALSE. Received: %s", key, value))
    }
    identical(value_lower, "true")
  }

  parsed_args <- list()
  if ("input_csv" %in% names(parsed)) parsed_args$input_csv <- parsed$input_csv
  if ("output_dir" %in% names(parsed)) parsed_args$output_dir <- parsed$output_dir
  if ("B" %in% names(parsed)) parsed_args$B <- parse_integer(parsed$B, "B")
  if ("n_cores" %in% names(parsed)) parsed_args$n_cores <- parse_integer(parsed$n_cores, "n_cores")
  if ("seed" %in% names(parsed)) parsed_args$seed <- parse_integer(parsed$seed, "seed")
  if ("M_value" %in% names(parsed)) parsed_args$M_value <- parse_integer(parsed$M_value, "M_value")
  if ("ks_t_points" %in% names(parsed)) parsed_args$ks_t_points <- parse_integer(parsed$ks_t_points, "ks_t_points")
  if ("distance_type" %in% names(parsed)) parsed_args$distance_type <- parsed$distance_type
  if ("statistics" %in% names(parsed)) parsed_args$statistics <- normalize_sunspots_statistics(parsed$statistics)
  if ("statistic" %in% names(parsed)) parsed_args$statistics <- normalize_sunspots_statistics(parsed$statistic)
  if ("keep_bootstrap_thetas" %in% names(parsed)) {
    parsed_args$keep_bootstrap_thetas <- parse_logical(parsed$keep_bootstrap_thetas, "keep_bootstrap_thetas")
  }
  if ("verbose_timing" %in% names(parsed)) {
    parsed_args$verbose_timing <- parse_logical(parsed$verbose_timing, "verbose_timing")
  }

  parsed_args
}

run_sunspots_cycle24_small_circle_symmetric_mixture_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle24_small_circle_symmetric_mixture"),
    statistics = c("ks", "cvm"),
    B = 5000L,
    n_cores = 6L,
    seed = 20260603L,
    M_value = 60L,
    ks_t_points = 200L,
    distance_type = "geodesic",
    keep_bootstrap_thetas = FALSE,
    verbose_timing = TRUE,
    control = list(
      small_circle_symmetric_mixture2_profile_method = "legendre",
      small_circle_symmetric_mixture2_L_max = 200L,
      small_circle_symmetric_mixture2_quad_n = 400L,
      small_circle_symmetric_mixture2_tol = 1e-10,
      small_circle_symmetric_mixture2_optim_control = list(maxit = 400L, reltol = 1e-9)
    )) {
  statistics <- normalize_sunspots_statistics(statistics)

  if (!file.exists(input_csv)) {
    source(prep_path_sunspots_gof)
  }
  if (!file.exists(input_csv)) {
    stop(sprintf("Input CSV not found after running preparation script: %s", input_csv))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_lines_if_possible_sunspots(capture.output(sessionInfo()), file.path(output_dir, "sunspots_cycle24_sessionInfo.txt"))
  block_size <- as.integer(control$small_circle_symmetric_mixture2_cvm_block_size %||% 128L)

  timing_log <- function(stage, event, extra = NULL) {
    if (!isTRUE(verbose_timing)) {
      return(invisible(NULL))
    }
    suffix <- if (is.null(extra) || identical(extra, "")) "" else paste0(" | ", extra)
    message(sprintf(
      "[sunspots-cycle24] %s %s | statistic=%s | B=%d | n_cores=%d | block_size=%d%s",
      stage,
      event,
      paste(statistics, collapse = ","),
      as.integer(B),
      as.integer(n_cores),
      block_size,
      suffix
    ))
  }

  sunspots_df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  x <- as.matrix(sunspots_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)

  timing_log("MLE observed", "start")
  mle_start <- proc.time()[["elapsed"]]
  theta_hat <- small_circle_symmetric_mixture2_mle_s2_weighted(
    x = x,
    control = control
  )
  mle_elapsed <- proc.time()[["elapsed"]] - mle_start
  timing_log(
    "MLE observed",
    "end",
    sprintf(
      "elapsed=%.3fs | mu=(%.6f, %.6f, %.6f) | kappa=%.6f | nu=%.6f",
      mle_elapsed,
      theta_hat$mu[[1L]],
      theta_hat$mu[[2L]],
      theta_hat$mu[[3L]],
      theta_hat$kappa,
      theta_hat$nu
    )
  )

  prep_elapsed <- 0
  cvm_observed_prep <- NULL
  if ("cvm" %in% statistics) {
    timing_log("CvM observed prep", "start")
    prep_start <- proc.time()[["elapsed"]]
    spec_for_prep <- make_small_circle_symmetric_mixture2_spec(distance_type = distance_type)
    cvm_observed_prep <- spec_for_prep$cvm_prepare(x, theta_hat, control = control)
    prep_elapsed <- proc.time()[["elapsed"]] - prep_start
    timing_log(
      "CvM observed prep",
      "end",
      sprintf("elapsed=%.3fs | observed_stat=%.6f", prep_elapsed, cvm_observed_prep$statistic)
    )
  }

  log_density <- d_sph_small_circle_symmetric_mixture2_s2(
    x = x,
    mu = theta_hat$mu,
    kappa = theta_hat$kappa,
    nu = theta_hat$nu,
    log = TRUE
  )
  loglik <- sum(log_density)
  n <- nrow(x)
  k <- 4L
  aic <- 2 * k - 2 * loglik
  bic <- log(n) * k - 2 * loglik

  north_pole <- c(0, 0, 1)
  z <- pmin(pmax(as.numeric(x %*% north_pole), -1), 1)
  projected_fit <- sunspots_projected_gof(
    z = z,
    cdf_fun = function(z_grid) {
      1 - distance_profile_small_circle_symmetric_mixture2(
        omega = north_pole,
        t_values = acos(pmin(pmax(z_grid, -1), 1)),
        mu = theta_hat$mu,
        kappa = theta_hat$kappa,
        nu = theta_hat$nu,
        distance_type = "geodesic",
        method = control$small_circle_symmetric_mixture2_profile_method %||% "legendre",
        l_max = as.integer(control$small_circle_symmetric_mixture2_L_max %||% 200L),
        quad_n = as.integer(control$small_circle_symmetric_mixture2_quad_n %||% 400L),
        tol = as.numeric(control$small_circle_symmetric_mixture2_tol %||% 1e-10)
      )
    }
  )

  cdf_grid_df <- data.frame(
    z_grid = projected_fit$z_grid,
    cdf = projected_fit$fitted_cdf,
    empirical_cdf = projected_fit$empirical_cdf,
    stringsAsFactors = FALSE
  )

  spec <- make_small_circle_symmetric_mixture2_spec(distance_type = distance_type)
  ks_grid <- make_sample_unique_distance_ks_grid()

  timing_log("bootstrap", "start")
  bootstrap_start <- proc.time()[["elapsed"]]
  gof_result <- multiplier_bootstrap_gof(
    data = x,
    spec = spec,
    null = list(type = "composite"),
    statistics = statistics,
    ks_grid = if ("ks" %in% statistics) ks_grid else NULL,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    observed_theta_hat = theta_hat,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = isTRUE(keep_bootstrap_thetas)
    ),
    control = control
  )
  bootstrap_elapsed <- proc.time()[["elapsed"]] - bootstrap_start
  timing_log("bootstrap", "end", sprintf("elapsed=%.3fs", bootstrap_elapsed))

  total_elapsed <- mle_elapsed + prep_elapsed + bootstrap_elapsed
  timing_log("total", "end", sprintf("elapsed=%.3fs", total_elapsed))

  result_rows <- do.call(rbind, lapply(statistics, function(stat_name) {
    data.frame(
      cycle = 24L,
      dataset = "cycle24_all",
      model = "small_circle_symmetric_mixture2",
      mixture_components = 2L,
      statistic_type = stat_name,
      n = n,
      loglik = loglik,
      aic = aic,
      bic = bic,
      test_statistic = extract_sunspots_stat_inference(gof_result, stat_name, "observed"),
      p_value_raw = extract_sunspots_stat_inference(gof_result, stat_name, "p_value"),
      mu_1 = theta_hat$mu[[1L]],
      mu_2 = theta_hat$mu[[2L]],
      mu_3 = theta_hat$mu[[3L]],
      kappa_hat = theta_hat$kappa,
      nu_hat = theta_hat$nu,
      stringsAsFactors = FALSE
    )
  }))

  theta_df <- data.frame(
    mu_1 = theta_hat$mu[[1L]],
    mu_2 = theta_hat$mu[[2L]],
    mu_3 = theta_hat$mu[[3L]],
    kappa_hat = theta_hat$kappa,
    nu_hat = theta_hat$nu,
    loglik = loglik,
    aic = aic,
    bic = bic,
    stringsAsFactors = FALSE
  )

  utils::write.csv(result_rows, file = file.path(output_dir, "sunspots_cycle24_small_circle_symmetric_mixture_gof_results.csv"), row.names = FALSE)
  utils::write.csv(theta_df, file = file.path(output_dir, "sunspots_cycle24_theta_hat.csv"), row.names = FALSE)
  utils::write.csv(cdf_grid_df, file = file.path(output_dir, "sunspots_cycle24_north_pole_projected_cdf_grid.csv"), row.names = FALSE)
  utils::write.csv(data.frame(z = z), file = file.path(output_dir, "sunspots_cycle24_north_pole_projected_data.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(
      mle_observed_seconds = mle_elapsed,
      cvm_observed_prep_seconds = prep_elapsed,
      bootstrap_seconds = bootstrap_elapsed,
      total_seconds = total_elapsed,
      B = as.integer(B),
      n_cores = as.integer(n_cores),
      block_size = block_size,
      statistic_type = paste(statistics, collapse = ","),
      stringsAsFactors = FALSE
    ),
    file = file.path(output_dir, "sunspots_cycle24_timing_summary.csv"),
    row.names = FALSE
  )
  if (isTRUE(keep_bootstrap_thetas) && !is.null(gof_result$bootstrap$theta_star)) {
    theta_star_df <- do.call(rbind, lapply(seq_along(gof_result$bootstrap$theta_star), function(i) {
      theta_star <- gof_result$bootstrap$theta_star[[i]]
      data.frame(
        replicate = i,
        mu_1 = theta_star$mu[[1L]],
        mu_2 = theta_star$mu[[2L]],
        mu_3 = theta_star$mu[[3L]],
        kappa_hat = theta_star$kappa,
        nu_hat = theta_star$nu,
        stringsAsFactors = FALSE
      )
    }))
    utils::write.csv(theta_star_df, file = file.path(output_dir, "sunspots_cycle24_theta_star.csv"), row.names = FALSE)
  }

  plot_sunspots_projected_diagnostics(
    z = z,
    cdf_grid_df = cdf_grid_df,
    output_dir = output_dir
  )

  attr(result_rows, "timing_summary") <- list(
    mle_observed_seconds = mle_elapsed,
    cvm_observed_prep_seconds = prep_elapsed,
    bootstrap_seconds = bootstrap_elapsed,
    total_seconds = total_elapsed,
    B = as.integer(B),
    n_cores = as.integer(n_cores),
    block_size = block_size
  )
  attr(result_rows, "theta_hat") <- theta_hat
  attr(result_rows, "theta_star") <- gof_result$bootstrap$theta_star %||% NULL

  result_rows
}

if (sys.nframe() == 0L) {
  cli_args <- parse_sunspots_gof_args()
  print(do.call(run_sunspots_cycle24_small_circle_symmetric_mixture_gof, cli_args))
}
