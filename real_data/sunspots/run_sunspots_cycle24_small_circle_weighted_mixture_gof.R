#!/usr/bin/env Rscript

resolve_sunspots_weighted_gof_path <- function(...) {
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

utils_path_sunspots_weighted_gof <- resolve_sunspots_weighted_gof_path("utils.R")
model_specs_path_sunspots_weighted_gof <- resolve_sunspots_weighted_gof_path("bootstrap", "model_specs.R")
bootstrap_path_sunspots_weighted_gof <- resolve_sunspots_weighted_gof_path("bootstrap", "multiplier_bootstrap.R")
spec_path_sunspots_weighted_gof <- resolve_sunspots_weighted_gof_path("bootstrap", "small_circle_weighted_mixture2_model_spec.R")
prep_path_sunspots_weighted_gof <- resolve_sunspots_weighted_gof_path("real_data", "sunspots", "sunspots.R")

source(utils_path_sunspots_weighted_gof)
source(model_specs_path_sunspots_weighted_gof)
source(spec_path_sunspots_weighted_gof)
source(bootstrap_path_sunspots_weighted_gof)

normalize_sunspots_weighted_statistics <- function(statistics) {
  statistics <- unique(tolower(trimws(unlist(strsplit(paste(statistics, collapse = ","), ",", fixed = TRUE), use.names = FALSE))))
  statistics <- statistics[nzchar(statistics)]
  if (length(statistics) == 0L) {
    stop("At least one statistic must be supplied.")
  }
  if ("cvm" %in% statistics) {
    stop("CvM not implemented for weighted mixture yet; use KS.")
  }
  if (!all(statistics %in% "ks")) {
    stop("Weighted sunspots runner currently supports only KS.")
  }
  statistics
}

parse_sunspots_weighted_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) {
    return(list())
  }

  help_text <- paste(
    "Usage:",
    "  Rscript real_data/sunspots/run_sunspots_cycle24_small_circle_weighted_mixture_gof.R [options]",
    "",
    "Options:",
    "  --help",
    "  --input_csv=PATH",
    "  --output_dir=PATH",
    "  --B=INTEGER",
    "  --n_cores=INTEGER",
    "  --seed=INTEGER",
    "  --statistics=KS",
    "  --statistic=KS",
    "  --distance_type=geodesic|chordal",
    sep = "\n"
  )

  if (any(args %in% c("--help", "-h"))) {
    cat(help_text, sep = "\n")
    quit(save = "no", status = 0L)
  }

  parsed <- list()
  for (arg in args) {
    key_value <- sub("^--", "", arg)
    parts <- strsplit(key_value, "=", fixed = TRUE)[[1L]]
    key <- parts[[1L]]
    value <- if (length(parts) >= 2L) paste(parts[-1L], collapse = "=") else NULL
    if (is.null(value)) {
      stop(sprintf("Option '--%s' requires the form --%s=value.", key, key))
    }
    parsed[[key]] <- value
  }

  parse_integer <- function(value, key) {
    if (!grepl("^-?[0-9]+$", value)) {
      stop(sprintf("Option '--%s' must be an integer.", key))
    }
    as.integer(value)
  }

  out <- list()
  if ("input_csv" %in% names(parsed)) out$input_csv <- parsed$input_csv
  if ("output_dir" %in% names(parsed)) out$output_dir <- parsed$output_dir
  if ("B" %in% names(parsed)) out$B <- parse_integer(parsed$B, "B")
  if ("n_cores" %in% names(parsed)) out$n_cores <- parse_integer(parsed$n_cores, "n_cores")
  if ("seed" %in% names(parsed)) out$seed <- parse_integer(parsed$seed, "seed")
  if ("distance_type" %in% names(parsed)) out$distance_type <- parsed$distance_type
  if ("statistics" %in% names(parsed)) out$statistics <- normalize_sunspots_weighted_statistics(parsed$statistics)
  if ("statistic" %in% names(parsed)) out$statistics <- normalize_sunspots_weighted_statistics(parsed$statistic)
  out
}

sunspots_weighted_default_symmetric_start <- function() {
  small_circle_symmetric_mixture2_canonicalize_theta(list(
    mu = c(0.002746, 0.002837, -0.999992),
    kappa = 25.045825,
    nu = 0.251825
  ))
}

run_sunspots_cycle24_small_circle_weighted_mixture_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle24_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle24_small_circle_weighted_mixture"),
    dataset_label = "cycle24_all",
    statistics = "ks",
    B = 5000L,
    n_cores = 6L,
    seed = 123L,
    distance_type = "geodesic",
    control = list(
      small_circle_weighted_mixture2_profile_method = "legendre",
      small_circle_weighted_mixture2_L_max = 200L,
      small_circle_weighted_mixture2_quad_n = 400L,
      small_circle_weighted_mixture2_tol = 1e-10,
      small_circle_weighted_mixture2_optim_control = list(maxit = 400L, reltol = 1e-9),
      small_circle_weighted_mixture2_n_starts = 8L
    )) {
  statistics <- normalize_sunspots_weighted_statistics(statistics)

  if (!file.exists(input_csv)) {
    source(prep_path_sunspots_weighted_gof)
  }
  if (!file.exists(input_csv)) {
    stop(sprintf("Input CSV not found: %s", input_csv))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sunspots_cycle24_weighted_sessionInfo.txt"))

  sunspots_df <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  x <- as.matrix(sunspots_df[, c("x1", "x2", "x3")])
  x <- jp_normalize_unit_matrix(x, arg_name = "`x`", min_ncol = 3L)

  control$small_circle_weighted_mixture2_symmetric_start_theta <- sunspots_weighted_default_symmetric_start()

  theta_hat <- small_circle_weighted_mixture2_mle_s2_weighted(x = x, control = control)
  log_density <- d_sph_small_circle_weighted_mixture2_s2(
    x = x,
    mu = theta_hat$mu,
    pi = theta_hat$pi,
    kappa1 = theta_hat$kappa1,
    nu1 = theta_hat$nu1,
    kappa2 = theta_hat$kappa2,
    nu2 = theta_hat$nu2,
    log = TRUE
  )
  loglik <- sum(log_density)
  n <- nrow(x)
  k <- 7L
  aic <- 2 * k - 2 * loglik
  bic <- log(n) * k - 2 * loglik

  spec <- make_small_circle_weighted_mixture2_spec(distance_type = distance_type)
  ks_grid <- make_sample_unique_distance_ks_grid()
  gof_result <- multiplier_bootstrap_gof(
    data = x,
    spec = spec,
    null = list(type = "composite"),
    statistics = statistics,
    ks_grid = ks_grid,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    observed_theta_hat = theta_hat,
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = TRUE
    ),
    control = control
  )

  result_rows <- data.frame(
    cycle = 24L,
    dataset = dataset_label,
    model = "small_circle_weighted_mixture2",
    mixture_components = 2L,
    statistic_type = "ks",
    n = n,
    loglik = loglik,
    aic = aic,
    bic = bic,
    test_statistic = gof_result$inference$ks$observed,
    p_value_raw = gof_result$inference$ks$p_value,
    pi_hat = theta_hat$pi,
    mu_1 = theta_hat$mu[[1L]],
    mu_2 = theta_hat$mu[[2L]],
    mu_3 = theta_hat$mu[[3L]],
    kappa1_hat = theta_hat$kappa1,
    nu1_hat = theta_hat$nu1,
    kappa2_hat = theta_hat$kappa2,
    nu2_hat = theta_hat$nu2,
    stringsAsFactors = FALSE
  )

  theta_df <- data.frame(
    pi_hat = theta_hat$pi,
    mu_1 = theta_hat$mu[[1L]],
    mu_2 = theta_hat$mu[[2L]],
    mu_3 = theta_hat$mu[[3L]],
    kappa1_hat = theta_hat$kappa1,
    nu1_hat = theta_hat$nu1,
    kappa2_hat = theta_hat$kappa2,
    nu2_hat = theta_hat$nu2,
    loglik = loglik,
    aic = aic,
    bic = bic,
    stringsAsFactors = FALSE
  )

  utils::write.csv(result_rows, file.path(output_dir, "sunspots_cycle24_small_circle_weighted_mixture_gof_results.csv"), row.names = FALSE)
  utils::write.csv(theta_df, file.path(output_dir, "sunspots_cycle24_theta_hat.csv"), row.names = FALSE)

  message(sprintf("n = %d", n))
  message(sprintf("loglik = %.6f | AIC = %.6f | BIC = %.6f", loglik, aic, bic))
  message(sprintf("pi_hat = %.6f", theta_hat$pi))
  message(sprintf("mu_hat = (%.6f, %.6f, %.6f)", theta_hat$mu[[1L]], theta_hat$mu[[2L]], theta_hat$mu[[3L]]))
  message(sprintf("kappa1_hat = %.6f | nu1_hat = %.6f | kappa2_hat = %.6f | nu2_hat = %.6f",
                  theta_hat$kappa1, theta_hat$nu1, theta_hat$kappa2, theta_hat$nu2))
  message(sprintf("KS statistic = %.6f | p-value = %.6f", result_rows$test_statistic[[1L]], result_rows$p_value_raw[[1L]]))

  result_rows
}

if (sys.nframe() == 0L) {
  cli_args <- parse_sunspots_weighted_args()
  print(do.call(run_sunspots_cycle24_small_circle_weighted_mixture_gof, cli_args))
}
