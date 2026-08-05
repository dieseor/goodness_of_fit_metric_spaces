#!/usr/bin/env Rscript

resolve_sunspots_temporal_uniform_beta_path <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_sunspots_temporal_uniform_beta_path(
  "real_data",
  "sunspots",
  "sunspots_cycle23_joint_time_models_parsimonious.R"
))

normalize_sunspots_temporal_statistics <- function(
    statistics = c("ks", "cvm")) {
  statistics <- unique(tolower(as.character(statistics)))
  if (length(statistics) == 0L ||
      any(!statistics %in% c("ks", "cvm"))) {
    stop("`statistics` must contain one or both of 'ks' and 'cvm'.")
  }
  statistics
}

sunspots_temporal_uniform_beta_statistics <- function(
    s,
    eta,
    control = list()) {
  s <- sort(as.numeric(s))
  if (length(s) < 2L ||
      any(!is.finite(s)) ||
      any(s <= 0 | s >= 1)) {
    stop("`s` must contain at least two finite values in (0, 1).")
  }

  eta <- sunspots_joint_parsimonious_time_validate_eta(
    eta,
    time_model = "uniform_beta",
    control = control
  )
  fitted_cdf <- sunspots_joint_parsimonious_time_cdf(
    s,
    eta,
    time_model = "uniform_beta",
    control = control
  )
  if (any(!is.finite(fitted_cdf))) {
    stop("The fitted temporal CDF contains non-finite values.")
  }

  n <- length(s)
  index <- seq_len(n)
  ks_plus <- max(index / n - fitted_cdf)
  ks_minus <- max(fitted_cdf - (index - 1) / n)
  ks <- sqrt(n) * max(ks_plus, ks_minus)
  cvm <- 1 / (12 * n) + sum(
    (fitted_cdf - (2 * index - 1) / (2 * n))^2
  )

  list(
    ks = as.numeric(ks),
    cvm = as.numeric(cvm),
    sorted_s = s,
    fitted_cdf = fitted_cdf,
    ks_plus = as.numeric(ks_plus),
    ks_minus = as.numeric(ks_minus)
  )
}

sunspots_temporal_uniform_beta_refit <- function(
    s,
    control = list()) {
  first <- suppressWarnings(try(
    fit_sunspots_joint_parsimonious_time(
      s,
      time_model = "uniform_beta",
      control = control
    ),
    silent = TRUE
  ))
  first_ok <- !inherits(first, "try-error") &&
    isTRUE(first$selected_converged) &&
    is.finite(first$loglik)

  if (first_ok) return(first)

  retry_control <- utils::modifyList(
    control,
    list(
      parsimonious_time_n_starts = as.integer(
        control$temporal_gof_retry_n_starts %||% 9L
      ),
      parsimonious_time_nelder_mead_control =
        control$temporal_gof_retry_nelder_mead_control %||%
        list(maxit = 6000L, reltol = 1e-11),
      parsimonious_time_optim_control =
        control$temporal_gof_retry_optim_control %||%
        list(maxit = 3000L, factr = 1e7, pgtol = 1e-9)
    )
  )
  second <- suppressWarnings(try(
    fit_sunspots_joint_parsimonious_time(
      s,
      time_model = "uniform_beta",
      control = retry_control
    ),
    silent = TRUE
  ))
  second_ok <- !inherits(second, "try-error") &&
    isTRUE(second$selected_converged) &&
    is.finite(second$loglik)

  if (!second_ok) {
    first_message <- if (inherits(first, "try-error")) {
      as.character(first)
    } else {
      as.character(first$opt$message %||% "non-converged first fit")
    }
    second_message <- if (inherits(second, "try-error")) {
      as.character(second)
    } else {
      as.character(second$opt$message %||% "non-converged retry")
    }
    stop(
      sprintf(
        paste(
          "Uniform+Beta temporal refit failed.",
          "First attempt: %s. Retry: %s."
        ),
        first_message,
        second_message
      ),
      call. = FALSE
    )
  }

  second
}

sunspots_temporal_uniform_beta_bootstrap <- function(
    s,
    eta_hat,
    B = 1000L,
    n_cores = 1L,
    seed = 20260714L,
    control = list()) {
  s <- as.numeric(s)
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  seed <- as.integer(seed)

  if (length(B) != 1L || !is.finite(B) || B < 1L) {
    stop("`B` must be a positive integer.")
  }
  if (length(n_cores) != 1L ||
      !is.finite(n_cores) ||
      n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be one finite integer.")
  }

  eta_hat <- sunspots_joint_parsimonious_time_validate_eta(
    eta_hat,
    time_model = "uniform_beta",
    control = control
  )
  n <- length(s)

  replicate_seeds <- sunspots_joint_with_seed(
    seed,
    sample.int(.Machine$integer.max, B, replace = FALSE)
  )

  worker <- function(replicate_index) {
    set.seed(replicate_seeds[[replicate_index]])
    simulated <- sample_sunspots_joint_parsimonious_time(
      n,
      eta_hat,
      time_model = "uniform_beta",
      control = control
    )
    fit_star <- sunspots_temporal_uniform_beta_refit(
      simulated,
      control = control
    )
    statistic_star <- sunspots_temporal_uniform_beta_statistics(
      simulated,
      fit_star,
      control = control
    )
    fit_summary <- sunspots_joint_parsimonious_time_summary(
      fit_star,
      time_model = "uniform_beta",
      control = control
    )

    data.frame(
      replicate = replicate_index,
      seed = replicate_seeds[[replicate_index]],
      ks = statistic_star$ks,
      cvm = statistic_star$cvm,
      beta_weight = fit_summary$temporal_beta_weight,
      uniform_weight = fit_summary$temporal_uniform_weight,
      alpha = fit_summary$temporal_alpha,
      beta = fit_summary$temporal_beta,
      loglik = fit_star$loglik,
      convergence = as.integer(fit_star$opt$convergence),
      selected_converged = isTRUE(fit_star$selected_converged),
      boundary_weight =
        isTRUE(fit_summary$temporal_boundary_weight),
      boundary_shape_lower =
        isTRUE(fit_summary$temporal_boundary_shape_lower),
      boundary_shape_upper =
        isTRUE(fit_summary$temporal_boundary_shape_upper),
      identification_failure =
        isTRUE(fit_summary$temporal_identification_failure),
      fast_regular =
        isTRUE(fit_summary$temporal_fast_regular),
      stringsAsFactors = FALSE
    )
  }

  results <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(
      seq_len(B),
      worker,
      mc.cores = min(n_cores, B),
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(B), worker)
  }

  failed <- vapply(results, inherits, logical(1L), what = "try-error")
  if (any(failed)) {
    stop(
      sprintf(
        "Temporal bootstrap replicate %d failed: %s",
        which(failed)[[1L]],
        as.character(results[[which(failed)[[1L]]]])
      ),
      call. = FALSE
    )
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

plot_sunspots_temporal_uniform_beta_diagnostics <- function(
    s,
    eta_hat,
    output_path,
    control = list()) {
  s <- as.numeric(s)
  grid <- seq(1e-5, 1 - 1e-5, length.out = 1001L)
  fitted_density <- sunspots_joint_parsimonious_time_density(
    grid,
    eta_hat,
    time_model = "uniform_beta",
    control = control
  )
  fitted_cdf <- sunspots_joint_parsimonious_time_cdf(
    grid,
    eta_hat,
    time_model = "uniform_beta",
    control = control
  )
  pit <- sunspots_joint_parsimonious_time_cdf(
    s,
    eta_hat,
    time_model = "uniform_beta",
    control = control
  )

  grDevices::png(output_path, width = 1800, height = 600, res = 150)
  old_par <- graphics::par(
    mfrow = c(1, 3),
    mar = c(4, 4, 3, 1)
  )
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::hist(
    s,
    breaks = 40,
    freq = FALSE,
    col = "grey85",
    border = "white",
    xlab = "Normalized dequantized time",
    main = "Temporal density"
  )
  graphics::lines(grid, fitted_density, lwd = 2)

  graphics::plot(
    stats::ecdf(s),
    verticals = TRUE,
    do.points = FALSE,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "Normalized dequantized time",
    ylab = "CDF",
    main = "Empirical and fitted CDF"
  )
  graphics::lines(grid, fitted_cdf, lwd = 2)

  stats::qqplot(
    stats::ppoints(length(pit)),
    sort(pit),
    xlim = c(0, 1),
    ylim = c(0, 1),
    pch = 16,
    cex = 0.45,
    xlab = "Uniform quantiles",
    ylab = "Temporal PIT quantiles",
    main = "Temporal PIT Q-Q plot"
  )
  graphics::abline(0, 1, lwd = 2)

  invisible(list(
    grid = grid,
    fitted_density = fitted_density,
    fitted_cdf = fitted_cdf,
    pit = pit
  ))
}

run_sunspots_cycle23_temporal_uniform_beta_gof <- function(
    input_csv = file.path(
      "real_data",
      "sunspots",
      "output",
      "sunspots_cycle23_s2_all.csv"
    ),
    output_dir = file.path(
      "real_data",
      "sunspots",
      "output",
      "cycle23_temporal_uniform_beta_gof"
    ),
    start_date = "1997-06-01",
    end_date = "2006-01-01",
    statistics = c("ks", "cvm"),
    B = 5000L,
    n_cores = 8L,
    dequantization_seed = 20260712L,
    bootstrap_seed = 20260714L,
    control = list()) {
  statistics <- normalize_sunspots_temporal_statistics(statistics)
  B <- as.integer(B)
  n_cores <- as.integer(n_cores)

  total_start <- proc.time()[["elapsed"]]
  retained <- prepare_sunspots_cycle23_joint_time_space_data(
    input_csv = input_csv,
    start_date = start_date,
    end_date = end_date,
    dequantization_seed = dequantization_seed
  )
  s <- retained$s

  mle_start <- proc.time()[["elapsed"]]
  eta_hat <- sunspots_temporal_uniform_beta_refit(
    s,
    control = control
  )
  mle_seconds <- proc.time()[["elapsed"]] - mle_start
  eta_summary <- sunspots_joint_parsimonious_time_summary(
    eta_hat,
    time_model = "uniform_beta",
    control = control
  )

  if (!isTRUE(eta_hat$selected_converged)) {
    stop("The observed uniform+Beta temporal MLE did not converge.")
  }
  if (!isTRUE(eta_summary$temporal_fast_regular)) {
    stop(
      paste(
        "The observed uniform+Beta temporal MLE is on a boundary or",
        "in a nonidentified region; regular GOF calibration is not used."
      ),
      call. = FALSE
    )
  }

  observed <- sunspots_temporal_uniform_beta_statistics(
    s,
    eta_hat,
    control = control
  )

  bootstrap_start <- proc.time()[["elapsed"]]
  bootstrap <- sunspots_temporal_uniform_beta_bootstrap(
    s = s,
    eta_hat = eta_hat,
    B = B,
    n_cores = n_cores,
    seed = bootstrap_seed,
    control = control
  )
  bootstrap_seconds <- proc.time()[["elapsed"]] - bootstrap_start
  total_seconds <- proc.time()[["elapsed"]] - total_start

  temporal_criteria <- sunspots_joint_time_information_criteria(
    loglik = eta_hat$loglik,
    n = length(s),
    n_parameters = 3L
  )

  summary_rows <- lapply(statistics, function(statistic_name) {
    bootstrap_values <- bootstrap[[statistic_name]]
    observed_value <- observed[[statistic_name]]
    data.frame(
      model = "temporal_uniform_beta",
      statistic_type = if (identical(statistic_name, "ks")) {
        "ks_fitted_temporal_cdf"
      } else {
        "cvm_fitted_temporal_cdf"
      },
      bootstrap_method = "reestimated_parametric",
      n = length(s),
      B = B,
      n_cores = n_cores,
      start_date = start_date,
      end_date_exclusive = end_date,
      dequantization_seed = as.integer(dequantization_seed),
      bootstrap_seed = as.integer(bootstrap_seed),
      statistic = observed_value,
      critical_value_0.95 = as.numeric(stats::quantile(
        bootstrap_values,
        probs = 0.95,
        names = FALSE,
        type = 8
      )),
      p_value = (
        1 + sum(bootstrap_values >= observed_value)
      ) / (B + 1),
      temporal_uniform_weight =
        eta_summary$temporal_uniform_weight,
      temporal_beta_weight =
        eta_summary$temporal_beta_weight,
      temporal_alpha = eta_summary$temporal_alpha,
      temporal_beta = eta_summary$temporal_beta,
      temporal_mean = eta_summary$temporal_mean,
      temporal_loglik = eta_hat$loglik,
      temporal_aic = temporal_criteria$aic,
      temporal_bic = temporal_criteria$bic,
      temporal_convergence =
        as.integer(eta_hat$opt$convergence),
      temporal_selected_converged =
        isTRUE(eta_hat$selected_converged),
      temporal_fast_regular =
        isTRUE(eta_summary$temporal_fast_regular),
      temporal_identification_distance =
        eta_summary$temporal_identification_distance,
      bootstrap_boundary_weight_count =
        sum(bootstrap$boundary_weight),
      bootstrap_boundary_shape_count =
        sum(
          bootstrap$boundary_shape_lower |
            bootstrap$boundary_shape_upper
        ),
      bootstrap_identification_failure_count =
        sum(bootstrap$identification_failure),
      bootstrap_nonregular_count =
        sum(!bootstrap$fast_regular),
      temporal_mle_seconds = mle_seconds,
      bootstrap_seconds = bootstrap_seconds,
      total_seconds = total_seconds,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  diagnostics <- plot_sunspots_temporal_uniform_beta_diagnostics(
    s,
    eta_hat,
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_diagnostics.png"
    ),
    control = control
  )

  retained$temporal_fitted_density <-
    sunspots_joint_parsimonious_time_density(
      s,
      eta_hat,
      time_model = "uniform_beta",
      control = control
    )
  retained$temporal_pit <- diagnostics$pit

  utils::write.csv(
    summary,
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_gof_results.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    bootstrap,
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_bootstrap.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    retained,
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_retained_data.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    eta_hat$boundary_diagnostics,
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_boundary_diagnostics.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      temporal_mle_seconds = mle_seconds,
      bootstrap_seconds = bootstrap_seconds,
      total_seconds = total_seconds
    ),
    file.path(
      output_dir,
      "cycle23_temporal_uniform_beta_timing.csv"
    ),
    row.names = FALSE
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
  )

  invisible(list(
    summary = summary,
    fit = eta_hat,
    observed = observed,
    bootstrap = bootstrap,
    retained = retained,
    output_dir = output_dir
  ))
}

parse_sunspots_temporal_uniform_beta_args <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())

  out <- list()
  integer_keys <- c(
    "B",
    "n_cores",
    "dequantization_seed",
    "bootstrap_seed"
  )
  character_keys <- c(
    "input_csv",
    "output_dir",
    "start_date",
    "end_date"
  )

  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --statistics=ks,cvm ",
        "--B=INTEGER --n_cores=INTEGER ",
        "--dequantization_seed=INTEGER ",
        "--bootstrap_seed=INTEGER ",
        "--start_date=YYYY-MM-DD ",
        "--end_date=YYYY-MM-DD ",
        "--input_csv=PATH --output_dir=PATH\n"
      ))
      quit(save = "no", status = 0L)
    }

    parts <- strsplit(
      sub("^--", "", arg),
      "=",
      fixed = TRUE
    )[[1L]]
    if (length(parts) != 2L) {
      stop(sprintf("Invalid option: %s", arg))
    }
    key <- parts[[1L]]
    value <- parts[[2L]]

    if (key %in% integer_keys) {
      out[[key]] <- as.integer(value)
    } else if (identical(key, "statistics")) {
      out[[key]] <- strsplit(
        tolower(value),
        ",",
        fixed = TRUE
      )[[1L]]
    } else if (key %in% character_keys) {
      out[[key]] <- value
    } else {
      stop(sprintf("Unknown option: --%s", key))
    }
  }

  out
}

if (sys.nframe() == 0L) {
  do.call(
    run_sunspots_cycle23_temporal_uniform_beta_gof,
    parse_sunspots_temporal_uniform_beta_args()
  )
}
