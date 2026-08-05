#!/usr/bin/env Rscript

# Composite-null GOF for the temporal two-beta mixture on cycle-23 sunspots.
# The null is composite: each bootstrap sample is re-fitted before computing
# KS, CvM, and AD statistics.

resolve_sunspots_temporal_gof_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) return(candidate)
  }
  stop(sprintf("Could not resolve path: %s", file.path(...)))
}

source(resolve_sunspots_temporal_gof_path(
  "real_data", "sunspots", "sunspots_cycle23_joint_time_space.R"
))

sunspots_temporal_gof_with_seed <- function(seed, expr) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || !is.finite(seed)) stop("`seed` must be one finite integer.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  eval.parent(substitute(expr))
}

sunspots_temporal_cycle23_full_window <- function(input_csv) {
  if (!file.exists(input_csv)) stop(sprintf("Input CSV not found: %s", input_csv))
  input <- utils::read.csv(input_csv, stringsAsFactors = FALSE)
  if (!all(c("cycle", "date") %in% names(input))) {
    stop("`input_csv` must contain `cycle` and `date` columns.")
  }
  if (any(input$cycle != 23L)) stop("This runner is restricted to cycle 23.")
  calendar_day <- as.Date(as.POSIXct(input$date, tz = "UTC"), tz = "UTC")
  if (anyNA(calendar_day)) stop("The input `date` column contains invalid timestamps.")
  list(
    start_date = as.character(min(calendar_day)),
    end_date = as.character(max(calendar_day) + 1)
  )
}

sunspots_temporal_u_values <- function(s, eta_hat) {
  u <- sunspots_joint_time_cdf(s, eta_hat)
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}

sunspots_temporal_eta_from_vector <- function(param, control = list()) {
  param <- as.numeric(param)
  if (length(param) != 5L || any(!is.finite(param))) {
    stop("Temporal beta-mixture parameters must contain five finite entries.")
  }
  sunspots_joint_time_canonicalize_eta(list(
    weight1 = param[[1L]],
    alpha1 = param[[2L]],
    beta1 = param[[3L]],
    alpha2 = param[[4L]],
    beta2 = param[[5L]]
  ), control = control)
}

sunspots_temporal_require_rgof <- function() {
  if (!requireNamespace("Rgof", quietly = TRUE)) {
    stop(
      paste(
        "The Rgof package is required for `engine = 'rgof'`.",
        "Install it with install.packages('Rgof')."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sunspots_temporal_rgof_composite_refit_audit <- function(B = 50L,
                                                          sample_size = 80L,
                                                          seed = 20260805L,
                                                          gof_test_fun = NULL) {
  B <- as.integer(B)
  sample_size <- as.integer(sample_size)
  if (!is.finite(B) || B < 1L) stop("`B` must be a positive integer.")
  if (!is.finite(sample_size) || sample_size < 10L) stop("`sample_size` must be an integer >= 10.")
  if (is.null(gof_test_fun)) {
    sunspots_temporal_require_rgof()
    gof_test_fun <- Rgof::gof_test
  }

  phat_calls <- 0L
  phat_means <- numeric(0L)
  phat <- function(x) {
    phat_calls <<- phat_calls + 1L
    phat_means <<- c(phat_means, mean(x))
    c(mean(x), stats::sd(x))
  }
  pnull <- function(x, param) {
    sigma <- max(as.numeric(param[[2L]]), 1e-8)
    stats::pnorm(x, mean = as.numeric(param[[1L]]), sd = sigma)
  }
  rnull <- function(param) {
    sigma <- max(as.numeric(param[[2L]]), 1e-8)
    stats::rnorm(sample_size, mean = as.numeric(param[[1L]]), sd = sigma)
  }

  x <- sunspots_temporal_gof_with_seed(seed, stats::rnorm(sample_size, mean = 0.1, sd = 1.3))
  fit <- sunspots_temporal_gof_with_seed(seed + 1L, gof_test_fun(
    x = x,
    vals = NA,
    pnull = pnull,
    rnull = rnull,
    phat = phat,
    B = B,
    doMethods = c("KS", "CvM", "AD"),
    maxProcessor = 1L
  ))

  unique_means <- length(unique(signif(phat_means, 10L)))
  required_min_calls <- B + 1L
  list(
    B = B,
    sample_size = sample_size,
    phat_calls = phat_calls,
    required_min_calls = required_min_calls,
    unique_phat_means = unique_means,
    refit_pass = isTRUE(phat_calls >= required_min_calls && unique_means > 1L),
    statistics_names = names(fit$statistics %||% numeric(0L)),
    pvalue_names = names(fit$p.values %||% numeric(0L))
  )
}

sunspots_temporal_rgof_gof <- function(s, B, control = list(), n_cores = 1L,
                                       audit_B = 50L, audit_seed = 20260805L) {
  sunspots_temporal_require_rgof()
  audit <- sunspots_temporal_rgof_composite_refit_audit(
    B = as.integer(audit_B),
    sample_size = max(40L, min(200L, length(s))),
    seed = as.integer(audit_seed)
  )
  if (!isTRUE(audit$refit_pass)) {
    stop(
      sprintf(
        paste(
          "Rgof composite-refit audit failed:",
          "phat_calls=%d, required_min_calls=%d, unique_phat_means=%d.",
          "The run is stopped to avoid reporting a non-composite result."
        ),
        audit$phat_calls,
        audit$required_min_calls,
        audit$unique_phat_means
      ),
      call. = FALSE
    )
  }

  phat <- function(x) {
    fit <- fit_sunspots_joint_time_beta_mixture2(x, control = control)
    c(fit$weight1, fit$alpha1, fit$beta1, fit$alpha2, fit$beta2)
  }
  pnull <- function(x, param) {
    eta <- sunspots_temporal_eta_from_vector(param, control = control)
    sunspots_joint_time_cdf(x, eta, control = control)
  }
  rnull <- function(param) {
    eta <- sunspots_temporal_eta_from_vector(param, control = control)
    sample_sunspots_joint_time_beta_mixture2(length(s), eta, control = control)
  }

  fit <- Rgof::gof_test(
    x = s,
    vals = NA,
    pnull = pnull,
    rnull = rnull,
    phat = phat,
    doMethods = c("KS", "CvM", "AD"),
    B = as.integer(B),
    maxProcessor = as.integer(max(1L, n_cores))
  )
  list(result = fit, audit = audit)
}

sunspots_temporal_gof_statistics_from_u <- function(u) {
  u <- sort(as.numeric(u))
  n <- length(u)
  if (n < 1L) stop("`u` must contain at least one value.")
  i <- seq_len(n)

  ks <- max(c(i / n - u, u - (i - 1) / n))
  cvm <- 1 / (12 * n) + sum((u - (2 * i - 1) / (2 * n))^2)
  ad <- -n - mean((2 * i - 1) * (log(u) + log(1 - rev(u))))

  c(ks = ks, cvm = cvm, ad = ad)
}

sunspots_temporal_fit_row <- function(eta_hat, n, fit_type) {
  criteria <- sunspots_joint_time_information_criteria(eta_hat$loglik, n = n, n_parameters = 5L)
  data.frame(
    fit_type = fit_type,
    n = as.integer(n),
    weight1 = eta_hat$weight1,
    alpha1 = eta_hat$alpha1,
    beta1 = eta_hat$beta1,
    alpha2 = eta_hat$alpha2,
    beta2 = eta_hat$beta2,
    mean1 = eta_hat$mean1,
    mean2 = eta_hat$mean2,
    loglik = eta_hat$loglik,
    aic = criteria$aic,
    bic = criteria$bic,
    convergence = as.integer(eta_hat$opt$convergence %||% NA_integer_),
    convergence_message = as.character(eta_hat$opt$message %||% ""),
    n_starts = as.integer(eta_hat$n_starts),
    n_finite_fits = as.integer(eta_hat$n_successful_starts),
    n_converged_fits = as.integer(eta_hat$n_converged_starts),
    selected_converged = isTRUE(eta_hat$selected_converged),
    any_boundary = any(unlist(eta_hat$boundary_flags, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

run_sunspots_cycle23_temporal_beta_composite_gof <- function(
    input_csv = file.path("real_data", "sunspots", "output", "sunspots_cycle23_s2_all.csv"),
    output_dir = file.path("real_data", "sunspots", "output", "cycle23_temporal_beta_composite_gof_full"),
    dequantization_seed = 20260712L,
    bootstrap_seed = 20260805L,
    engine = c("manual", "rgof"),
    B = 999L,
    n_cores = 1L,
    rgof_audit_B = 50L,
    rgof_audit_seed = 20260805L,
    control = list()) {
  engine <- tolower(as.character(engine[[1L]]))
  if (!engine %in% c("manual", "rgof")) {
    stop("`engine` must be either 'manual' or 'rgof'.")
  }
  window <- sunspots_temporal_cycle23_full_window(input_csv)
  data <- prepare_sunspots_joint_time_data(
    input_csv = input_csv,
    start_date = window$start_date,
    end_date = window$end_date,
    dequantization_seed = dequantization_seed
  )
  n <- nrow(data)
  if (n < 20L) stop("Too few observations for a stable composite-null GOF run.")

  eta_hat <- fit_sunspots_joint_time_beta_mixture2(data$s, control = control)
  u_obs <- sunspots_temporal_u_values(data$s, eta_hat)
  stats_obs <- sunspots_temporal_gof_statistics_from_u(u_obs)

  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("`B` must be a positive integer.")

  bootstrap <- sunspots_temporal_gof_with_seed(bootstrap_seed, lapply(seq_len(B), function(b) {
    s_boot <- sample_sunspots_joint_time_beta_mixture2(n, eta_hat, control = control)
    fit_boot <- fit_sunspots_joint_time_beta_mixture2(s_boot, control = control)
    u_boot <- sunspots_temporal_u_values(s_boot, fit_boot)
    stats_boot <- sunspots_temporal_gof_statistics_from_u(u_boot)
    c(
      replicate = b,
      ks = stats_boot[["ks"]],
      cvm = stats_boot[["cvm"]],
      ad = stats_boot[["ad"]],
      convergence = as.integer(fit_boot$opt$convergence %||% NA_integer_),
      n_converged_fits = as.integer(fit_boot$n_converged_starts)
    )
  }))
  bootstrap <- do.call(rbind, bootstrap)

  p_values <- c(
    ks = (1 + sum(bootstrap[, "ks"] >= stats_obs[["ks"]])) / (B + 1),
    cvm = (1 + sum(bootstrap[, "cvm"] >= stats_obs[["cvm"]])) / (B + 1),
    ad = (1 + sum(bootstrap[, "ad"] >= stats_obs[["ad"]])) / (B + 1)
  )

  rgof_observed <- c(ks = NA_real_, cvm = NA_real_, ad = NA_real_)
  rgof_p_values <- c(ks = NA_real_, cvm = NA_real_, ad = NA_real_)
  rgof_audit <- NULL
  if (identical(engine, "rgof")) {
    rgof_out <- sunspots_temporal_gof_with_seed(bootstrap_seed, sunspots_temporal_rgof_gof(
      s = data$s,
      B = B,
      control = control,
      n_cores = as.integer(n_cores),
      audit_B = as.integer(rgof_audit_B),
      audit_seed = as.integer(rgof_audit_seed)
    ))
    rgof_audit <- rgof_out$audit
    rgof_observed <- c(
      ks = as.numeric(rgof_out$result$statistics[["KS"]] %||% NA_real_),
      cvm = as.numeric(rgof_out$result$statistics[["CvM"]] %||% NA_real_),
      ad = as.numeric(rgof_out$result$statistics[["AD"]] %||% NA_real_)
    )
    rgof_p_values <- c(
      ks = as.numeric(rgof_out$result$p.values[["KS"]] %||% NA_real_),
      cvm = as.numeric(rgof_out$result$p.values[["CvM"]] %||% NA_real_),
      ad = as.numeric(rgof_out$result$p.values[["AD"]] %||% NA_real_)
    )
  }

  summary <- data.frame(
    sample = "full_cycle",
    n = n,
    start_date = window$start_date,
    end_date_exclusive = window$end_date,
    dequantization_seed = as.integer(dequantization_seed),
    bootstrap_seed = as.integer(bootstrap_seed),
    engine = engine,
    B = B,
    ks_observed = stats_obs[["ks"]],
    cvm_observed = stats_obs[["cvm"]],
    ad_observed = stats_obs[["ad"]],
    ks_p_value = p_values[["ks"]],
    cvm_p_value = p_values[["cvm"]],
    ad_p_value = p_values[["ad"]],
    rgof_ks_observed = rgof_observed[["ks"]],
    rgof_cvm_observed = rgof_observed[["cvm"]],
    rgof_ad_observed = rgof_observed[["ad"]],
    rgof_ks_p_value = rgof_p_values[["ks"]],
    rgof_cvm_p_value = rgof_p_values[["cvm"]],
    rgof_ad_p_value = rgof_p_values[["ad"]],
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fit_row <- sunspots_temporal_fit_row(eta_hat, n, fit_type = "observed")
  utils::write.csv(fit_row, file.path(output_dir, "temporal_beta_observed_fit.csv"), row.names = FALSE)
  utils::write.csv(eta_hat$boundary_diagnostics, file.path(output_dir, "temporal_beta_observed_boundary_diagnostics.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "temporal_beta_composite_gof_summary.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(bootstrap), file.path(output_dir, "temporal_beta_composite_gof_bootstrap_statistics.csv"), row.names = FALSE)
  utils::write.csv(data, file.path(output_dir, "temporal_beta_composite_gof_retained_data.csv"), row.names = FALSE)
  if (!is.null(rgof_audit)) {
    utils::write.csv(data.frame(
      B = rgof_audit$B,
      sample_size = rgof_audit$sample_size,
      phat_calls = rgof_audit$phat_calls,
      required_min_calls = rgof_audit$required_min_calls,
      unique_phat_means = rgof_audit$unique_phat_means,
      refit_pass = rgof_audit$refit_pass,
      statistics_names = paste(rgof_audit$statistics_names, collapse = ","),
      pvalue_names = paste(rgof_audit$pvalue_names, collapse = ","),
      stringsAsFactors = FALSE
    ), file.path(output_dir, "temporal_beta_composite_gof_rgof_audit.csv"), row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  message("Composite-null GOF completed.")
  message(sprintf("Output directory: %s", output_dir))
  invisible(list(summary = summary, observed_fit = fit_row, output_dir = output_dir))
}

parse_sunspots_temporal_beta_composite_gof_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) return(list())
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(paste0(
        "Options: --input_csv=PATH --output_dir=PATH --dequantization_seed=INTEGER ",
        "--bootstrap_seed=INTEGER --engine=manual|rgof --B=INTEGER --n_cores=INTEGER ",
        "--rgof_audit_B=INTEGER --rgof_audit_seed=INTEGER\n"
      ))
      quit(save = "no", status = 0L)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) stop(sprintf("Invalid option: %s", arg))
    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key %in% c("input_csv", "output_dir")) out[[key]] <- value
    if (identical(key, "engine")) out[[key]] <- tolower(value)
    if (key %in% c("dequantization_seed", "bootstrap_seed", "B", "n_cores", "rgof_audit_B", "rgof_audit_seed")) {
      out[[key]] <- as.integer(value)
    }
  }
  out
}

if (sys.nframe() == 0L) {
  do.call(run_sunspots_cycle23_temporal_beta_composite_gof, parse_sunspots_temporal_beta_composite_gof_args())
}
