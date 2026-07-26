#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory diagnostic script.  It does not establish
# a production bootstrap or quadrature rule, and its conclusions require
# mathematical review and explicit approval before any use in the paper.

# Audit the stored multiplier-bootstrap distributions behind selected real-data
# KS/CvM p-values.  This script does not modify the manuscript.
#
# It uses the exact PogoJump result that feeds the paper table, contrasts it
# with the stored fully re-estimated run, and reruns the paper's Risoe
# November--December 125m case while retaining every bootstrap replicate.

`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

parse_arguments <- function(args) {
  out <- list(
    output_dir = file.path("real_data", "bootstrap_audit", "ks_cvm_real_data_20260723"),
    wind_B = 5000L,
    wind_seed = 2026072301L,
    wind_cores = 2L,
    rerun_wind = TRUE,
    matched_B = 5000L,
    matched_seed = 2026072401L,
    matched_cores = 8L,
    rerun_matched = FALSE
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    if (identical(key, "output_dir")) out$output_dir <- value
    if (identical(key, "wind_B")) out$wind_B <- as.integer(value)
    if (identical(key, "wind_seed")) out$wind_seed <- as.integer(value)
    if (identical(key, "wind_cores")) out$wind_cores <- as.integer(value)
    if (identical(key, "rerun_wind")) out$rerun_wind <- tolower(value) %in% c("true", "t", "1", "yes")
    if (identical(key, "matched_B")) out$matched_B <- as.integer(value)
    if (identical(key, "matched_seed")) out$matched_seed <- as.integer(value)
    if (identical(key, "matched_cores")) out$matched_cores <- as.integer(value)
    if (identical(key, "rerun_matched")) out$rerun_matched <- tolower(value) %in% c("true", "t", "1", "yes")
  }

  if (!is.finite(out$wind_B) || out$wind_B < 1L) stop("`wind_B` must be a positive integer.")
  if (!is.finite(out$wind_seed)) stop("`wind_seed` must be finite.")
  if (!is.finite(out$wind_cores) || out$wind_cores < 1L) stop("`wind_cores` must be a positive integer.")
  if (!is.finite(out$matched_B) || out$matched_B < 1L) stop("`matched_B` must be a positive integer.")
  if (!is.finite(out$matched_seed)) stop("`matched_seed` must be finite.")
  if (!is.finite(out$matched_cores) || out$matched_cores < 1L) stop("`matched_cores` must be a positive integer.")
  out
}

wilson_interval <- function(k, B, level = 0.95) {
  if (B <= 0L) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- k / B
  den <- 1 + z^2 / B
  center <- (p + z^2 / (2 * B)) / den
  radius <- z * sqrt((p * (1 - p) + z^2 / (4 * B)) / B) / den
  c(lower = max(0, center - radius), upper = min(1, center + radius))
}

safe_lag1_cor <- function(x) {
  if (length(x) < 3L || stats::sd(x) == 0) return(NA_real_)
  stats::cor(x[-length(x)], x[-1L])
}

bootstrap_stat_summary <- function(case_id,
                                   domain,
                                   method,
                                   statistic,
                                   observed,
                                   values,
                                   stored_p_value = NA_real_,
                                   manuscript_p_value = NA_real_) {
  values <- as.numeric(values)
  B <- length(values)
  finite <- is.finite(values)
  values_finite <- values[finite]
  exceedances <- sum(values_finite >= observed)
  p_value <- (1 + exceedances) / (B + 1)
  p_raw <- exceedances / B
  mc_se <- sqrt(p_raw * (1 - p_raw) / B)
  ci <- wilson_interval(exceedances, B)
  qs <- stats::quantile(values_finite, c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99), type = 8, names = FALSE)

  n_batches <- min(20L, B)
  batch_id <- cut(seq_len(B), breaks = n_batches, labels = FALSE)
  batch_p <- vapply(split(values_finite, batch_id), function(v) {
    (1 + sum(v >= observed)) / (1 + length(v))
  }, numeric(1))

  data.frame(
    case = case_id,
    domain = domain,
    method = method,
    statistic = statistic,
    B = B,
    observed = observed,
    finite_replicates = sum(finite),
    nonfinite_replicates = sum(!finite),
    negative_replicates = sum(values_finite < 0),
    exceedances = exceedances,
    p_value = p_value,
    stored_p_value = stored_p_value,
    manuscript_p_value = manuscript_p_value,
    p_reconstruction_error = if (is.finite(stored_p_value)) p_value - stored_p_value else NA_real_,
    mc_se = mc_se,
    mc_wilson_lower = ci[["lower"]],
    mc_wilson_upper = ci[["upper"]],
    bootstrap_mean = mean(values_finite),
    bootstrap_sd = stats::sd(values_finite),
    bootstrap_cv = stats::sd(values_finite) / mean(values_finite),
    bootstrap_min = min(values_finite),
    q01 = qs[[1L]], q05 = qs[[2L]], q25 = qs[[3L]], q50 = qs[[4L]],
    q75 = qs[[5L]], q95 = qs[[6L]], q99 = qs[[7L]],
    bootstrap_max = max(values_finite),
    observed_z = (observed - mean(values_finite)) / stats::sd(values_finite),
    observed_ecdf = mean(values_finite <= observed),
    lag1_correlation = safe_lag1_cor(values_finite),
    batch_p_mean = mean(batch_p),
    batch_p_sd = stats::sd(batch_p),
    batch_p_min = min(batch_p),
    batch_p_max = max(batch_p),
    stringsAsFactors = FALSE
  )
}

extract_bootstrap_payload <- function(result, fallback_method = NA_character_) {
  raw <- result$bootstrap$raw_result %||% result
  values <- raw$bootstrap$statistics
  observed <- vapply(c("ks", "cvm"), function(statistic) {
    raw$observed[[statistic]]$statistic
  }, numeric(1))
  p_values <- vapply(c("ks", "cvm"), function(statistic) {
    raw$inference[[statistic]]$p_value
  }, numeric(1))
  list(
    raw = raw,
    ks = as.numeric(values$ks),
    cvm = as.numeric(values$cvm),
    observed = observed,
    p_values = p_values,
    method = raw$diagnostics$effective_bootstrap_method %||% fallback_method,
    diagnostics = raw$diagnostics
  )
}

save_case_replicates <- function(payload, case_id, output_dir) {
  replicates <- data.frame(
    replicate = seq_along(payload$ks),
    ks = payload$ks,
    cvm = payload$cvm
  )
  utils::write.csv(
    replicates,
    file.path(output_dir, paste0(case_id, "_bootstrap_replicates.csv")),
    row.names = FALSE
  )
  saveRDS(replicates, file.path(output_dir, paste0(case_id, "_bootstrap_replicates.rds")))
  replicates
}

plot_case_bootstrap <- function(payload, case_id, output_dir) {
  path <- file.path(output_dir, "plots", paste0(case_id, "_bootstrap_distributions.png"))
  grDevices::png(path, width = 2100, height = 720, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3.2, 1))

  graphics::hist(payload$ks, breaks = "FD", main = paste0(case_id, ": KS"),
                 xlab = "Bootstrap KS statistic", col = "grey85", border = "white")
  graphics::abline(v = payload$observed[["ks"]], col = "firebrick3", lwd = 3)

  graphics::hist(payload$cvm, breaks = "FD", main = paste0(case_id, ": CvM"),
                 xlab = "Bootstrap CvM statistic", col = "grey85", border = "white")
  graphics::abline(v = payload$observed[["cvm"]], col = "firebrick3", lwd = 3)

  graphics::plot(payload$ks, payload$cvm, pch = 16, cex = 0.35,
                 col = grDevices::adjustcolor("grey25", alpha.f = 0.30),
                 xlab = "Bootstrap KS statistic", ylab = "Bootstrap CvM statistic",
                 main = paste0(case_id, ": paired replicates"))
  graphics::points(payload$observed[["ks"]], payload$observed[["cvm"]],
                   pch = 8, cex = 2, lwd = 2, col = "firebrick3")
  invisible(path)
}

run_wind_fast_bootstrap <- function(B, seed, n_cores) {
  source(file.path("bootstrap", "model_specs.R"))
  source(file.path("bootstrap", "multiplier_bootstrap.R"))
  data_path <- file.path(
    "real_data", "wind", "month_diagnostics", "windows_contiguous_two_month_step4_b1000_seed20260528",
    "risoe_clean_125m_nov_dec_start4_1996_2003_hvmf.csv"
  )
  if (!file.exists(data_path)) stop("The stored Risoe paper-case CSV was not found: ", data_path)
  df_case <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  config <- list(
    dataset_id = "risoe_paper_nov_dec_125m_start4",
    window = "nov_dec",
    months = c(11L, 12L),
    pattern = "start4",
    day_pattern = seq.int(4L, 30L, by = 4L)
  )
  if (nrow(df_case) != 101L) stop("The stored Risoe paper case must contain 101 observations.")
  X <- as.matrix(df_case[, c("x0", "x1", "x2")])
  result <- multiplier_bootstrap_hvmf(
    data = X,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 1000L,
      derivative_mc_seed = 20260613L,
      fast_multiplier_cvm_block_size = 50L
    ),
    unknown_param = "both"
  )
  list(
    result = result,
    n = nrow(df_case),
    data_path = data_path,
    config = config
  )
}

audit_logistic_gaussian_hbe_profiles <- function(
    results_dir,
    selected_slugs = c("pogojump", "coxite", "animalvegetation"),
    output_dir = NULL) {
  # The paper screening runs explicitly request the HBE approximation.  This
  # audit evaluates the same fitted profiles with Farebrother's numerical CDF
  # instead, and compares the resulting observed statistics against the
  # *stored* fast-multiplier replicates.  The latter are based on indicators
  # and score Monte Carlo, rather than on the HBE profile evaluation.
  source(file.path("bootstrap", "model_specs.R"))
  result_files <- list.files(results_dir, pattern = "_results[.]rds$", full.names = TRUE)
  selected_files <- paste0(selected_slugs, "_results.rds")
  result_files <- result_files[basename(result_files) %in% selected_files]

  rows <- lapply(result_files, function(path) {
    result <- readRDS(path)
    raw <- result$bootstrap$raw_result
    if (is.null(raw$observed$cvm$process_matrix) ||
        is.null(raw$bootstrap$statistics$ks) ||
        is.null(result$data_prep$X_closed)) {
      return(NULL)
    }

    observed <- raw$observed$cvm
    theta_hat <- raw$observed$theta_hat
    n <- nrow(observed$process_matrix)
    normalized <- normalize_logistic_gaussian_data(result$data_prep$X_closed)
    shift_matrix <- matrix(
      rep(theta_hat$mu_ilr, each = n),
      nrow = n,
      ncol = length(theta_hat$mu_ilr)
    ) - normalized$ilr

    exact_profile <- tryCatch(
      suppressWarnings(evaluate_mvnorm_distance_profile_matrix(
        shift_matrix = shift_matrix,
        t_matrix = observed$distance_matrix,
        eigenvalues_full = theta_hat$eigenvalues_full,
        eigenvectors_full = theta_hat$eigenvectors_full,
        positive_idx = theta_hat$positive_idx,
        control = list(logistic_gaussian_quadform_method = "farebrother")
      )),
      error = function(e) e
    )
    if (inherits(exact_profile, "error")) {
      return(data.frame(
        dataset = result$dataset_name,
        n = n,
        status = paste("farebrother_error:", conditionMessage(exact_profile)),
        stringsAsFactors = FALSE
      ))
    }

    exact_process <- sqrt(n) * (observed$empirical_profile - exact_profile)
    exact_ks <- max(abs(exact_process))
    exact_cvm <- mean(exact_process^2)
    boot_ks <- as.numeric(raw$bootstrap$statistics$ks)
    boot_cvm <- as.numeric(raw$bootstrap$statistics$cvm)
    B <- length(boot_ks)
    profile_error <- abs(observed$theoretical_profile - exact_profile)
    max_index <- arrayInd(which.max(profile_error), dim(profile_error))

    if (!is.null(output_dir)) {
      slug <- tolower(gsub("[^a-zA-Z0-9]+", "_", result$dataset_name))
      pointwise <- data.frame(
        center_index = rep.int(seq_len(n), n),
        threshold_index = rep(seq_len(n), each = n),
        distance = as.vector(observed$distance_matrix),
        empirical_profile = as.vector(observed$empirical_profile),
        hbe_profile = as.vector(observed$theoretical_profile),
        farebrother_profile = as.vector(exact_profile),
        hbe_process = as.vector(observed$process_matrix),
        farebrother_process = as.vector(exact_process),
        abs_profile_difference = as.vector(profile_error)
      )
      utils::write.csv(
        pointwise,
        file.path(output_dir, paste0(slug, "_hbe_vs_farebrother_distance_profiles.csv")),
        row.names = FALSE
      )

      plot_path <- file.path(output_dir, "plots", paste0(slug, "_hbe_vs_farebrother_profiles.png"))
      grDevices::png(plot_path, width = 1800, height = 720, res = 150)
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
      graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1))
      graphics::plot(
        pointwise$farebrother_profile, pointwise$hbe_profile,
        pch = 16, cex = 0.45, col = grDevices::adjustcolor("grey25", alpha.f = 0.35),
        xlab = "Farebrother profile", ylab = "HBE profile",
        main = paste0(result$dataset_name, ": theoretical profiles")
      )
      graphics::abline(0, 1, col = "firebrick3", lwd = 2)
      graphics::plot(
        pointwise$farebrother_process, pointwise$hbe_process,
        pch = 16, cex = 0.45, col = grDevices::adjustcolor("grey25", alpha.f = 0.35),
        xlab = "Farebrother empirical process", ylab = "HBE empirical process",
        main = paste0(result$dataset_name, ": observed process")
      )
      graphics::abline(0, 1, col = "firebrick3", lwd = 2)
    }

    data.frame(
      dataset = result$dataset_name,
      n = n,
      status = "ok",
      hbe_ks = raw$observed$ks$statistic,
      hbe_cvm = raw$observed$cvm$statistic,
      hbe_p_ks = raw$inference$ks$p_value,
      hbe_p_cvm = raw$inference$cvm$p_value,
      farebrother_ks = exact_ks,
      farebrother_cvm = exact_cvm,
      farebrother_p_ks_against_stored_fast_bootstrap = (1 + sum(boot_ks >= exact_ks)) / (B + 1),
      farebrother_p_cvm_against_stored_fast_bootstrap = (1 + sum(boot_cvm >= exact_cvm)) / (B + 1),
      max_abs_profile_difference = max(profile_error),
      mean_abs_profile_difference = mean(profile_error),
      q95_abs_profile_difference = as.numeric(stats::quantile(profile_error, 0.95, names = FALSE)),
      fraction_hbe_zero_farebrother_gt_001 = mean(observed$theoretical_profile == 0 & exact_profile > 0.01),
      max_difference_center_index = max_index[[1L]],
      max_difference_threshold_index = max_index[[2L]],
      max_difference_distance = observed$distance_matrix[max_index[[1L]], max_index[[2L]]],
      hbe_profile_at_max_difference = observed$theoretical_profile[max_index[[1L]], max_index[[2L]]],
      farebrother_profile_at_max_difference = exact_profile[max_index[[1L]], max_index[[2L]]],
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

audit_hvmf_tabulated_profile <- function(data_path, raw_result, output_dir = NULL) {
  # The wind code uses a spline-tabulated exact profile for CvM.  Verify its
  # numerical error directly against the integration implementation.
  source("utils.R")
  df_case <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  data <- as.matrix(df_case[, c("x0", "x1", "x2")])
  observed <- raw_result$observed$cvm
  theta_hat <- raw_result$observed$theta_hat
  n <- nrow(data)
  exact_profile <- t(vapply(seq_len(n), function(i) {
    theoretical_distance_profile_hvmf(
      omega = data[i, ], mu = theta_hat$mu, kappa = theta_hat$kappa,
      t_values = observed$distance_matrix[i, ]
    )
  }, numeric(n)))
  profile_error <- abs(observed$theoretical_profile - exact_profile)
  exact_process <- sqrt(n) * (observed$empirical_profile - exact_profile)
  if (!is.null(output_dir)) {
    pointwise <- data.frame(
      center_index = rep.int(seq_len(n), n),
      threshold_index = rep(seq_len(n), each = n),
      distance = as.vector(observed$distance_matrix),
      empirical_profile = as.vector(observed$empirical_profile),
      tabulated_profile = as.vector(observed$theoretical_profile),
      exact_profile = as.vector(exact_profile),
      tabulated_process = as.vector(observed$process_matrix),
      exact_process = as.vector(exact_process),
      abs_profile_difference = as.vector(profile_error)
    )
    utils::write.csv(
      pointwise,
      file.path(output_dir, "risoe_hvmf_tabulated_vs_exact_distance_profiles.csv"),
      row.names = FALSE
    )
  }

  data.frame(
    case = "risoe_nov_dec_125m_start4",
    n = n,
    tabulated_ks = raw_result$observed$ks$statistic,
    tabulated_cvm = raw_result$observed$cvm$statistic,
    exact_ks = max(abs(exact_process)),
    exact_cvm = mean(exact_process^2),
    max_abs_profile_difference = max(profile_error),
    mean_abs_profile_difference = mean(profile_error),
    q95_abs_profile_difference = as.numeric(stats::quantile(profile_error, 0.95, names = FALSE)),
    stringsAsFactors = FALSE
  )
}

make_fixed_exponential_multiplier_spec <- function(raw_matrix) {
  raw_matrix <- as.matrix(raw_matrix)
  cursor <- 0L
  list(
    name = "fixed_Exp1_matrix",
    mean = 1,
    sd = 1,
    generator = function(n) {
      cursor <<- cursor + 1L
      if (cursor > nrow(raw_matrix)) stop("Fixed multiplier matrix was exhausted.")
      if (as.integer(n) != ncol(raw_matrix)) stop("Fixed multiplier row has incompatible length.")
      as.numeric(raw_matrix[cursor, ])
    }
  )
}

make_raw_exponential_multiplier_matrix <- function(B, n, seed) {
  set.seed(as.integer(seed))
  matrix(stats::rexp(as.integer(B) * as.integer(n)), nrow = as.integer(B), ncol = as.integer(n), byrow = TRUE)
}

tail_probability <- function(values, observed) {
  values <- as.numeric(values)
  (1 + sum(values >= observed)) / (length(values) + 1)
}

run_matched_bootstrap_pair <- function(wrapper,
                                       data,
                                       control,
                                       B,
                                       seed,
                                       n_cores,
                                       wrapper_label) {
  n <- nrow(data)
  raw_multipliers <- make_raw_exponential_multiplier_matrix(B = B, n = n, seed = seed)
  common_args <- list(
    data = data,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = make_sample_unique_distance_ks_grid(),
    B = as.integer(B),
    alpha = 0.05,
    n_cores = as.integer(n_cores),
    seed = NULL,
    keep = list(
      observed_process = TRUE,
      bootstrap_statistics = TRUE,
      bootstrap_thetas = FALSE
    ),
    control = control,
    unknown_param = "both",
    distance_profile_backend = "r"
  )
  fast <- do.call(wrapper, c(
    common_args,
    list(
      multipliers = make_fixed_exponential_multiplier_spec(raw_multipliers),
      bootstrap_method = "fast_multiplier"
    )
  ))
  reestimated <- do.call(wrapper, c(
    common_args,
    list(
      multipliers = make_fixed_exponential_multiplier_spec(raw_multipliers),
      bootstrap_method = "reestimated"
    )
  ))

  observed_fast <- vapply(c("ks", "cvm"), function(name) fast$observed[[name]]$statistic, numeric(1))
  observed_reestimated <- vapply(c("ks", "cvm"), function(name) reestimated$observed[[name]]$statistic, numeric(1))
  if (!isTRUE(all.equal(observed_fast, observed_reestimated, tolerance = 1e-12))) {
    stop("Fast and re-estimated runs disagree on the observed statistics despite identical data and code.")
  }

  list(
    label = wrapper_label,
    B = as.integer(B),
    seed = as.integer(seed),
    n_cores = as.integer(n_cores),
    raw_multipliers = raw_multipliers,
    fast = fast,
    reestimated = reestimated
  )
}

summarize_matched_bootstrap_pair <- function(pair, case_id, domain) {
  stats <- c("ks", "cvm")
  do.call(rbind, lapply(stats, function(statistic) {
    fast_values <- as.numeric(pair$fast$bootstrap$statistics[[statistic]])
    reestimated_values <- as.numeric(pair$reestimated$bootstrap$statistics[[statistic]])
    observed <- pair$fast$observed[[statistic]]$statistic
    data.frame(
      case = case_id,
      domain = domain,
      statistic = toupper(statistic),
      B = length(fast_values),
      observed = observed,
      p_fast = tail_probability(fast_values, observed),
      p_reestimated = tail_probability(reestimated_values, observed),
      paired_correlation = stats::cor(fast_values, reestimated_values),
      mean_fast = mean(fast_values),
      mean_reestimated = mean(reestimated_values),
      mean_difference = mean(fast_values - reestimated_values),
      sd_difference = stats::sd(fast_values - reestimated_values),
      q95_abs_difference = as.numeric(stats::quantile(abs(fast_values - reestimated_values), 0.95, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
}

matched_pair_replicates <- function(pair) {
  data.frame(
    replicate = seq_len(pair$B),
    fast_ks = as.numeric(pair$fast$bootstrap$statistics$ks),
    reestimated_ks = as.numeric(pair$reestimated$bootstrap$statistics$ks),
    fast_cvm = as.numeric(pair$fast$bootstrap$statistics$cvm),
    reestimated_cvm = as.numeric(pair$reestimated$bootstrap$statistics$cvm)
  )
}

plot_matched_bootstrap_pair <- function(pair, case_id, output_dir) {
  path <- file.path(output_dir, "plots", paste0(case_id, "_matched_fast_vs_reestimated.png"))
  grDevices::png(path, width = 1500, height = 720, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1))
  for (statistic in c("ks", "cvm")) {
    fast_values <- as.numeric(pair$fast$bootstrap$statistics[[statistic]])
    reestimated_values <- as.numeric(pair$reestimated$bootstrap$statistics[[statistic]])
    limits <- range(c(fast_values, reestimated_values))
    graphics::plot(
      reestimated_values, fast_values,
      pch = 16, cex = 0.35, col = grDevices::adjustcolor("grey25", alpha.f = 0.25),
      xlim = limits, ylim = limits,
      xlab = paste0("Re-estimated ", toupper(statistic)),
      ylab = paste0("Fast ", toupper(statistic)),
      main = paste0(case_id, ": ", toupper(statistic))
    )
    graphics::abline(0, 1, col = "firebrick3", lwd = 2)
  }
  invisible(path)
}

localize_observed_maxima <- function(raw_result, case_id, domain, top_k = 25L) {
  observed <- raw_result$observed$cvm
  process <- observed$process_matrix
  n <- nrow(process)
  top_k <- min(as.integer(top_k), length(process))
  linear_index <- order(abs(process), decreasing = TRUE)[seq_len(top_k)]
  index <- arrayInd(linear_index, dim(process))
  total_cvm_sum <- sum(process^2)
  output <- data.frame(
    case = case_id,
    domain = domain,
    rank = seq_len(top_k),
    center_index = index[, 1L],
    threshold_index = index[, 2L],
    distance = rep(NA_real_, top_k),
    empirical_profile = rep(NA_real_, top_k),
    theoretical_profile = rep(NA_real_, top_k),
    process = rep(NA_real_, top_k),
    abs_process = rep(NA_real_, top_k),
    contribution_to_cvm = rep(NA_real_, top_k),
    cumulative_share_cvm = rep(NA_real_, top_k),
    stringsAsFactors = FALSE
  )
  output$distance <- observed$distance_matrix[linear_index]
  output$empirical_profile <- observed$empirical_profile[linear_index]
  output$theoretical_profile <- observed$theoretical_profile[linear_index]
  output$process <- process[linear_index]
  output$abs_process <- abs(process[linear_index])
  output$contribution_to_cvm <- process[linear_index]^2 / (n * n)
  output$cumulative_share_cvm <- cumsum(process[linear_index]^2) / total_cvm_sum
  output
}

write_report <- function(summary_df, stored_method_comparison, output_dir) {
  lines <- c(
    "# KS/CvM bootstrap audit for selected real-data cases",
    "",
    "The table reports the conditional Monte Carlo uncertainty of each bootstrap p-value. `bootstrap_sd` is the standard deviation of the statistic over bootstrap replicates; `mc_se` is the Monte Carlo standard error of the tail probability estimate, not a sampling standard error for the data-generating model.",
    "",
    "## Bootstrap-statistic summaries",
    ""
  )

  summary_columns <- c("case", "method", "statistic", "B", "observed", "exceedances", "p_value", "mc_se", "mc_wilson_lower", "mc_wilson_upper", "bootstrap_mean", "bootstrap_sd", "observed_z", "lag1_correlation")
  summary_text <- utils::capture.output(print(summary_df[, summary_columns, drop = FALSE], row.names = FALSE, digits = 7))
  lines <- c(lines, "```text", summary_text, "```", "", "## Stored fast versus stored fully re-estimated multiplier bootstrap (PogoJump)", "")
  comparison_text <- utils::capture.output(print(stored_method_comparison, row.names = FALSE, digits = 7))
  lines <- c(lines, "```text", comparison_text, "```", "")
  lines <- c(lines,
    "Interpretation notes:",
    "",
    "- A large difference between the calibrated KS and CvM p-values is not, by itself, an error: KS is a supremum functional, while CvM integrates squared departures. Their null distributions are different, so their p-values admit no general ordering.",
    "- The `p_reconstruction_error` column must be zero up to floating-point precision whenever the stored raw replicates are available.",
    "- The Risoe run is saved in this directory because the original paper-table CSV retained only p-values, not the bootstrap-statistic vectors."
  )
  writeLines(lines, file.path(output_dir, "README.md"))
}

main <- function() {
  options <- parse_arguments(commandArgs(trailingOnly = TRUE))
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(options$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  pogo_fast_path <- file.path(
    "real_data", "logistic_gaussian", "screening", "fast",
    "paper_table_B5000_sampleks_fast_rerun_20260718", "pogojump_results.rds"
  )
  pogo_reestimated_path <- file.path(
    "real_data", "logistic_gaussian", "screening", "fast",
    "paper_results_B5000_consolidated", "pogojump_results.rds"
  )
  coxite_path <- file.path(
    "real_data", "logistic_gaussian", "screening", "fast",
    "paper_table_B5000_sampleks_fast_rerun_20260718", "coxite_results.rds"
  )
  animalvegetation_path <- file.path(
    "real_data", "logistic_gaussian", "screening", "fast",
    "paper_table_B5000_sampleks_fast_rerun_20260718", "animalvegetation_results.rds"
  )
  stopifnot(file.exists(pogo_fast_path), file.exists(pogo_reestimated_path),
            file.exists(coxite_path), file.exists(animalvegetation_path))

  pogo_fast <- extract_bootstrap_payload(readRDS(pogo_fast_path), fallback_method = "fast_multiplier")
  pogo_reestimated <- extract_bootstrap_payload(readRDS(pogo_reestimated_path), fallback_method = "reestimated")

  manuscript_p <- c(ks = 0.0980, cvm = 0.5229)
  payloads <- list(
    pogojump_fast = list(
      payload = pogo_fast, domain = "compositional", manuscript_p = manuscript_p,
      source = pogo_fast_path
    ),
    pogojump_reestimated = list(
      payload = pogo_reestimated, domain = "compositional", manuscript_p = c(ks = NA_real_, cvm = NA_real_),
      source = pogo_reestimated_path
    )
  )

  wind_path <- file.path(options$output_dir, "risoe_nov_dec_125m_start4_fast_result.rds")
  if (isTRUE(options$rerun_wind) || !file.exists(wind_path)) {
    wind_run <- run_wind_fast_bootstrap(options$wind_B, options$wind_seed, options$wind_cores)
    saveRDS(wind_run, wind_path)
  } else {
    wind_run <- readRDS(wind_path)
  }
  wind_payload <- extract_bootstrap_payload(wind_run$result, fallback_method = "fast_multiplier")
  payloads$risoe_nov_dec_125m_start4_fast <- list(
    payload = wind_payload,
    domain = "wind",
    manuscript_p = c(ks = 0.4385, cvm = 0.7233),
    source = wind_path
  )

  matched_path <- file.path(options$output_dir, "matched_fast_vs_reestimated_coxite_wind.rds")
  if (isTRUE(options$rerun_matched) || !file.exists(matched_path)) {
    source(file.path("bootstrap", "model_specs.R"))
    source(file.path("bootstrap", "multiplier_bootstrap.R"))
    coxite_result <- readRDS(coxite_path)
    coxite_pair <- run_matched_bootstrap_pair(
      wrapper = multiplier_bootstrap_logistic_gaussian,
      data = coxite_result$data_prep$X_closed,
      control = list(
        logistic_gaussian_quadform_method = "hbe",
        derivative_method = "score_mc",
        derivative_mc_size = 1000L,
        derivative_mc_seed = 20260613L,
        fast_multiplier_cvm_block_size = 50L
      ),
      B = options$matched_B,
      seed = options$matched_seed,
      n_cores = options$matched_cores,
      wrapper_label = "logistic_gaussian_coxite"
    )
    wind_data <- utils::read.csv(wind_run$data_path, stringsAsFactors = FALSE)
    wind_pair <- run_matched_bootstrap_pair(
      wrapper = multiplier_bootstrap_hvmf,
      data = as.matrix(wind_data[, c("x0", "x1", "x2")]),
      control = list(
        hvmf_profile_method = "tabulated",
        derivative_method = "score_mc",
        derivative_mc_size = 1000L,
        derivative_mc_seed = 20260613L,
        fast_multiplier_cvm_block_size = 50L
      ),
      B = options$matched_B,
      seed = options$matched_seed + 1L,
      n_cores = options$matched_cores,
      wrapper_label = "hvmf_risoe_nov_dec_125m_start4"
    )
    matched_pairs <- list(coxite = coxite_pair, risoe_nov_dec_125m_start4 = wind_pair)
    saveRDS(matched_pairs, matched_path)
  } else {
    matched_pairs <- readRDS(matched_path)
  }

  summary_rows <- list()
  pair_rows <- list()
  source_rows <- list()
  index <- 1L
  for (case_id in names(payloads)) {
    entry <- payloads[[case_id]]
    payload <- entry$payload
    save_case_replicates(payload, case_id, options$output_dir)
    plot_case_bootstrap(payload, case_id, options$output_dir)
    source_rows[[index]] <- data.frame(
      case = case_id,
      domain = entry$domain,
      method = payload$method,
      source = entry$source,
      stringsAsFactors = FALSE
    )
    for (statistic in c("ks", "cvm")) {
      summary_rows[[length(summary_rows) + 1L]] <- bootstrap_stat_summary(
        case_id = case_id,
        domain = entry$domain,
        method = payload$method,
        statistic = toupper(statistic),
        observed = payload$observed[[statistic]],
        values = payload[[statistic]],
        stored_p_value = payload$p_values[[statistic]],
        manuscript_p_value = entry$manuscript_p[[statistic]]
      )
    }
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      case = case_id,
      domain = entry$domain,
      method = payload$method,
      correlation_ks_cvm = stats::cor(payload$ks, payload$cvm),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  summary_df <- do.call(rbind, summary_rows)
  pair_df <- do.call(rbind, pair_rows)
  source_df <- do.call(rbind, source_rows)

  fast_summary <- summary_df[summary_df$case == "pogojump_fast", , drop = FALSE]
  reestimated_summary <- summary_df[summary_df$case == "pogojump_reestimated", , drop = FALSE]
  stored_method_comparison <- merge(
    fast_summary[, c("statistic", "p_value", "mc_se", "observed", "bootstrap_mean", "bootstrap_sd")],
    reestimated_summary[, c("statistic", "p_value", "mc_se", "observed", "bootstrap_mean", "bootstrap_sd")],
    by = "statistic", suffixes = c("_fast", "_reestimated")
  )
  stored_method_comparison$p_value_difference <- stored_method_comparison$p_value_fast - stored_method_comparison$p_value_reestimated
  stored_method_comparison$combined_mc_se <- sqrt(stored_method_comparison$mc_se_fast^2 + stored_method_comparison$mc_se_reestimated^2)
  stored_method_comparison$difference_in_mc_se_units <- stored_method_comparison$p_value_difference / stored_method_comparison$combined_mc_se

  matched_summary <- rbind(
    summarize_matched_bootstrap_pair(matched_pairs$coxite, "Coxite", "compositional"),
    summarize_matched_bootstrap_pair(matched_pairs$risoe_nov_dec_125m_start4, "Risoe Nov-Dec 125m", "wind")
  )
  coxite_matched_replicates <- matched_pair_replicates(matched_pairs$coxite)
  wind_matched_replicates <- matched_pair_replicates(matched_pairs$risoe_nov_dec_125m_start4)
  plot_matched_bootstrap_pair(matched_pairs$coxite, "coxite", options$output_dir)
  plot_matched_bootstrap_pair(matched_pairs$risoe_nov_dec_125m_start4, "risoe_nov_dec_125m_start4", options$output_dir)
  maxima_locations <- rbind(
    localize_observed_maxima(pogo_fast$raw, "PogoJump", "compositional"),
    localize_observed_maxima(readRDS(coxite_path)$bootstrap$raw_result, "Coxite", "compositional"),
    localize_observed_maxima(readRDS(animalvegetation_path)$bootstrap$raw_result, "AnimalVegetation", "compositional"),
    localize_observed_maxima(matched_pairs$risoe_nov_dec_125m_start4$fast, "Risoe Nov-Dec 125m", "wind")
  )

  logistic_profile_audit <- audit_logistic_gaussian_hbe_profiles(
    results_dir = dirname(pogo_fast_path),
    output_dir = options$output_dir
  )
  hvmf_profile_audit <- audit_hvmf_tabulated_profile(
    data_path = wind_run$data_path,
    raw_result = wind_run$result,
    output_dir = options$output_dir
  )

  utils::write.csv(summary_df, file.path(options$output_dir, "bootstrap_statistic_summary.csv"), row.names = FALSE)
  utils::write.csv(pair_df, file.path(options$output_dir, "paired_ks_cvm_correlation.csv"), row.names = FALSE)
  utils::write.csv(source_df, file.path(options$output_dir, "source_manifest.csv"), row.names = FALSE)
  utils::write.csv(stored_method_comparison, file.path(options$output_dir, "pogojump_stored_fast_vs_reestimated.csv"), row.names = FALSE)
  utils::write.csv(matched_summary, file.path(options$output_dir, "matched_fast_vs_reestimated_summary.csv"), row.names = FALSE)
  utils::write.csv(coxite_matched_replicates, file.path(options$output_dir, "coxite_matched_fast_vs_reestimated_replicates.csv"), row.names = FALSE)
  utils::write.csv(wind_matched_replicates, file.path(options$output_dir, "risoe_matched_fast_vs_reestimated_replicates.csv"), row.names = FALSE)
  utils::write.csv(maxima_locations, file.path(options$output_dir, "observed_process_top_maxima.csv"), row.names = FALSE)
  utils::write.csv(logistic_profile_audit, file.path(options$output_dir, "logistic_gaussian_hbe_profile_audit.csv"), row.names = FALSE)
  utils::write.csv(hvmf_profile_audit, file.path(options$output_dir, "risoe_hvmf_profile_audit.csv"), row.names = FALSE)
  saveRDS(list(
    summary = summary_df,
    paired = pair_df,
    stored_method_comparison = stored_method_comparison,
    matched_bootstrap_comparison = matched_summary,
    observed_process_top_maxima = maxima_locations,
    logistic_profile_audit = logistic_profile_audit,
    hvmf_profile_audit = hvmf_profile_audit
  ),
          file.path(options$output_dir, "audit_summary.rds"))
  write_report(summary_df, stored_method_comparison, options$output_dir)

  message("Audit artifacts written to: ", normalizePath(options$output_dir))
  print(summary_df[, c("case", "statistic", "B", "observed", "p_value", "mc_se", "bootstrap_sd", "observed_z")], row.names = FALSE)
  print(stored_method_comparison, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
