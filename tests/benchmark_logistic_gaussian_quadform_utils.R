Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_lg_quadform_path <- function(...) {
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

required_lg_quadform_functions <- c(
  "normalize_logistic_gaussian_theta",
  "logistic_gaussian_ilr_matrix",
  "make_logistic_gaussian_simple_calibration_scenario",
  "default_logistic_gaussian_simple_calibration_scenarios",
  "default_logistic_gaussian_composite_calibration_scenarios"
)

if (any(!vapply(required_lg_quadform_functions, exists, logical(1), mode = "function"))) {
  source(resolve_lg_quadform_path("bootstrap", "calibration_study.R"))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

capture_benchmark_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

is_valid_compquad_probability <- function(prob) {
  is.numeric(prob) &&
    length(prob) == 1L &&
    is.finite(prob) &&
    prob >= 0 &&
    prob <= 1
}

format_case_id <- function(prefix, idx) {
  sprintf("%s_%04d", prefix, as.integer(idx))
}

make_lambda_vector <- function(dim, family) {
  dim <- as.integer(dim)
  family <- as.character(family)

  if (identical(family, "balanced")) {
    lambda <- rep.int(1, dim)
  } else if (identical(family, "geometric_mild")) {
    lambda <- exp(seq(0, log(1e-2), length.out = dim))
  } else if (identical(family, "geometric_severe")) {
    lambda <- exp(seq(0, log(1e-6), length.out = dim))
  } else if (identical(family, "one_spike")) {
    lambda <- c(1, rep.int(1e-4, dim - 1L))
  } else {
    stop(sprintf("Unknown lambda family '%s'.", family))
  }

  lambda / mean(lambda)
}

make_delta_vector <- function(dim, family, seed) {
  dim <- as.integer(dim)
  family <- as.character(family)
  set.seed(seed)

  if (identical(family, "central")) {
    delta <- rep.int(0, dim)
  } else if (identical(family, "moderate_random")) {
    delta <- abs(stats::rnorm(dim, mean = 1, sd = 0.35))
  } else if (identical(family, "large_random")) {
    delta <- stats::rlnorm(dim, meanlog = log(4), sdlog = 0.8)
  } else if (identical(family, "tiny_eig_concentrated")) {
    delta <- rep.int(0, dim)
    delta[[dim]] <- 1e4
  } else if (identical(family, "diffuse_extreme")) {
    delta <- rep.int(75, dim)
  } else {
    stop(sprintf("Unknown delta family '%s'.", family))
  }

  as.numeric(delta)
}

make_stress_quadform_cases <- function() {
  dims <- c(2L, 3L, 5L, 10L, 13L)
  lambda_families <- c("balanced", "geometric_mild", "geometric_severe", "one_spike")
  delta_families <- c(
    "central",
    "moderate_random",
    "large_random",
    "tiny_eig_concentrated",
    "diffuse_extreme"
  )
  q_ratios <- c(0.02, 0.10, 0.25, 0.50, 1.00, 2.00, 4.00)

  idx <- 0L
  cases <- vector("list", length(dims) * length(lambda_families) * length(delta_families) * length(q_ratios))

  for (dim in dims) {
    for (lambda_family in lambda_families) {
      lambda <- make_lambda_vector(dim = dim, family = lambda_family)
      for (delta_family in delta_families) {
        seed <- dim * 1000L + match(lambda_family, lambda_families) * 100L + match(delta_family, delta_families)
        delta <- make_delta_vector(dim = dim, family = delta_family, seed = seed)
        mean_stat <- sum(lambda * (1 + delta))
        for (q_ratio in q_ratios) {
          idx <- idx + 1L
          q <- q_ratio * mean_stat
          cases[[idx]] <- data.frame(
            case_id = format_case_id("stress", idx),
            source = "stress",
            regime = paste(dim, lambda_family, delta_family, sep = "::"),
            dim = dim,
            lambda_family = lambda_family,
            delta_family = delta_family,
            q_family = sprintf("ratio_%.2f", q_ratio),
            q_ratio = q_ratio,
            q = q,
            mean_stat = mean_stat,
            lambda_min = min(lambda),
            lambda_max = max(lambda),
            condition_number = max(lambda) / min(lambda),
            delta_max = max(delta),
            delta_sum = sum(delta),
            delta_min_lambda = delta[[which.min(lambda)]],
            lambda_json = paste(signif(lambda, 16), collapse = ";"),
            delta_json = paste(signif(delta, 16), collapse = ";"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  do.call(rbind, cases)
}

make_realistic_logistic_gaussian_cases <- function() {
  scenarios <- c(
    default_logistic_gaussian_simple_calibration_scenarios(),
    default_logistic_gaussian_composite_calibration_scenarios()
  )

  cases <- list()
  idx <- 0L

  for (scenario in scenarios) {
    theta <- normalize_logistic_gaussian_theta(
      list(
        mu_ilr = scenario$sample_params$mu_ilr,
        Sigma_ilr = scenario$sample_params$Sigma_ilr
      ),
      ambient_dim = length(scenario$sample_params$mu_ilr) + 1L
    )

    # Calibration KS statistics intentionally use every observed centre and
    # distance.  A quadratic-form benchmark needs a deterministic finite set
    # of nonzero radii instead, so construct its fixed diagnostic grid from
    # the scenario's null parameter rather than assuming an obsolete explicit
    # `scenario$ks_grid` payload.
    diagnostic_grid <- make_logistic_gaussian_ks_grid(
      mu_ilr = scenario$sample_params$mu_ilr,
      Sigma_ilr = scenario$sample_params$Sigma_ilr,
      ambient_dim = theta$ambient_dim
    )
    omega_grid <- diagnostic_grid$omega_grid
    t_grid <- diagnostic_grid$t_grid
    omega_ilr <- logistic_gaussian_ilr_matrix(omega_grid)

    for (i in seq_len(nrow(omega_ilr))) {
      shift <- theta$mu_ilr - omega_ilr[i, ]
      nu <- as.vector(crossprod(theta$eigenvectors_full, shift))
      delta <- nu^2 / theta$eigenvalues_full
      mean_stat <- sum(theta$eigenvalues_full * (1 + delta))

      for (t_value in t_grid[t_grid > 0]) {
        idx <- idx + 1L
        q <- t_value^2
        cases[[idx]] <- data.frame(
          case_id = format_case_id("realistic", idx),
          source = "realistic",
          regime = scenario$id,
          dim = theta$ilr_dim,
          lambda_family = "from_theta",
          delta_family = sprintf("omega_%02d", i),
          q_family = "ks_grid",
          q_ratio = q / mean_stat,
          q = q,
          mean_stat = mean_stat,
          lambda_min = min(theta$eigenvalues_full),
          lambda_max = max(theta$eigenvalues_full),
          condition_number = max(theta$eigenvalues_full) / min(theta$eigenvalues_full),
          delta_max = max(delta),
          delta_sum = sum(delta),
          delta_min_lambda = delta[[which.min(theta$eigenvalues_full)]],
          lambda_json = paste(signif(theta$eigenvalues_full, 16), collapse = ";"),
          delta_json = paste(signif(delta, 16), collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, cases)
}

make_random_stress_quadform_cases <- function(n_cases = 240L, seed = 123L) {
  set.seed(seed)
  dims <- sample(c(2L, 3L, 5L, 10L, 13L), size = n_cases, replace = TRUE)
  families <- sample(c("balanced", "geometric_mild", "geometric_severe", "one_spike"), size = n_cases, replace = TRUE)

  cases <- vector("list", n_cases)
  for (i in seq_len(n_cases)) {
    dim <- dims[[i]]
    lambda <- make_lambda_vector(dim = dim, family = families[[i]])
    random_direction <- stats::rnorm(dim)
    random_direction <- random_direction / sqrt(sum(random_direction^2))
    magnitude <- stats::rlnorm(1L, meanlog = log(5), sdlog = 1)
    delta <- magnitude * random_direction^2 / mean(random_direction^2)
    if (runif(1) < 0.25) {
      delta[[which.min(lambda)]] <- delta[[which.min(lambda)]] + 1e4
    }
    mean_stat <- sum(lambda * (1 + delta))
    q_ratio <- exp(stats::runif(1L, min = log(0.01), max = log(8)))
    q <- q_ratio * mean_stat

    cases[[i]] <- data.frame(
      case_id = format_case_id("random", i),
      source = "random",
      regime = sprintf("random_dim_%02d", dim),
      dim = dim,
      lambda_family = families[[i]],
      delta_family = "random",
      q_family = "random_ratio",
      q_ratio = q_ratio,
      q = q,
      mean_stat = mean_stat,
      lambda_min = min(lambda),
      lambda_max = max(lambda),
      condition_number = max(lambda) / min(lambda),
      delta_max = max(delta),
      delta_sum = sum(delta),
      delta_min_lambda = delta[[which.min(lambda)]],
      lambda_json = paste(signif(lambda, 16), collapse = ";"),
      delta_json = paste(signif(delta, 16), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, cases)
}

default_quadform_benchmark_cases <- function(random_n = 240L,
                                             random_seed = 123L,
                                             stress = TRUE,
                                             realistic = TRUE,
                                             random_cases = TRUE) {
  parts <- list()

  if (isTRUE(stress)) {
    parts[[length(parts) + 1L]] <- make_stress_quadform_cases()
  }
  if (isTRUE(realistic)) {
    parts[[length(parts) + 1L]] <- make_realistic_logistic_gaussian_cases()
  }
  if (isTRUE(random_cases) && random_n > 0L) {
    parts[[length(parts) + 1L]] <- make_random_stress_quadform_cases(n_cases = random_n, seed = random_seed)
  }

  do.call(rbind, parts)
}

parse_case_vector <- function(text) {
  as.numeric(strsplit(text, ";", fixed = TRUE)[[1L]])
}

evaluate_quadform_method_once <- function(lambda,
                                          delta,
                                          q,
                                          method) {
  start <- proc.time()[[3L]]
  warnings <- character(0)
  status <- "ok"
  error_message <- NA_character_
  value <- NA_real_
  ifault <- NA_integer_
  qq <- NA_real_
  abserr <- NA_real_

  captured <- tryCatch(
    capture_benchmark_warnings({
      if (identical(method, "farebrother")) {
        CompQuadForm::farebrother(
          q = q,
          lambda = lambda,
          h = rep.int(1, length(lambda)),
          delta = delta,
          maxit = 100000L,
          eps = 1e-8
        )
      } else if (identical(method, "davies")) {
        CompQuadForm::davies(
          q = q,
          lambda = lambda,
          h = rep.int(1, length(lambda)),
          delta = delta,
          lim = 20000L,
          acc = 1e-8
        )
      } else if (identical(method, "imhof")) {
        CompQuadForm::imhof(
          q = q,
          lambda = lambda,
          h = rep.int(1, length(lambda)),
          delta = delta,
          epsabs = 1e-8,
          epsrel = 1e-8,
          limit = 20000L
        )
      } else if (identical(method, "sphunif_sw")) {
        sphunif::p_wschisq(
          x = q,
          weights = lambda,
          dfs = rep.int(1, length(lambda)),
          ncps = delta,
          method = "SW"
        )
      } else if (identical(method, "sphunif_hbe")) {
        sphunif::p_wschisq(
          x = q,
          weights = lambda,
          dfs = rep.int(1, length(lambda)),
          ncps = delta,
          method = "HBE"
        )
      } else {
        stop(sprintf("Unknown method '%s'.", method))
      }
    }),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  elapsed <- proc.time()[[3L]] - start

  if (!is.null(captured)) {
    warnings <- captured$warnings
    result <- captured$value

    if (identical(method, "farebrother") || identical(method, "davies")) {
      ifault <- as.integer(result$ifault %||% NA_integer_)
      qq <- as.numeric(result$Qq %||% NA_real_)
      value <- 1 - qq
      if (!is.na(ifault) && ifault != 0L) {
        status <- "ifault"
      }
    } else if (identical(method, "imhof")) {
      qq <- as.numeric(result$Qq %||% NA_real_)
      abserr <- as.numeric(result$abserr %||% NA_real_)
      value <- 1 - qq
    } else {
      value <- as.numeric(result)[[1L]]
    }
  }

  if (!is_valid_compquad_probability(value)) {
    if (identical(status, "ok")) {
      status <- "invalid_value"
    }
  }

  data.frame(
    method = method,
    elapsed = elapsed,
    status = status,
    value = value,
    ifault = ifault,
    qq = qq,
    abserr = abserr,
    n_warnings = length(warnings),
    warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

benchmark_quadform_case <- function(case_row,
                                    methods = c("farebrother", "davies", "imhof", "sphunif_sw", "sphunif_hbe")) {
  lambda <- parse_case_vector(case_row$lambda_json)
  delta <- parse_case_vector(case_row$delta_json)
  q <- as.numeric(case_row$q)

  per_method <- lapply(methods, function(method) {
    cbind(case_row, evaluate_quadform_method_once(lambda = lambda, delta = delta, q = q, method = method))
  })

  do.call(rbind, per_method)
}

derive_reference_values <- function(results) {
  split_by_case <- split(results, results$case_id)

  reference_rows <- lapply(split_by_case, function(df_case) {
    imhof_row <- df_case[df_case$method == "imhof", , drop = FALSE]
    fare_row <- df_case[df_case$method == "farebrother", , drop = FALSE]

    ref_method <- NA_character_
    ref_value <- NA_real_
    agreement_abs <- NA_real_

    imhof_ok <- nrow(imhof_row) == 1L && identical(imhof_row$status, "ok") && is_valid_compquad_probability(imhof_row$value)
    fare_ok <- nrow(fare_row) == 1L && identical(fare_row$status, "ok") && is_valid_compquad_probability(fare_row$value)

    if (imhof_ok && fare_ok) {
      agreement_abs <- abs(imhof_row$value - fare_row$value)
      ref_method <- "imhof"
      ref_value <- imhof_row$value
    } else if (imhof_ok) {
      ref_method <- "imhof"
      ref_value <- imhof_row$value
    } else if (fare_ok) {
      ref_method <- "farebrother"
      ref_value <- fare_row$value
    }

    data.frame(
      case_id = df_case$case_id[[1L]],
      reference_method = ref_method,
      reference_value = ref_value,
      exact_agreement_abs = agreement_abs,
      stringsAsFactors = FALSE
    )
  })

  ref_df <- do.call(rbind, reference_rows)
  merged <- merge(results, ref_df, by = "case_id", all.x = TRUE, sort = FALSE)
  merged$abs_error_vs_reference <- abs(merged$value - merged$reference_value)
  merged$rel_error_vs_reference <- merged$abs_error_vs_reference / pmax(abs(merged$reference_value), 1e-12)
  merged
}

summarize_quadform_benchmark <- function(results) {
  split_methods <- split(results, results$method)

  method_summary <- do.call(rbind, lapply(split_methods, function(df_method) {
    valid_mask <- is.finite(df_method$value)
    ok_mask <- identical(df_method$status, "ok")
    valid_exact_ref <- is.finite(df_method$abs_error_vs_reference)

    data.frame(
      method = df_method$method[[1L]],
      n = nrow(df_method),
      n_status_ok = sum(df_method$status == "ok", na.rm = TRUE),
      n_status_ifault = sum(df_method$status == "ifault", na.rm = TRUE),
      n_status_error = sum(df_method$status == "error", na.rm = TRUE),
      n_status_invalid = sum(df_method$status == "invalid_value", na.rm = TRUE),
      median_elapsed = stats::median(df_method$elapsed, na.rm = TRUE),
      p95_elapsed = as.numeric(stats::quantile(df_method$elapsed, probs = 0.95, na.rm = TRUE, names = FALSE)),
      max_elapsed = max(df_method$elapsed, na.rm = TRUE),
      median_abs_error = stats::median(df_method$abs_error_vs_reference[valid_exact_ref], na.rm = TRUE),
      p95_abs_error = as.numeric(stats::quantile(df_method$abs_error_vs_reference[valid_exact_ref], probs = 0.95, na.rm = TRUE, names = FALSE)),
      max_abs_error = max(df_method$abs_error_vs_reference[valid_exact_ref], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  farebrother_rows <- results[results$method == "farebrother", , drop = FALSE]
  slow_threshold <- as.numeric(stats::quantile(farebrother_rows$elapsed, probs = 0.95, na.rm = TRUE, names = FALSE))
  farebrother_rows$slow_case <- farebrother_rows$elapsed >= slow_threshold
  farebrother_rows$cond_bin <- cut(
    farebrother_rows$condition_number,
    breaks = c(-Inf, 1e2, 1e4, Inf),
    labels = c("<=1e2", "(1e2,1e4]", ">1e4")
  )
  farebrother_rows$delta_bin <- cut(
    farebrother_rows$delta_max,
    breaks = c(-Inf, 1, 10, 1e3, Inf),
    labels = c("<=1", "(1,10]", "(10,1e3]", ">1e3")
  )
  farebrother_rows$q_ratio_bin <- cut(
    farebrother_rows$q_ratio,
    breaks = c(-Inf, 0.1, 0.5, 1.5, Inf),
    labels = c("very_low", "low", "central", "upper")
  )

  by_regime <- stats::aggregate(
    slow_case ~ cond_bin + delta_bin + q_ratio_bin,
    data = farebrother_rows,
    FUN = mean
  )
  by_regime_n <- stats::aggregate(
    slow_case ~ cond_bin + delta_bin + q_ratio_bin,
    data = farebrother_rows,
    FUN = length
  )

  regime_summary <- data.frame(
    cond_bin = by_regime$cond_bin,
    delta_bin = by_regime$delta_bin,
    q_ratio_bin = by_regime$q_ratio_bin,
    slow_rate = by_regime$slow_case,
    n = by_regime_n$slow_case
  )

  list(
    method_summary = method_summary,
    farebrother_regime_summary = regime_summary,
    farebrother_slow_threshold = slow_threshold
  )
}

summarize_dispatcher_rules <- function(results) {
  clip_probability <- function(x) pmin(pmax(x, 0), 1)

  wide <- reshape(
    results[, c(
      "case_id", "method", "elapsed", "status", "value", "reference_value",
      "condition_number", "delta_max", "dim", "q_ratio"
    )],
    idvar = "case_id",
    timevar = "method",
    direction = "wide"
  )

  wide$cond <- wide$condition_number.farebrother
  wide$delta <- wide$delta_max.farebrother
  wide$dim_case <- wide$dim.farebrother
  wide$qratio <- wide$q_ratio.farebrother
  wide$fare_ok <- wide$status.farebrother == "ok"
  wide$imhof_clipped <- clip_probability(wide$value.imhof)
  wide$reference <- wide$reference_value.farebrother

  score_rule <- function(use_imhof, rule_name) {
    chosen <- ifelse(
      use_imhof,
      wide$imhof_clipped,
      ifelse(wide$fare_ok, wide$value.farebrother, wide$imhof_clipped)
    )
    elapsed <- ifelse(
      use_imhof,
      wide$elapsed.imhof,
      ifelse(wide$fare_ok, wide$elapsed.farebrother, wide$elapsed.farebrother + wide$elapsed.imhof)
    )
    abs_error <- abs(chosen - wide$reference)

    data.frame(
      rule = rule_name,
      share_imhof = mean(use_imhof),
      total_time = sum(elapsed, na.rm = TRUE),
      median_time = stats::median(elapsed, na.rm = TRUE),
      p95_time = as.numeric(stats::quantile(elapsed, probs = 0.95, na.rm = TRUE, names = FALSE)),
      max_abs_error = max(abs_error, na.rm = TRUE),
      p95_abs_error = as.numeric(stats::quantile(abs_error, probs = 0.95, na.rm = TRUE, names = FALSE)),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  rules <- list(
    fallback_after_farebrother = rep(FALSE, nrow(wide)),
    cond_gt_1e4 = wide$cond > 1e4,
    cond_gt_1e4_or_dim10_delta10 = wide$cond > 1e4 | (wide$dim_case >= 10 & wide$delta > 10),
    cond_gt_1e4_or_dim5_delta10_qhigh = wide$cond > 1e4 | (wide$dim_case >= 5 & wide$delta > 10 & wide$qratio > 0.1)
  )

  out <- do.call(
    rbind,
    lapply(names(rules), function(name) score_rule(rules[[name]], rule_name = name))
  )

  out[order(out$total_time, out$p95_time), , drop = FALSE]
}

write_quadform_benchmark_report <- function(output_dir, summary_list, results) {
  method_summary <- summary_list$method_summary
  slow_threshold <- summary_list$farebrother_slow_threshold
  dispatcher_summary <- summary_list$dispatcher_summary
  best_rule <- dispatcher_summary[1, , drop = FALSE]

  imhof_row <- method_summary[method_summary$method == "imhof", , drop = FALSE]
  fare_row <- method_summary[method_summary$method == "farebrother", , drop = FALSE]
  davies_row <- method_summary[method_summary$method == "davies", , drop = FALSE]
  sw_row <- method_summary[method_summary$method == "sphunif_sw", , drop = FALSE]
  hbe_row <- method_summary[method_summary$method == "sphunif_hbe", , drop = FALSE]

  lines <- c(
    "# Logistic Gaussian Quadratic-Form Benchmark",
    "",
    "This benchmark compares candidate methods for evaluating weighted sums of noncentral chi-squared variables arising in the logistic Gaussian distance-profile calculations.",
    "",
    "Reference rule:",
    "- use `imhof` as the reference when both `imhof` and `farebrother` are valid;",
    "- otherwise fall back to the available exact method among `imhof` and `farebrother`.",
    "",
    sprintf("Farebrother slow-case threshold: %.6f seconds (empirical 95th percentile).", slow_threshold),
    "",
    "Headline findings:",
    sprintf(
      "- `farebrother`: median %.6f s, 95th percentile %.6f s, max abs. error %.3e, %d non-ok cases.",
      fare_row$median_elapsed,
      fare_row$p95_elapsed,
      fare_row$max_abs_error,
      fare_row$n - fare_row$n_status_ok
    ),
    sprintf(
      "- `imhof`: median %.6f s, 95th percentile %.6f s, used as exact reference.",
      imhof_row$median_elapsed,
      imhof_row$p95_elapsed
    ),
    sprintf(
      "- `davies`: %d non-ok cases out of %d, so it is not suitable as the default exact route without explicit failure handling.",
      davies_row$n - davies_row$n_status_ok,
      davies_row$n
    ),
    sprintf(
      "- `sphunif_sw`: median abs. error %.3e; `sphunif_hbe`: median abs. error %.3e.",
      sw_row$median_abs_error,
      hbe_row$median_abs_error
    ),
    sprintf(
      "- Best dispatcher among the benchmarked simple rules: `%s`, with total benchmark time %.2f s, 95th percentile %.5f s, and maximum absolute deviation %.3e versus the benchmark reference.",
      best_rule$rule,
      best_rule$total_time,
      best_rule$p95_time,
      best_rule$max_abs_error
    ),
    "",
    "Interpretation:",
    "- `farebrother` is the fastest exact method whenever it succeeds and remains numerically stable.",
    "- `imhof` is slower but provides a robust exact fallback and a reliable benchmark reference.",
    "- `davies` can fail with non-zero `ifault` even in moderate cases, so it should not be the default route.",
    "- `SW` and `HBE` are useful only as fast approximations, not as exact replacements.",
    sprintf(
      "- The benchmark favours the dispatcher `%s`: send clearly ill-conditioned or strongly noncentral medium/high-dimensional cases directly to `imhof`, and otherwise try `farebrother` first.",
      best_rule$rule
    )
  )

  writeLines(lines, con = file.path(output_dir, "recommendation.md"))
}

plot_quadform_benchmark <- function(results, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  ggplot2 <- asNamespace("ggplot2")

  p_time <- ggplot2$ggplot(results, ggplot2$aes(x = elapsed, colour = method)) +
    ggplot2$stat_ecdf(linewidth = 0.7) +
    ggplot2$scale_x_log10() +
    ggplot2$labs(
      title = "Empirical runtime distribution by method",
      x = "Elapsed time per scalar CDF evaluation (seconds, log scale)",
      y = "ECDF"
    ) +
    ggplot2$theme_minimal()

  fare_df <- results[results$method == "farebrother", , drop = FALSE]
  p_fare <- ggplot2$ggplot(
    fare_df,
    ggplot2$aes(x = condition_number, y = elapsed, colour = delta_max, shape = source)
  ) +
    ggplot2$geom_point(alpha = 0.8, size = 2) +
    ggplot2$scale_x_log10() +
    ggplot2$scale_y_log10() +
    ggplot2$labs(
      title = "Farebrother runtime across parameter regimes",
      x = "Condition number of lambda (log scale)",
      y = "Elapsed time (seconds, log scale)",
      colour = "max delta"
    ) +
    ggplot2$theme_minimal()

  approx_df <- results[results$method %in% c("sphunif_sw", "sphunif_hbe"), , drop = FALSE]
  p_error <- ggplot2$ggplot(
    approx_df,
    ggplot2$aes(x = elapsed, y = abs_error_vs_reference, colour = method)
  ) +
    ggplot2$geom_point(alpha = 0.6, size = 1.5) +
    ggplot2$scale_x_log10() +
    ggplot2$scale_y_log10() +
    ggplot2$labs(
      title = "Approximation error versus runtime",
      x = "Elapsed time (seconds, log scale)",
      y = "Absolute error versus exact reference (log scale)"
    ) +
    ggplot2$theme_minimal()

  ggplot2$ggsave(file.path(output_dir, "runtime_ecdf.png"), plot = p_time, width = 8, height = 5, dpi = 160)
  ggplot2$ggsave(file.path(output_dir, "farebrother_runtime_regimes.png"), plot = p_fare, width = 8, height = 5, dpi = 160)
  ggplot2$ggsave(file.path(output_dir, "approximation_error_vs_runtime.png"), plot = p_error, width = 8, height = 5, dpi = 160)
}

run_logistic_gaussian_quadform_benchmark <- function(output_dir = file.path("tests", "benchmark_outputs", "logistic_gaussian_quadform"),
                                                     n_cores = 12L,
                                                     random_n = 240L,
                                                     random_seed = 123L,
                                                     max_cases = NULL,
                                                     stress = TRUE,
                                                     realistic = TRUE,
                                                     random_cases = TRUE,
                                                     methods = c("farebrother", "davies", "imhof", "sphunif_sw", "sphunif_hbe"),
                                                     save_plots = TRUE) {
  n_cores <- max(1L, min(12L, as.integer(n_cores)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cases <- default_quadform_benchmark_cases(
    random_n = random_n,
    random_seed = random_seed,
    stress = stress,
    realistic = realistic,
    random_cases = random_cases
  )
  if (!is.null(max_cases)) {
    max_cases <- as.integer(max_cases)
    cases <- utils::head(cases, max_cases)
  }
  split_cases <- split(cases, seq_len(nrow(cases)))

  worker_fun <- function(case_df) {
    benchmark_quadform_case(case_row = case_df, methods = methods)
  }

  results_list <- parallel::mclapply(
    split_cases,
    worker_fun,
    mc.cores = n_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )

  results <- do.call(rbind, results_list)
  results <- derive_reference_values(results)
  summary_list <- summarize_quadform_benchmark(results)
  summary_list$dispatcher_summary <- summarize_dispatcher_rules(results)

  utils::write.csv(cases, file.path(output_dir, "quadform_cases.csv"), row.names = FALSE)
  utils::write.csv(results, file.path(output_dir, "quadform_results.csv"), row.names = FALSE)
  utils::write.csv(summary_list$method_summary, file.path(output_dir, "quadform_method_summary.csv"), row.names = FALSE)
  utils::write.csv(summary_list$farebrother_regime_summary, file.path(output_dir, "farebrother_regime_summary.csv"), row.names = FALSE)
  utils::write.csv(summary_list$dispatcher_summary, file.path(output_dir, "dispatcher_rule_summary.csv"), row.names = FALSE)

  slow_cases <- results[
    results$method == "farebrother" &
      results$elapsed >= summary_list$farebrother_slow_threshold,
    ,
    drop = FALSE
  ]
  slow_cases <- slow_cases[order(-slow_cases$elapsed), ]
  utils::write.csv(slow_cases, file.path(output_dir, "farebrother_slow_cases.csv"), row.names = FALSE)

  write_quadform_benchmark_report(output_dir = output_dir, summary_list = summary_list, results = results)

  if (save_plots) {
    plot_quadform_benchmark(results = results, output_dir = output_dir)
  }

  sink(file.path(output_dir, "sessionInfo.txt"))
  print(utils::sessionInfo())
  sink()

  saveRDS(
    list(
      cases = cases,
      results = results,
      summaries = summary_list,
      config = list(
        output_dir = output_dir,
        n_cores = n_cores,
        random_n = random_n,
        random_seed = random_seed,
        max_cases = max_cases,
        stress = stress,
        realistic = realistic,
        random_cases = random_cases,
        methods = methods
      )
    ),
    file.path(output_dir, "quadform_benchmark.rds")
  )

  invisible(
    list(
      cases = cases,
      results = results,
      summaries = summary_list,
      output_dir = output_dir
    )
  )
}

mc_reference_for_quadform_case <- function(case_row,
                                           methods = c("imhof", "sphunif_hbe", "sphunif_sw", "farebrother"),
                                           M = 2000000L,
                                           chunk_size = 200000L,
                                           seed = 1L) {
  lambda <- parse_case_vector(case_row$lambda_json)
  delta <- parse_case_vector(case_row$delta_json)
  q <- as.numeric(case_row$q)
  df_vec <- rep.int(1L, length(lambda))

  set.seed(seed)
  remaining <- as.integer(M)
  hits <- 0L
  while (remaining > 0L) {
    current_chunk <- min(remaining, as.integer(chunk_size))
    sims <- vapply(
      seq_along(lambda),
      function(j) stats::rchisq(current_chunk, df = df_vec[[j]], ncp = delta[[j]]),
      numeric(current_chunk)
    )
    sims <- as.matrix(sims)
    stat_values <- as.vector(sims %*% lambda)
    hits <- hits + sum(stat_values <= q)
    remaining <- remaining - current_chunk
  }

  p_mc <- hits / M
  se_mc <- sqrt(p_mc * (1 - p_mc) / M)

  benchmark_rows <- benchmark_quadform_case(case_row = case_row, methods = methods)
  benchmark_rows$mc_reference <- p_mc
  benchmark_rows$mc_se <- se_mc
  benchmark_rows$abs_error_vs_mc <- abs(benchmark_rows$value - p_mc)
  benchmark_rows
}

run_quadform_mc_reference_study <- function(output_dir,
                                            benchmark_results_csv,
                                            n_cores = 12L,
                                            M = 2000000L,
                                            chunk_size = 200000L,
                                            seed = 123L) {
  n_cores <- max(1L, min(12L, as.integer(n_cores)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- utils::read.csv(benchmark_results_csv, stringsAsFactors = FALSE)
  case_catalog <- unique(results[, c(
    "case_id", "source", "regime", "dim", "lambda_family", "delta_family",
    "q_family", "q_ratio", "q", "mean_stat", "lambda_min", "lambda_max",
    "condition_number", "delta_max", "delta_sum", "delta_min_lambda",
    "lambda_json", "delta_json", "reference_method"
  )])

  unresolved_cases <- case_catalog[is.na(case_catalog$reference_method), , drop = FALSE]
  if (nrow(unresolved_cases) == 0L) {
    stop("No unresolved cases found in benchmark results.")
  }

  case_list <- split(unresolved_cases, seq_len(nrow(unresolved_cases)))
  rows <- parallel::mclapply(
    seq_along(case_list),
    function(i) {
      mc_reference_for_quadform_case(
        case_row = case_list[[i]],
        M = M,
        chunk_size = chunk_size,
        seed = seed + i
      )
    },
    mc.cores = n_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )

  mc_results <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(mc_results, mc_results$method), function(df) {
    data.frame(
      method = df$method[[1L]],
      n = nrow(df),
      median_elapsed = stats::median(df$elapsed, na.rm = TRUE),
      p95_elapsed = as.numeric(stats::quantile(df$elapsed, probs = 0.95, na.rm = TRUE, names = FALSE)),
      median_abs_error_vs_mc = stats::median(df$abs_error_vs_mc, na.rm = TRUE),
      p95_abs_error_vs_mc = as.numeric(stats::quantile(df$abs_error_vs_mc, probs = 0.95, na.rm = TRUE, names = FALSE)),
      max_abs_error_vs_mc = max(df$abs_error_vs_mc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  utils::write.csv(mc_results, file.path(output_dir, "quadform_mc_reference_results.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "quadform_mc_reference_summary.csv"), row.names = FALSE)
  sink(file.path(output_dir, "sessionInfo.txt"))
  print(utils::sessionInfo())
  sink()

  saveRDS(
    list(
      mc_results = mc_results,
      summary = summary,
      config = list(
        benchmark_results_csv = benchmark_results_csv,
        n_cores = n_cores,
        M = M,
        chunk_size = chunk_size,
        seed = seed
      )
    ),
    file.path(output_dir, "quadform_mc_reference.rds")
  )

  invisible(list(mc_results = mc_results, summary = summary, output_dir = output_dir))
}
