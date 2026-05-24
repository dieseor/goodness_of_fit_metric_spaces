resolve_hvmf_gig_study_path <- function(...) {
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

utils_path_hvmf_gig <- normalizePath(resolve_hvmf_gig_study_path("utils.R"), winslash = "/", mustWork = TRUE)
model_specs_path_hvmf_gig <- normalizePath(resolve_hvmf_gig_study_path("bootstrap", "model_specs.R"), winslash = "/", mustWork = TRUE)
multiplier_bootstrap_path_hvmf_gig <- normalizePath(resolve_hvmf_gig_study_path("bootstrap", "multiplier_bootstrap.R"), winslash = "/", mustWork = TRUE)

if (!exists("minkowski_inner_product", mode = "function")) {
  source(utils_path_hvmf_gig)
}
if (!exists("make_hvmf_spec", mode = "function")) {
  source(model_specs_path_hvmf_gig)
}
if (!exists("multiplier_bootstrap_hvmf", mode = "function")) {
  source(multiplier_bootstrap_path_hvmf_gig)
}

hvmf_gig_fixed_mu <- function() {
  t0 <- stats::qnorm(0.25, mean = 0, sd = 1 / 4)
  c(
    cosh(abs(t0)),
    sinh(abs(t0)) * sign(t0) / sqrt(2),
    sinh(abs(t0)) * sign(t0) / sqrt(2)
  )
}

format_hvmf_gig_kappa_tag <- function(kappa) {
  kappa <- as.numeric(kappa)
  if (!is.finite(kappa) || length(kappa) != 1L) {
    stop("`kappa` must be a finite scalar.")
  }

  tag <- format(kappa, trim = TRUE, scientific = FALSE)
  tag <- gsub("\\.", "p", tag, fixed = FALSE)
  tag <- gsub("-", "m", tag, fixed = TRUE)
  tag
}

format_hvmf_gig_warning_block <- function(warnings) {
  warnings <- unique(as.character(warnings))
  if (length(warnings) == 0L) {
    return("warnings: none")
  }

  c("warnings:", paste0("- ", warnings))
}

capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(warning_condition) {
      warnings <<- c(warnings, conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = warnings)
}

setup_hvmf_gig_cluster <- function(n_cores) {
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }

  if (n_cores == 1L) {
    return(NULL)
  }

  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(
    cl,
    c("utils_path_hvmf_gig", "model_specs_path_hvmf_gig", "multiplier_bootstrap_path_hvmf_gig"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    source(utils_path_hvmf_gig)
    source(model_specs_path_hvmf_gig)
    source(multiplier_bootstrap_path_hvmf_gig)
    NULL
  })

  cl
}

run_hvmf_gig_mle_monte_carlo_study <- function(output_csv = file.path("bootstrap", "results", "hvmf_gig_mle_monte_carlo_summary.csv"),
                                               kappas = c(50, 200),
                                               n_values = c(50, 100, 200),
                                               n_replicates = 1000,
                                               seed = 123,
                                               mu = hvmf_gig_fixed_mu()) {
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  mu <- as.numeric(mu)

  summary_rows <- list()
  base_seed <- as.integer(seed)

  for (kappa in kappas) {
    for (n in n_values) {
      scenario_seed <- base_seed + 10000L * as.integer(kappa) + as.integer(n)
      set.seed(scenario_seed)
      start_time <- Sys.time()

      estimates <- replicate(as.integer(n_replicates), {
        x <- rhvmf_h2_gig(n = n, mu = mu, kappa = kappa, check = TRUE)
        fit <- hvmf_mle_h2(x)
        c(
          kappa_hat = fit$kappa,
          angular_error = hyperbolic_geodesic_distance_h2(fit$mu, mu)
        )
      })

      estimates <- as.data.frame(t(estimates))
      elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      theoretical_mean_kappa_hat <- as.numeric(kappa) * as.integer(n) / (as.integer(n) - 2L)
      mean_kappa_hat <- mean(estimates$kappa_hat)

      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        n = as.integer(n),
        kappa = as.numeric(kappa),
        n_replicates = as.integer(n_replicates),
        mean_kappa_hat = mean_kappa_hat,
        sd_kappa_hat = stats::sd(estimates$kappa_hat),
        theoretical_mean_kappa_hat = theoretical_mean_kappa_hat,
        relative_bias_vs_kappa = (mean_kappa_hat - kappa) / kappa,
        relative_error_vs_theoretical_mean = abs(mean_kappa_hat - theoretical_mean_kappa_hat) / theoretical_mean_kappa_hat,
        mean_angular_error = mean(estimates$angular_error),
        elapsed_seconds = elapsed_seconds,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  utils::write.csv(summary_df, output_csv, row.names = FALSE)
  summary_df
}

run_hvmf_gig_cvm_sample <- function(data,
                                   null_theta,
                                   B = 1000,
                                   n_cores = 3,
                                   seed = NULL,
                                   control = list(),
                                   cluster = NULL) {
  spec <- make_hvmf_spec(unknown_param = "both")
  null <- list(type = "simple", theta = null_theta)
  data_normalized <- normalize_hvmf_data(data, control)
  theta_hat <- spec$fit_theta(
    data = data_normalized,
    weights = NULL,
    null = null,
    control = control
  )

  cvm_prep <- prepare_cvm_observed_data(
    data = data_normalized,
    spec = spec,
    theta_hat = theta_hat,
    control = control
  )

  multiplier_spec <- resolve_multiplier_spec(NULL)
  scale_factor <- multiplier_spec$mean / multiplier_spec$sd
  raw_multiplier_matrix <- generate_multiplier_matrix(
    B = as.integer(B),
    n = nrow(data_normalized),
    multiplier_spec = multiplier_spec,
    seed = seed
  )
  normalized_multiplier_matrix <- raw_multiplier_matrix / rowMeans(raw_multiplier_matrix)

  n_cores_effective <- min(as.integer(n_cores), as.integer(B))
  chunk_ids <- split(seq_len(as.integer(B)), rep(seq_len(n_cores_effective), length.out = as.integer(B)))
  weight_chunks <- lapply(chunk_ids, function(indices) {
    normalized_multiplier_matrix[indices, , drop = FALSE]
  })

  if (n_cores_effective == 1L || is.null(cluster)) {
    chunk_results <- lapply(weight_chunks, run_bootstrap_chunk,
      spec = spec,
      data = data_normalized,
      null = null,
      control = control,
      scale_factor = scale_factor,
      ks_prep = NULL,
      cvm_prep = cvm_prep,
      want_ks = FALSE,
      want_cvm = TRUE
    )
  } else {
    chunk_results <- parallel::parLapply(cluster, weight_chunks, function(chunk) {
      run_bootstrap_chunk(
        weight_chunk = chunk,
        spec = spec,
        data = data,
        null = null,
        control = control,
        scale_factor = scale_factor,
        ks_prep = NULL,
        cvm_prep = cvm_prep,
        want_ks = FALSE,
        want_cvm = TRUE
      )
    })
  }

  bootstrap_statistics <- unlist(lapply(chunk_results, `[[`, "cvm"), use.names = FALSE)
  inference <- compute_inference_summary(
    observed_statistics = list(cvm = cvm_prep$statistic),
    bootstrap_statistics = list(cvm = bootstrap_statistics),
    alpha = 0.05
  )

  list(
    theta_hat = theta_hat,
    cvm_prep = cvm_prep,
    bootstrap_statistics = bootstrap_statistics,
    inference = inference$cvm,
    n = nrow(data_normalized)
  )
}

run_hvmf_gig_simple_cvm_calibration_study <- function(output_dir = file.path("bootstrap", "results"),
                                                      scenarios = list(
                                                        list(n = 100L, kappa = 1.5, M = 1000L, B = 1000L, n_cores = 3L)
                                                      ),
                                                      seed = 123,
                                                      mu = hvmf_gig_fixed_mu()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mu <- as.numeric(mu)

  cluster <- setup_hvmf_gig_cluster(max(vapply(scenarios, function(scenario) scenario$n_cores %||% 1L, integer(1))))
  if (!is.null(cluster)) {
    on.exit(parallel::stopCluster(cluster), add = TRUE)
  }

  summary_rows <- list()
  raw_rows <- list()
  log_lines <- c(
    sprintf("start_time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("sampler: rhvmf_h2_gig"),
    sprintf("mu: %s", paste(signif(mu, 6), collapse = ", "))
  )

  for (scenario_index in seq_along(scenarios)) {
    scenario <- scenarios[[scenario_index]]
    n <- as.integer(scenario$n)
    kappa <- as.numeric(scenario$kappa)
    M <- as.integer(scenario$M %||% 1000L)
    B <- as.integer(scenario$B %||% 1000L)
    n_cores <- as.integer(scenario$n_cores %||% 3L)
    scenario_seed <- as.integer(seed) + 1000L * scenario_index
    null_theta <- list(mu = mu, kappa = kappa)
    scenario_raw_rows <- list()

    log_lines <- c(
      log_lines,
      sprintf("scenario_%d: n=%d kappa=%s M=%d B=%d n_cores=%d", scenario_index, n, format(kappa, trim = TRUE, scientific = FALSE), M, B, n_cores),
      sprintf("scenario_%d_seed: %d", scenario_index, scenario_seed)
    )

    warnings_all <- character(0)
    p_values <- numeric(M)
    observed_statistics <- numeric(M)
    replicate_elapsed <- numeric(M)

    scenario_start <- Sys.time()
    for (m in seq_len(M)) {
      replicate_seed <- scenario_seed + m
      sample_result <- capture_warnings(
        rhvmf_h2_gig(n = n, mu = mu, kappa = kappa, check = TRUE)
      )
      warnings_all <- c(warnings_all, sample_result$warnings)

      bootstrap_result <- capture_warnings(
        run_hvmf_gig_cvm_sample(
          data = sample_result$value,
          null_theta = null_theta,
          B = B,
          n_cores = n_cores,
          seed = replicate_seed,
          cluster = cluster
        )
      )
      warnings_all <- c(warnings_all, bootstrap_result$warnings)

      p_values[m] <- bootstrap_result$value$inference$p_value
      observed_statistics[m] <- bootstrap_result$value$cvm_prep$statistic
      replicate_elapsed[m] <- NA_real_

      raw_rows[[length(raw_rows) + 1L]] <- data.frame(
        scenario_index = scenario_index,
        n = n,
        kappa = kappa,
        replicate = m,
        p_value = p_values[m],
        observed_cvm = observed_statistics[m],
        reject_5pct = as.integer(p_values[m] <= 0.05),
        reject_10pct = as.integer(p_values[m] <= 0.10),
        stringsAsFactors = FALSE
      )
      scenario_raw_rows[[length(scenario_raw_rows) + 1L]] <- raw_rows[[length(raw_rows)]]
    }

    elapsed_seconds <- as.numeric(difftime(Sys.time(), scenario_start, units = "secs"))
    p_value_quantiles <- stats::quantile(p_values, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE, type = 8)
    output_rds <- file.path(
      output_dir,
      sprintf(
        "hvmf_simple_gig_cvm_calibration_n%d_kappa%s.rds",
        n,
        format_hvmf_gig_kappa_tag(kappa)
      )
    )

    scenario_result <- list(
      config = list(
        n = n,
        kappa = kappa,
        M = M,
        B = B,
        n_cores = n_cores,
        seed = scenario_seed,
        mu = mu,
        sampler = "gig"
      ),
      raw_results = do.call(rbind, scenario_raw_rows),
      p_values = p_values,
      observed_statistics = observed_statistics,
      warnings = unique(warnings_all)
    )
    saveRDS(scenario_result, output_rds)

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      n = n,
      kappa = kappa,
      M = M,
      B = B,
      n_cores = n_cores,
      rejection_rate_5pct = mean(p_values <= 0.05),
      rejection_rate_10pct = mean(p_values <= 0.10),
      mean_p_value = mean(p_values),
      p_value_quantiles = paste(
        sprintf(
          c("q05=%.6f", "q25=%.6f", "q50=%.6f", "q75=%.6f", "q95=%.6f"),
          p_value_quantiles
        ),
        collapse = "; "
      ),
      elapsed_seconds = elapsed_seconds,
      output_rds = output_rds,
      stringsAsFactors = FALSE
    )

    log_lines <- c(
      log_lines,
      sprintf("scenario_%d_elapsed_seconds: %.3f", scenario_index, elapsed_seconds),
      format_hvmf_gig_warning_block(warnings_all)
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  raw_df <- do.call(rbind, raw_rows)
  summary_csv <- file.path(output_dir, "hvmf_simple_gig_cvm_calibration_summary.csv")
  raw_csv <- file.path(output_dir, "hvmf_simple_gig_cvm_calibration_raw.csv")
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)
  utils::write.csv(raw_df, raw_csv, row.names = FALSE)

  log_lines <- c(
    log_lines,
    sprintf("end_time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  log_txt <- file.path(output_dir, "hvmf_simple_gig_cvm_calibration.log")
  writeLines(log_lines, con = log_txt)

  list(
    summary = summary_df,
    raw_results = raw_df,
    summary_csv = summary_csv,
    raw_csv = raw_csv,
    log_txt = log_txt,
    output_dir = output_dir
  )
}