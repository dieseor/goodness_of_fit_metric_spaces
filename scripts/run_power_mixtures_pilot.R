resolve_power_mixtures_path <- function(...) {
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

calibration_study_script_power_mixtures <- resolve_power_mixtures_path(
  "bootstrap",
  "calibration_study.R"
)
source(calibration_study_script_power_mixtures)

logistic_gaussian_screening_utils_script_power_mixtures <- resolve_power_mixtures_path(
  "real_data",
  "logistic_gaussian",
  "utils_logistic_gaussian_screening.R"
)
source(logistic_gaussian_screening_utils_script_power_mixtures)

logistic_gaussian_distance_profile_script_power_mixtures <- resolve_power_mixtures_path(
  "distance_profiles",
  "logistic_gaussian_distance_profile_analysis.R"
)

load_named_functions_from_file <- function(path, function_names, envir = parent.frame()) {
  expressions <- parse(file = path, keep.source = FALSE)
  found <- character(0)

  for (expr in expressions) {
    if (is.call(expr) &&
        identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) &&
        as.character(expr[[2L]]) %in% function_names) {
      eval(expr, envir = envir)
      found <- c(found, as.character(expr[[2L]]))
    }
  }

  missing_names <- setdiff(function_names, unique(found))
  if (length(missing_names) > 0L) {
    stop(sprintf(
      "Could not load the following functions from '%s': %s",
      path,
      paste(missing_names, collapse = ", ")
    ))
  }

  invisible(unique(found))
}

mvrnorm <- MASS::mvrnorm
load_named_functions_from_file(
  path = logistic_gaussian_distance_profile_script_power_mixtures,
  function_names = c("softmax", "generate_logistic_gaussian_samples"),
  envir = environment()
)

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

parse_named_args_power_mixtures <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  if (length(args) == 0L) {
    return(list())
  }

  output <- vector("list", length(args))
  names(output) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    output[[i]] <- value
    names(output)[[i]] <- key
  }

  output
}

parse_integer_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) {
    return(as.integer(default))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
}

parse_numeric_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) {
    return(as.numeric(default))
  }
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

parse_character_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) {
    return(as.character(default))
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

parse_block_filter <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return("all")
  }

  value <- tolower(trimws(as.character(value)))
  if (!value %in% c("all", "logistic_gaussian", "vmf")) {
    stop("`block` must be one of 'all', 'logistic_gaussian', or 'vmf'.")
  }

  value
}

parse_logical_flag <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(isTRUE(default))
  }

  value_chr <- tolower(trimws(as.character(value)))
  if (value_chr %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value_chr %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop(sprintf("Could not parse logical flag from '%s'.", value))
}

safe_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

format_elapsed_compact <- function(seconds) {
  seconds <- max(0, as.numeric(seconds))
  if (!is.finite(seconds)) {
    return("NA")
  }

  total_seconds <- as.integer(round(seconds))
  hours <- total_seconds %/% 3600L
  minutes <- (total_seconds %% 3600L) %/% 60L
  secs <- total_seconds %% 60L

  if (hours > 0L) {
    sprintf("%02d:%02d:%02d", hours, minutes, secs)
  } else {
    sprintf("%02d:%02d", minutes, secs)
  }
}

update_console_progress <- function(progress_bar,
                                    completed,
                                    total,
                                    started_at_proc,
                                    prefix = "Progress") {
  if (!inherits(progress_bar, "txtProgressBar")) {
    return(invisible(NULL))
  }

  utils::setTxtProgressBar(progress_bar, completed)
  if (completed >= total) {
    cat("\n")
  }
  flush.console()
  invisible(NULL)
}

seed_from_job <- function(base_seed, design_id, rep, stream = 0L) {
  modulus <- 2147483647
  seed <- (as.numeric(base_seed) +
    1000003 * as.numeric(design_id) +
    1009 * as.numeric(rep) +
    10000019 * as.numeric(stream)) %% modulus
  seed <- as.integer(seed)
  if (seed <= 0L) {
    seed <- seed + 1L
  }
  seed
}

rdirichlet_symmetric <- function(n, a) {
  if (length(a) != 1L || !is.finite(a) || a <= 0) {
    stop("`a` must be a strictly positive finite scalar.")
  }
  gamma_draws <- matrix(stats::rgamma(3L * n, shape = a, rate = 1), nrow = n, ncol = 3L)
  gamma_draws / rowSums(gamma_draws)
}

generate_vmf_antipodal_sample <- function(n, beta, mu, kappa) {
  # Here beta is the mixture weight of the alternative component.
  choose_alt <- stats::runif(n) < beta
  sample <- matrix(0, nrow = n, ncol = length(mu))
  n_alt <- sum(choose_alt)
  n_null <- n - n_alt

  if (n_null > 0L) {
    sample[!choose_alt, ] <- rotasym::r_vMF(n = n_null, mu = mu, kappa = kappa)
  }
  if (n_alt > 0L) {
    sample[choose_alt, ] <- rotasym::r_vMF(n = n_alt, mu = -mu, kappa = kappa)
  }

  sample
}

generate_vmf_nonantipodal_sample <- function(n, beta, mu, nu, kappa) {
  # Here beta is the mixture weight of the alternative component.
  choose_alt <- stats::runif(n) < beta
  sample <- matrix(0, nrow = n, ncol = length(mu))
  n_alt <- sum(choose_alt)
  n_null <- n - n_alt

  if (n_null > 0L) {
    sample[!choose_alt, ] <- rotasym::r_vMF(n = n_null, mu = mu, kappa = kappa)
  }
  if (n_alt > 0L) {
    sample[choose_alt, ] <- rotasym::r_vMF(n = n_alt, mu = nu, kappa = kappa)
  }

  sample
}

generate_logistic_dirichlet_sample <- function(n, beta, dirichlet_alpha) {
  choose_alt <- stats::runif(n) < beta
  sample <- matrix(0, nrow = n, ncol = 3L)
  n_alt <- sum(choose_alt)
  n_null <- n - n_alt

  if (n_null > 0L) {
    sample[!choose_alt, ] <- generate_logistic_gaussian_samples(
      n = n_null,
      mu = c(0, 0, 0),
      Sigma = diag(3L)
    )
  }
  if (n_alt > 0L) {
    alpha_vec <- as.numeric(dirichlet_alpha)
    gamma_draws <- matrix(
      stats::rgamma(n_alt * length(alpha_vec), shape = rep(alpha_vec, each = n_alt), rate = 1),
      nrow = n_alt,
      ncol = length(alpha_vec)
    )
    sample[choose_alt, ] <- gamma_draws / rowSums(gamma_draws)
  }

  sample
}

make_vmf_ks_grid <- function(kappa = 2) {
  make_sample_unique_distance_ks_grid()
}

make_exchangeable_logistic_ks_grid <- function() {
  make_sample_unique_distance_ks_grid()
}

filter_design_by_block <- function(design, block) {
  block <- parse_block_filter(block)

  if (identical(block, "all")) {
    return(design)
  }
  if (identical(block, "logistic_gaussian")) {
    return(design[design$scenario == "logistic_gaussian_simplex_d3", , drop = FALSE])
  }

  design[design$scenario != "logistic_gaussian_simplex_d3", , drop = FALSE]
}

make_design_grid <- function(n_values, beta_values) {
  mu <- c(1, 0, 0)
  kappa <- 2
  gamma_values <- c(120)
  dirichlet_cases <- list(
    c(0.2, 0.2, 0.2),
    c(1, 1, 8)
  )

  rows <- list()
  idx <- 1L

  for (n in n_values) {
    for (beta in beta_values) {
      rows[[idx]] <- data.frame(
        scenario = "vmf_s2_antipodal",
        alternative = "antipodal_mixture",
        n = as.integer(n),
        beta = as.numeric(beta),
        gamma_deg = NA_real_,
        dirichlet_a = NA_real_,
        dirichlet_alpha1 = NA_real_,
        dirichlet_alpha2 = NA_real_,
        dirichlet_alpha3 = NA_real_,
        mu_x = mu[[1L]],
        mu_y = mu[[2L]],
        mu_z = mu[[3L]],
        kappa = kappa,
        nu_x = NA_real_,
        nu_y = NA_real_,
        nu_z = NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  for (gamma_deg in gamma_values) {
    gamma_rad <- gamma_deg * pi / 180
    nu <- c(cos(gamma_rad), sin(gamma_rad), 0)
    for (n in n_values) {
      for (beta in beta_values) {
        rows[[idx]] <- data.frame(
          scenario = "vmf_s2_non_antipodal",
          alternative = "non_antipodal_mixture",
          n = as.integer(n),
          beta = as.numeric(beta),
          gamma_deg = as.numeric(gamma_deg),
          dirichlet_a = NA_real_,
          dirichlet_alpha1 = NA_real_,
          dirichlet_alpha2 = NA_real_,
          dirichlet_alpha3 = NA_real_,
          mu_x = mu[[1L]],
          mu_y = mu[[2L]],
          mu_z = mu[[3L]],
          kappa = kappa,
          nu_x = nu[[1L]],
          nu_y = nu[[2L]],
          nu_z = nu[[3L]],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  for (dirichlet_alpha in dirichlet_cases) {
    for (n in n_values) {
      for (beta in beta_values) {
        rows[[idx]] <- data.frame(
          scenario = "logistic_gaussian_simplex_d3",
          alternative = "dirichlet_mixture",
          n = as.integer(n),
          beta = as.numeric(beta),
          gamma_deg = NA_real_,
          dirichlet_a = if (length(unique(dirichlet_alpha)) == 1L) as.numeric(dirichlet_alpha[[1L]]) else NA_real_,
          dirichlet_alpha1 = as.numeric(dirichlet_alpha[[1L]]),
          dirichlet_alpha2 = as.numeric(dirichlet_alpha[[2L]]),
          dirichlet_alpha3 = as.numeric(dirichlet_alpha[[3L]]),
          mu_x = NA_real_,
          mu_y = NA_real_,
          mu_z = NA_real_,
          kappa = NA_real_,
          nu_x = NA_real_,
          nu_y = NA_real_,
          nu_z = NA_real_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  design <- do.call(rbind, rows)
  design$design_id <- seq_len(nrow(design))
  design
}

run_single_power_job <- function(job_row,
                                 design_row,
                                 B,
                                 alpha_nominal,
                                 statistics,
                                 bootstrap_method,
                                 bootstrap_n_cores,
                                 base_seed,
                                 vmf_ks_grid,
                                 logistic_ks_grid,
                                 derivative_mc_size,
                                 fast_multiplier_cvm_block_size,
                                 logistic_gaussian_quadform_method) {
  rep_id <- as.integer(job_row$rep)
  design_id <- as.integer(design_row$design_id)
  data_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 0L)
  bootstrap_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 1L)
  derivative_seed <- seed_from_job(base_seed, design_id = design_id, rep = rep_id, stream = 2L)

  start_time <- proc.time()[["elapsed"]]

  out <- data.frame(
    scenario = design_row$scenario,
    alternative = design_row$alternative,
    n = as.integer(design_row$n),
    beta = as.numeric(design_row$beta),
    gamma_deg = as.numeric(design_row$gamma_deg %||% NA_real_),
    dirichlet_a = as.numeric(design_row$dirichlet_a %||% NA_real_),
    dirichlet_alpha1 = as.numeric(design_row$dirichlet_alpha1 %||% NA_real_),
    dirichlet_alpha2 = as.numeric(design_row$dirichlet_alpha2 %||% NA_real_),
    dirichlet_alpha3 = as.numeric(design_row$dirichlet_alpha3 %||% NA_real_),
    rep = rep_id,
    ks_stat = NA_real_,
    cvm_stat = NA_real_,
    ks_pvalue = NA_real_,
    cvm_pvalue = NA_real_,
    status = "ok",
    error_message = NA_character_,
    bootstrap_method_requested = bootstrap_method,
    bootstrap_method_effective = NA_character_,
    fallback_to_reestimated = NA,
    seed_data = data_seed,
    seed_bootstrap = bootstrap_seed,
    seed_derivative = derivative_seed,
    elapsed_seconds = NA_real_,
    stringsAsFactors = FALSE
  )

  result <- tryCatch(
    {
      set.seed(data_seed)

      if (identical(design_row$scenario, "vmf_s2_antipodal")) {
        x <- generate_vmf_antipodal_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          mu = c(design_row$mu_x, design_row$mu_y, design_row$mu_z),
          kappa = as.numeric(design_row$kappa)
        )

        multiplier_bootstrap_vmf(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = vmf_ks_grid,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            vmf_profile_method = "tabulated",
            vmf_profile_n_u = 4097L
          ),
          distance_type = "geodesic",
          unknown_param = "xi"
        )
      } else if (identical(design_row$scenario, "vmf_s2_non_antipodal")) {
        x <- generate_vmf_nonantipodal_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          mu = c(design_row$mu_x, design_row$mu_y, design_row$mu_z),
          nu = c(design_row$nu_x, design_row$nu_y, design_row$nu_z),
          kappa = as.numeric(design_row$kappa)
        )

        multiplier_bootstrap_vmf(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = vmf_ks_grid,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            vmf_profile_method = "tabulated",
            vmf_profile_n_u = 4097L
          ),
          distance_type = "geodesic",
          unknown_param = "xi"
        )
      } else if (identical(design_row$scenario, "logistic_gaussian_simplex_d3")) {
        x <- generate_logistic_dirichlet_sample(
          n = as.integer(design_row$n),
          beta = as.numeric(design_row$beta),
          dirichlet_alpha = c(
            as.numeric(design_row$dirichlet_alpha1),
            as.numeric(design_row$dirichlet_alpha2),
            as.numeric(design_row$dirichlet_alpha3)
          )
        )

        logistic_ks_grid_design <- make_exchangeable_logistic_ks_grid()

        multiplier_bootstrap_logistic_gaussian(
          data = x,
          null = list(type = "composite"),
          statistics = statistics,
          ks_grid = logistic_ks_grid_design,
          B = as.integer(B),
          alpha = alpha_nominal,
          n_cores = as.integer(bootstrap_n_cores),
          seed = bootstrap_seed,
          bootstrap_method = bootstrap_method,
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = as.integer(derivative_mc_size),
            derivative_mc_seed = as.integer(derivative_seed),
            fast_multiplier_cvm_block_size = as.integer(fast_multiplier_cvm_block_size),
            mvnormal_quadform_method = logistic_gaussian_quadform_method
          ),
          unknown_param = "both"
        )
      } else {
        stop(sprintf("Unsupported scenario '%s'.", design_row$scenario))
      }
    },
    error = function(e) e
  )

  out$elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  if (inherits(result, "error")) {
    out$status <- "error"
    out$error_message <- conditionMessage(result)
    out$bootstrap_method_effective <- NA_character_
    out$fallback_to_reestimated <- NA
    return(out)
  }

  out$ks_stat <- as.numeric(result$observed$ks$statistic %||% NA_real_)
  out$cvm_stat <- as.numeric(result$observed$cvm$statistic %||% NA_real_)
  out$ks_pvalue <- as.numeric(result$inference$ks$p_value %||% NA_real_)
  out$cvm_pvalue <- as.numeric(result$inference$cvm$p_value %||% NA_real_)
  out$bootstrap_method_effective <- as.character(result$diagnostics$effective_bootstrap_method %||% NA_character_)
  out$fallback_to_reestimated <- isTRUE(result$diagnostics$fallback_to_reestimated)
  out
}

summarize_power_results <- function(raw_results, alpha = 0.05) {
  ok <- raw_results[raw_results$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) {
    return(data.frame())
  }

  grouping_key <- paste(
    ok$scenario,
    ok$alternative,
    ok$n,
    ok$beta,
    ifelse(is.na(ok$gamma_deg), "NA", format(ok$gamma_deg, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha1), "NA", format(ok$dirichlet_alpha1, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha2), "NA", format(ok$dirichlet_alpha2, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha3), "NA", format(ok$dirichlet_alpha3, scientific = FALSE, trim = TRUE)),
    sep = "|"
  )

  split_rows <- split(ok, grouping_key, drop = TRUE)
  summary_rows <- lapply(split_rows, function(df) {
    first <- df[1L, , drop = FALSE]
    data.frame(
      scenario = first$scenario,
      alternative = first$alternative,
      n = as.integer(first$n),
      beta = as.numeric(first$beta),
      gamma_deg = as.numeric(first$gamma_deg),
      dirichlet_a = as.numeric(first$dirichlet_a),
      dirichlet_alpha1 = as.numeric(first$dirichlet_alpha1),
      dirichlet_alpha2 = as.numeric(first$dirichlet_alpha2),
      dirichlet_alpha3 = as.numeric(first$dirichlet_alpha3),
      n_success = nrow(df),
      n_fail = sum(raw_results$scenario == first$scenario &
        raw_results$alternative == first$alternative &
        raw_results$n == first$n &
        raw_results$beta == first$beta &
        ((is.na(raw_results$gamma_deg) & is.na(first$gamma_deg)) |
          (!is.na(raw_results$gamma_deg) & !is.na(first$gamma_deg) & raw_results$gamma_deg == first$gamma_deg)) &
        ((is.na(raw_results$dirichlet_alpha1) & is.na(first$dirichlet_alpha1)) |
          (!is.na(raw_results$dirichlet_alpha1) & !is.na(first$dirichlet_alpha1) & raw_results$dirichlet_alpha1 == first$dirichlet_alpha1)) &
        ((is.na(raw_results$dirichlet_alpha2) & is.na(first$dirichlet_alpha2)) |
          (!is.na(raw_results$dirichlet_alpha2) & !is.na(first$dirichlet_alpha2) & raw_results$dirichlet_alpha2 == first$dirichlet_alpha2)) &
        ((is.na(raw_results$dirichlet_alpha3) & is.na(first$dirichlet_alpha3)) |
          (!is.na(raw_results$dirichlet_alpha3) & !is.na(first$dirichlet_alpha3) & raw_results$dirichlet_alpha3 == first$dirichlet_alpha3)) &
        raw_results$status != "ok"),
      power_ks_005 = mean(df$ks_pvalue <= alpha, na.rm = TRUE),
      power_cvm_005 = mean(df$cvm_pvalue <= alpha, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)
  summary_df[order(
    summary_df$scenario,
    summary_df$alternative,
    summary_df$gamma_deg,
    summary_df$dirichlet_alpha1,
    summary_df$dirichlet_alpha2,
    summary_df$dirichlet_alpha3,
    summary_df$n,
    summary_df$beta
  ), , drop = FALSE]
}

summarize_alpha_curves <- function(raw_results, alpha_grid) {
  ok <- raw_results[raw_results$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) {
    return(data.frame())
  }

  grouping_key <- paste(
    ok$scenario,
    ok$alternative,
    ok$n,
    ok$beta,
    ifelse(is.na(ok$gamma_deg), "NA", format(ok$gamma_deg, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha1), "NA", format(ok$dirichlet_alpha1, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha2), "NA", format(ok$dirichlet_alpha2, scientific = FALSE, trim = TRUE)),
    ifelse(is.na(ok$dirichlet_alpha3), "NA", format(ok$dirichlet_alpha3, scientific = FALSE, trim = TRUE)),
    sep = "|"
  )

  split_rows <- split(ok, grouping_key, drop = TRUE)
  curve_rows <- vector("list", length(split_rows) * length(alpha_grid))
  idx <- 1L

  for (df in split_rows) {
    first <- df[1L, , drop = FALSE]
    for (alpha_value in alpha_grid) {
      curve_rows[[idx]] <- data.frame(
        scenario = first$scenario,
        alternative = first$alternative,
        n = as.integer(first$n),
        beta = as.numeric(first$beta),
        gamma_deg = as.numeric(first$gamma_deg),
        dirichlet_a = as.numeric(first$dirichlet_a),
        dirichlet_alpha1 = as.numeric(first$dirichlet_alpha1),
        dirichlet_alpha2 = as.numeric(first$dirichlet_alpha2),
        dirichlet_alpha3 = as.numeric(first$dirichlet_alpha3),
        alpha = as.numeric(alpha_value),
        rejection_prob_ks = mean(df$ks_pvalue <= alpha_value, na.rm = TRUE),
        rejection_prob_cvm = mean(df$cvm_pvalue <= alpha_value, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  curve_df <- do.call(rbind, curve_rows)
  curve_df[order(
    curve_df$scenario,
    curve_df$alternative,
    curve_df$gamma_deg,
    curve_df$dirichlet_alpha1,
    curve_df$dirichlet_alpha2,
    curve_df$dirichlet_alpha3,
    curve_df$n,
    curve_df$beta,
    curve_df$alpha
  ), , drop = FALSE]
}

plot_title_suffix <- function(df) {
  first <- df[1L, , drop = FALSE]

  if (identical(first$scenario, "vmf_s2_antipodal")) {
    return("vMF on S^2: antipodal mixture")
  }
  if (identical(first$scenario, "vmf_s2_non_antipodal")) {
    return(sprintf("vMF on S^2: non-antipodal mixture, gamma = %s deg", first$gamma_deg))
  }
  if (identical(first$scenario, "logistic_gaussian_simplex_d3")) {
    return(sprintf(
      "Simplex D = 3: LG vs Dirichlet (%s, %s, %s)",
      format(first$dirichlet_alpha1, trim = TRUE),
      format(first$dirichlet_alpha2, trim = TRUE),
      format(first$dirichlet_alpha3, trim = TRUE)
    ))
  }

  sprintf("%s / %s", first$scenario, first$alternative)
}

plot_group_slug <- function(df) {
  first <- df[1L, , drop = FALSE]
  pieces <- c(first$scenario, first$alternative)

  if (!is.na(first$gamma_deg)) {
    pieces <- c(pieces, sprintf("gamma_%s", first$gamma_deg))
  }
  if (!is.na(first$dirichlet_a)) {
    pieces <- c(pieces, sprintf("a_%s", first$dirichlet_a))
  } else if (!is.na(first$dirichlet_alpha1) && !is.na(first$dirichlet_alpha2) && !is.na(first$dirichlet_alpha3)) {
    pieces <- c(
      pieces,
      sprintf(
        "alpha_%s_%s_%s",
        format(first$dirichlet_alpha1, trim = TRUE),
        format(first$dirichlet_alpha2, trim = TRUE),
        format(first$dirichlet_alpha3, trim = TRUE)
      )
    )
  }

  safe_slug(paste(pieces, collapse = "_"))
}

save_power_plot <- function(summary_subset, stat_name, file_path) {
  n_values <- sort(unique(summary_subset$n))
  beta_values <- sort(unique(summary_subset$beta))
  colors <- viridisLite::viridis(length(n_values), option = "D", begin = 0.15, end = 0.9)
  y_column <- if (identical(stat_name, "ks")) "power_ks_005" else "power_cvm_005"
  calibration_band <- c(0.0365, 0.0635)

  grDevices::png(filename = file_path, width = 1800, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(6.2, 8.2, 2.2, 1.2))

  plot(
    NA,
    xlim = range(beta_values),
    ylim = c(0, 1),
    xlab = "",
    ylab = expression(P[M](p <= 0.05)),
    cex.axis = 1.5,
    cex.lab = 1.5
  )
  usr <- graphics::par("usr")
  x_pad <- 0.01 * diff(usr[1:2])
  grid(col = "grey78")
  graphics::rect(
    xleft = usr[[1L]] - x_pad,
    ybottom = calibration_band[[1L]],
    xright = usr[[2L]] + x_pad,
    ytop = calibration_band[[2L]],
    col = grDevices::adjustcolor("grey75", alpha.f = 0.35),
    border = NA
  )
  graphics::abline(h = 0.05, lty = 2, lwd = 2, col = "black")
  graphics::mtext(expression(beta), side = 1, line = 2.5, cex = 2)

  for (i in seq_along(n_values)) {
    n_i <- n_values[[i]]
    df_i <- summary_subset[summary_subset$n == n_i, , drop = FALSE]
    df_i <- df_i[order(df_i$beta), , drop = FALSE]
    lines(df_i$beta, df_i[[y_column]], col = colors[[i]], lwd = 3)
    points(df_i$beta, df_i[[y_column]], col = colors[[i]], pch = 19, cex = 1.1)
  }

  legend(
    "topleft",
    legend = sprintf("n = %d", n_values),
    col = colors,
    lwd = 3,
    pch = 19,
    bty = "n",
    cex = 2,
    pt.cex = 2
  )
}

save_alpha_curve_plot <- function(alpha_subset, stat_name, file_path) {
  n_values <- sort(unique(alpha_subset$n))
  colors <- viridisLite::viridis(length(n_values), option = "D", begin = 0.15, end = 0.9)
  y_column <- if (identical(stat_name, "ks")) "rejection_prob_ks" else "rejection_prob_cvm"
  beta_value <- unique(alpha_subset$beta)

  grDevices::png(filename = file_path, width = 1800, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(6.2, 8.2, 2.2, 1.2))

  plot(
    c(0, 1),
    c(0, 1),
    type = "n",
    xlab = "",
    ylab = expression(P[M](p <= alpha)),
    cex.axis = 1.5,
    cex.lab = 1.5
  )
  grid(col = "grey85")
  abline(0, 1, lty = 2, lwd = 2, col = "black")
  graphics::mtext(expression(alpha), side = 1, line = 2.5, cex = 2)

  for (i in seq_along(n_values)) {
    n_i <- n_values[[i]]
    df_i <- alpha_subset[alpha_subset$n == n_i, , drop = FALSE]
    df_i <- df_i[order(df_i$alpha), , drop = FALSE]
    lines(df_i$alpha, df_i[[y_column]], col = colors[[i]], lwd = 3)
  }

  legend(
    "bottomright",
    legend = sprintf("n = %d", n_values),
    col = colors,
    lwd = 3,
    bty = "n",
    cex = 2
  )
}

save_all_plots <- function(summary_df, alpha_curve_df, output_dir) {
  plots_dir <- file.path(output_dir, "plots")
  power_dir <- file.path(plots_dir, "power_vs_beta")
  alpha_dir <- file.path(plots_dir, "alpha_curves")

  dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(summary_df) > 0L) {
    group_key <- paste(
      summary_df$scenario,
      summary_df$alternative,
      ifelse(is.na(summary_df$gamma_deg), "NA", format(summary_df$gamma_deg, trim = TRUE)),
      ifelse(is.na(summary_df$dirichlet_a), "NA", format(summary_df$dirichlet_a, trim = TRUE)),
      sep = "|"
    )
    summary_groups <- split(summary_df, group_key, drop = TRUE)

    for (df in summary_groups) {
      group_slug <- plot_group_slug(df)
      save_power_plot(df, stat_name = "ks", file_path = file.path(power_dir, sprintf("%s_power_ks.png", group_slug)))
      save_power_plot(df, stat_name = "cvm", file_path = file.path(power_dir, sprintf("%s_power_cvm.png", group_slug)))
    }
  }

  if (nrow(alpha_curve_df) > 0L) {
    group_key <- paste(
      alpha_curve_df$scenario,
      alpha_curve_df$alternative,
      ifelse(is.na(alpha_curve_df$gamma_deg), "NA", format(alpha_curve_df$gamma_deg, trim = TRUE)),
      ifelse(is.na(alpha_curve_df$dirichlet_a), "NA", format(alpha_curve_df$dirichlet_a, trim = TRUE)),
      alpha_curve_df$beta,
      sep = "|"
    )
    alpha_groups <- split(alpha_curve_df, group_key, drop = TRUE)

    for (df in alpha_groups) {
      group_slug <- sprintf("%s_beta_%s", plot_group_slug(df), safe_slug(format(df$beta[[1L]], trim = TRUE)))
      save_alpha_curve_plot(df, stat_name = "ks", file_path = file.path(alpha_dir, sprintf("%s_alpha_ks.png", group_slug)))
      save_alpha_curve_plot(df, stat_name = "cvm", file_path = file.path(alpha_dir, sprintf("%s_alpha_cvm.png", group_slug)))
    }
  }
}

write_metadata_files <- function(output_dir,
                                 M,
                                 B,
                                 n_values,
                                 beta_values,
                                 statistics,
                                 bootstrap_method,
                                 bootstrap_n_cores,
                                 n_cores_outer,
                                 base_seed,
                                 derivative_mc_size,
                                 fast_multiplier_cvm_block_size,
                                 logistic_gaussian_quadform_method,
                                 alpha_grid,
                                 block) {
  metadata_lines <- c(
    sprintf("created_at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("M: %d", as.integer(M)),
    sprintf("B: %d", as.integer(B)),
    sprintf("n_values: %s", paste(as.integer(n_values), collapse = ", ")),
    sprintf("beta_values: %s", paste(format(as.numeric(beta_values), trim = TRUE), collapse = ", ")),
    sprintf("statistics: %s", paste(statistics, collapse = ", ")),
    sprintf("bootstrap_method_requested: %s", bootstrap_method),
    sprintf("bootstrap_n_cores: %d", as.integer(bootstrap_n_cores)),
    sprintf("n_cores_outer: %d", as.integer(n_cores_outer)),
    sprintf("block: %s", block),
    sprintf("base_seed: %d", as.integer(base_seed)),
    sprintf("derivative_mc_size: %d", as.integer(derivative_mc_size)),
    sprintf("fast_multiplier_cvm_block_size: %d", as.integer(fast_multiplier_cvm_block_size)),
    sprintf("logistic_gaussian_quadform_method: %s", logistic_gaussian_quadform_method),
    sprintf("alpha_grid: %s", paste(format(alpha_grid, trim = TRUE), collapse = ", ")),
    "logistic_gaussian_null_definition: the LG component is X = softmax(Y) with Y ~ N_3(0, I_3); the matching ilr parametrization used to build the KS grid is N_2(0, I_2)."
  )

  writeLines(metadata_lines, con = file.path(output_dir, "metadata.txt"))
  utils::capture.output(
    print(utils::sessionInfo()),
    file = file.path(output_dir, "sessionInfo.txt")
  )
}

run_jobs_with_progress <- function(total_jobs,
                                   run_job_by_index,
                                   n_cores_outer,
                                   show_progress = TRUE) {
  if (length(total_jobs) != 1L || !is.finite(total_jobs) || total_jobs <= 0L) {
    stop("`total_jobs` must be a strictly positive integer.")
  }

  total_jobs <- as.integer(total_jobs)
  n_cores_outer <- as.integer(n_cores_outer)
  started_at_proc <- proc.time()[["elapsed"]]
  progress_bar <- NULL

  if (isTRUE(show_progress)) {
    progress_bar <- utils::txtProgressBar(min = 0, max = total_jobs, style = 3)
    on.exit(close(progress_bar), add = TRUE)
  }

  if (n_cores_outer <= 1L || .Platform$OS.type == "windows") {
    results <- vector("list", total_jobs)
    for (i in seq_len(total_jobs)) {
      results[[i]] <- run_job_by_index(i)
      if (isTRUE(show_progress)) {
        update_console_progress(
          progress_bar = progress_bar,
          completed = i,
          total = total_jobs,
          started_at_proc = started_at_proc,
          prefix = "jobs"
        )
      }
    }
    return(results)
  }

  results <- vector("list", total_jobs)
  active_jobs <- list()
  active_index <- list()
  launched <- 0L
  completed <- 0L

  launch_next_job <- function(job_index) {
    child <- parallel::mcparallel(
      expr = try(run_job_by_index(job_index), silent = TRUE),
      silent = TRUE
    )
    pid_key <- as.character(child$pid)
    active_jobs[[pid_key]] <<- child
    active_index[[pid_key]] <<- job_index
  }

  while (completed < total_jobs) {
    while (launched < total_jobs && length(active_jobs) < n_cores_outer) {
      launched <- launched + 1L
      launch_next_job(launched)
    }

    collected <- parallel::mccollect(active_jobs, wait = FALSE, timeout = 0.2)
    if (is.null(collected) || length(collected) == 0L) {
      Sys.sleep(0.1)
      next
    }

    for (pid_key in names(collected)) {
      job_index <- active_index[[pid_key]]
      job_result <- collected[[pid_key]]

      if (inherits(job_result, "try-error")) {
        results[[job_index]] <- data.frame(
          scenario = NA_character_,
          alternative = NA_character_,
          n = NA_integer_,
          beta = NA_real_,
          gamma_deg = NA_real_,
          dirichlet_a = NA_real_,
          rep = NA_integer_,
          ks_stat = NA_real_,
          cvm_stat = NA_real_,
          ks_pvalue = NA_real_,
          cvm_pvalue = NA_real_,
          status = "error",
          error_message = as.character(job_result),
          bootstrap_method_requested = NA_character_,
          bootstrap_method_effective = NA_character_,
          fallback_to_reestimated = NA,
          seed_data = NA_integer_,
          seed_bootstrap = NA_integer_,
          seed_derivative = NA_integer_,
          elapsed_seconds = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        results[[job_index]] <- job_result
      }

      active_jobs[[pid_key]] <- NULL
      active_index[[pid_key]] <- NULL
      completed <- completed + 1L

      if (isTRUE(show_progress)) {
        update_console_progress(
          progress_bar = progress_bar,
          completed = completed,
          total = total_jobs,
          started_at_proc = started_at_proc,
          prefix = "jobs"
        )
      }
    }
  }

  results
}

replace_scenario_rows <- function(existing_df,
                                  new_df,
                                  scenario_ids,
                                  order_cols = NULL) {
  if (length(scenario_ids) == 0L) {
    return(new_df)
  }

  if (is.null(existing_df) || nrow(existing_df) == 0L) {
    combined <- new_df
  } else {
    missing_in_existing <- setdiff(names(new_df), names(existing_df))
    for (col_name in missing_in_existing) {
      existing_df[[col_name]] <- NA
    }

    missing_in_new <- setdiff(names(existing_df), names(new_df))
    for (col_name in missing_in_new) {
      new_df[[col_name]] <- NA
    }

    existing_df <- existing_df[, names(new_df), drop = FALSE]
    keep_existing <- !existing_df$scenario %in% scenario_ids
    combined <- rbind(existing_df[keep_existing, , drop = FALSE], new_df)
  }

  if (!is.null(order_cols)) {
    present_order_cols <- order_cols[order_cols %in% names(combined)]
    if (length(present_order_cols) > 0L) {
      ordering_args <- c(combined[present_order_cols], list(na.last = TRUE))
      combined <- combined[do.call(order, ordering_args), , drop = FALSE]
    }
  }

  rownames(combined) <- NULL
  combined
}

merge_with_existing_outputs <- function(output_dir,
                                        raw_results,
                                        summary_df,
                                        alpha_curve_df,
                                        scenario_ids,
                                        block) {
  if (identical(block, "all")) {
    return(list(
      raw_results = raw_results,
      summary = summary_df,
      alpha_curves = alpha_curve_df
    ))
  }

  raw_path <- file.path(output_dir, "power_replications_long.csv")
  summary_path <- file.path(output_dir, "power_summary_005.csv")
  alpha_path <- file.path(output_dir, "alpha_curve_summary.csv")

  existing_raw <- if (file.exists(raw_path)) utils::read.csv(raw_path, stringsAsFactors = FALSE) else NULL
  existing_summary <- if (file.exists(summary_path)) utils::read.csv(summary_path, stringsAsFactors = FALSE) else NULL
  existing_alpha <- if (file.exists(alpha_path)) utils::read.csv(alpha_path, stringsAsFactors = FALSE) else NULL

  list(
    raw_results = replace_scenario_rows(
      existing_df = existing_raw,
      new_df = raw_results,
      scenario_ids = scenario_ids,
      order_cols = c("scenario", "alternative", "gamma_deg", "dirichlet_alpha1", "dirichlet_alpha2", "dirichlet_alpha3", "n", "beta", "rep")
    ),
    summary = replace_scenario_rows(
      existing_df = existing_summary,
      new_df = summary_df,
      scenario_ids = scenario_ids,
      order_cols = c("scenario", "alternative", "gamma_deg", "dirichlet_alpha1", "dirichlet_alpha2", "dirichlet_alpha3", "n", "beta")
    ),
    alpha_curves = replace_scenario_rows(
      existing_df = existing_alpha,
      new_df = alpha_curve_df,
      scenario_ids = scenario_ids,
      order_cols = c("scenario", "alternative", "gamma_deg", "dirichlet_alpha1", "dirichlet_alpha2", "dirichlet_alpha3", "n", "beta", "alpha")
    )
  )
}

run_power_mixtures_pilot <- function(output_dir = file.path("simulation_results", "power_mixtures_pilot"),
                                     M = 250L,
                                     B = 250L,
                                     n_values = c(50L, 100L, 200L),
                                     beta_values = c(0, 0.25, 0.5, 1),
                                     statistics = c("ks", "cvm"),
                                     alpha_nominal = 0.05,
                                     alpha_grid = seq(0, 1, by = 0.01),
                                     bootstrap_method = "fast_multiplier",
                                     bootstrap_n_cores = 1L,
                                     n_cores_outer = 1L,
                                     base_seed = 20260617L,
                                     derivative_mc_size = 1000L,
                                     fast_multiplier_cvm_block_size = 50L,
                                     logistic_gaussian_quadform_method = "auto",
                                     block = "all",
                                     show_progress = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  design <- make_design_grid(
    n_values = as.integer(n_values),
    beta_values = as.numeric(beta_values)
  )
  design <- filter_design_by_block(design, block = block)
  if (nrow(design) == 0L) {
    stop("The selected `block` produced an empty design.")
  }

  jobs <- expand.grid(
    design_id = design$design_id,
    rep = seq_len(as.integer(M)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  vmf_ks_grid <- make_vmf_ks_grid(kappa = 2)
  logistic_ks_grid <- NULL

  if (isTRUE(show_progress)) {
    cat(sprintf(
      "Running %d Monte Carlo jobs with up to %d outer core(s) for block '%s'.\n",
      nrow(jobs),
      max(1L, as.integer(n_cores_outer)),
      block
    ))
  }

  run_job_by_index <- function(i) {
    job_row <- jobs[i, , drop = FALSE]
    design_row <- design[design$design_id == job_row$design_id, , drop = FALSE]
    run_single_power_job(
      job_row = job_row,
      design_row = design_row,
      B = B,
      alpha_nominal = alpha_nominal,
      statistics = statistics,
      bootstrap_method = bootstrap_method,
      bootstrap_n_cores = bootstrap_n_cores,
      base_seed = base_seed,
      vmf_ks_grid = vmf_ks_grid,
      logistic_ks_grid = logistic_ks_grid,
      derivative_mc_size = derivative_mc_size,
      fast_multiplier_cvm_block_size = fast_multiplier_cvm_block_size,
      logistic_gaussian_quadform_method = logistic_gaussian_quadform_method
    )
  }

  raw_rows <- run_jobs_with_progress(
    total_jobs = nrow(jobs),
    run_job_by_index = run_job_by_index,
    n_cores_outer = n_cores_outer,
    show_progress = show_progress
  )

  raw_results <- do.call(rbind, raw_rows)
  raw_results <- raw_results[order(
    raw_results$scenario,
    raw_results$alternative,
    raw_results$gamma_deg,
    raw_results$dirichlet_a,
    raw_results$n,
    raw_results$beta,
    raw_results$rep
  ), , drop = FALSE]

  summary_df <- summarize_power_results(raw_results, alpha = alpha_nominal)
  alpha_curve_df <- summarize_alpha_curves(raw_results, alpha_grid = alpha_grid)
  merged_outputs <- merge_with_existing_outputs(
    output_dir = output_dir,
    raw_results = raw_results,
    summary_df = summary_df,
    alpha_curve_df = alpha_curve_df,
    scenario_ids = unique(design$scenario),
    block = block
  )
  raw_results_out <- merged_outputs$raw_results
  summary_df_out <- merged_outputs$summary
  alpha_curve_df_out <- merged_outputs$alpha_curves

  utils::write.csv(raw_results_out, file.path(output_dir, "power_replications_long.csv"), row.names = FALSE)
  utils::write.csv(summary_df_out, file.path(output_dir, "power_summary_005.csv"), row.names = FALSE)
  utils::write.csv(alpha_curve_df_out, file.path(output_dir, "alpha_curve_summary.csv"), row.names = FALSE)

  save_all_plots(summary_df_out, alpha_curve_df_out, output_dir = output_dir)
  write_metadata_files(
    output_dir = output_dir,
    M = M,
    B = B,
    n_values = n_values,
    beta_values = beta_values,
    statistics = statistics,
    bootstrap_method = bootstrap_method,
    bootstrap_n_cores = bootstrap_n_cores,
    n_cores_outer = n_cores_outer,
    base_seed = base_seed,
    derivative_mc_size = derivative_mc_size,
    fast_multiplier_cvm_block_size = fast_multiplier_cvm_block_size,
    logistic_gaussian_quadform_method = logistic_gaussian_quadform_method,
    alpha_grid = alpha_grid,
    block = block
  )

  invisible(list(
    raw_results = raw_results_out,
    summary = summary_df_out,
    alpha_curves = alpha_curve_df_out
  ))
}

if (sys.nframe() == 0L) {
  args <- parse_named_args_power_mixtures(commandArgs(trailingOnly = TRUE))

  run_power_mixtures_pilot(
    output_dir = as.character(args$output_dir %||% file.path("simulation_results", "power_mixtures_pilot")),
    M = as.integer(args$M %||% 250L),
    B = as.integer(args$B %||% 250L),
    n_values = parse_integer_csv(args$n_values, c(50L, 100L, 200L)),
    beta_values = parse_numeric_csv(args$beta_values, c(0, 0.25, 0.5, 1)),
    statistics = parse_character_csv(args$statistics, c("ks", "cvm")),
    alpha_nominal = as.numeric(args$alpha_nominal %||% 0.05),
    alpha_grid = parse_numeric_csv(args$alpha_grid, seq(0, 1, by = 0.01)),
    bootstrap_method = as.character(args$bootstrap_method %||% "fast_multiplier"),
    bootstrap_n_cores = as.integer(args$bootstrap_n_cores %||% 1L),
    n_cores_outer = as.integer(args$n_cores_outer %||% 1L),
    base_seed = as.integer(args$seed %||% 20260617L),
    derivative_mc_size = as.integer(args$derivative_mc_size %||% 1000L),
    fast_multiplier_cvm_block_size = as.integer(args$fast_multiplier_cvm_block_size %||% 50L),
    logistic_gaussian_quadform_method = as.character(args$logistic_gaussian_quadform_method %||% "auto"),
    block = parse_block_filter(args$block %||% "all"),
    show_progress = parse_logical_flag(args$show_progress, default = TRUE)
  )
}
