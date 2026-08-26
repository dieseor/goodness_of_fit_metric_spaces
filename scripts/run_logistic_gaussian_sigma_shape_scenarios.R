#!/usr/bin/env Rscript

# Resumable pilot / production runner for the two Logistic-Gaussian scenarios:
#
# Null fitted family:
#   ilr(X) ~ N_d(mu, sigma^2 A_d),
#   A_d = diag(1, ..., 1, 1/(d+1)).
#
# Data-generating null (beta = 0):
#   LG(0, A_d), exactly M2 of Maphosa et al. (2026), expressed in the ilr
#   coordinates used by this repository.
#
# Scenario t4:
#   P_beta = (1-beta) LG(0, A_d) + beta S_{4,d},
# where S_{4,d} is M3 of Maphosa et al. expressed in ilr coordinates.
# By default the t uses scale A_d, reproducing M3 exactly. Set
# --t_standardized=TRUE to use scale ((nu-2)/nu) A_d so its covariance is A_d.
#
# Scenario dirichlet:
#   P_beta = (1-beta) LG(0, A_d) + beta Dir(alpha_d),
#   alpha_{d,k} = c * 10 k / ((d+1)(d+2)), k = 1,...,d+1.
# The default c=1 reproduces the Dirichlet component of M5 exactly. The scalar
# c is exposed only so a pilot can alter concentration without changing the
# relative alpha pattern if the default alternative saturates too early.
#
# No existing runner or model adapter is modified by this file.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <- if (length(fields) == 1L) "TRUE" else paste(fields[-1L], collapse = "=")
  }
  out
}

parse_csv <- function(value, default, integer = FALSE) {
  if (is.null(value) || !nzchar(value)) return(default)
  out <- as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!is.finite(out))) stop("Invalid comma-separated option.")
  if (isTRUE(integer)) {
    if (any(out != as.integer(out))) stop("Expected integer values.")
    out <- as.integer(out)
  }
  out
}

parse_character_csv <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

parse_bool <- function(value, default = FALSE) {
  if (is.null(value)) return(isTRUE(default))
  x <- tolower(trimws(as.character(value)))
  if (x %in% c("true", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "0", "no", "n")) return(FALSE)
  stop(sprintf("Could not parse logical value '%s'.", value))
}

lg_sigma_shape_betas <- c(0, 0.25, 0.5, 1)
lg_sigma_shape_scenarios <- c("t4", "dirichlet")

lg_sigma_shape_seed <- function(base_seed, design_id, replication, stream) {
  as.integer((as.numeric(base_seed) + 1000003 * design_id +
    1009 * replication + 10000019 * stream) %% 2147483647) + 1L
}

maphosa_dirichlet_alpha <- function(d, concentration_multiplier = 1) {
  d <- as.integer(d)
  concentration_multiplier <- as.numeric(concentration_multiplier)
  if (length(d) != 1L || !is.finite(d) || d < 1L) {
    stop("`d` must be a strictly positive integer.")
  }
  if (length(concentration_multiplier) != 1L ||
      !is.finite(concentration_multiplier) || concentration_multiplier <= 0) {
    stop("`concentration_multiplier` must be a strictly positive finite scalar.")
  }
  D <- d + 1L
  concentration_multiplier * 10 * seq_len(D) / (D * (D + 1L))
}

r_maphosa_t_logistic <- function(n, d, nu = 4, standardized = FALSE) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for the t-logistic scenario.")
  }
  if (nu <= 2) stop("`nu` must exceed 2.")
  shape <- logistic_gaussian_maphosa_shape(d)
  scale <- if (isTRUE(standardized)) ((nu - 2) / nu) * shape else shape
  z <- mvtnorm::rmvt(
    n = as.integer(n),
    sigma = scale,
    df = nu,
    delta = rep.int(0, d),
    type = "shifted"
  )
  logistic_gaussian_ilr_to_simplex(z, ambient_dim = d + 1L)
}

r_maphosa_dirichlet <- function(n, d, concentration_multiplier = 1) {
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop(paste(
      "Package 'gtools' is required only for the Dirichlet scenario.",
      "Install it locally with renv::install('gtools') (or install.packages('gtools'))."
    ))
  }
  alpha <- maphosa_dirichlet_alpha(d, concentration_multiplier)
  gtools::rdirichlet(n = as.integer(n), alpha = alpha)
}

generate_lg_sigma_shape_sample <- function(scenario, n, d, beta, nu = 4,
                                            t_standardized = FALSE,
                                            dirichlet_concentration_multiplier = 1) {
  scenario <- match.arg(scenario, lg_sigma_shape_scenarios)
  n <- as.integer(n)
  d <- as.integer(d)
  beta <- as.numeric(beta)
  if (n < 1L || d < 1L || !is.finite(beta) || beta < 0 || beta > 1) {
    stop("Invalid sample-generation settings.")
  }

  x <- rlogistic_gaussian_sigma_shape(
    n = n,
    mu_ilr = rep.int(0, d),
    sigma = 1
  )
  take_h1 <- stats::runif(n) < beta
  if (!any(take_h1)) return(x)

  replacement <- switch(
    scenario,
    t4 = r_maphosa_t_logistic(
      n = sum(take_h1),
      d = d,
      nu = nu,
      standardized = t_standardized
    ),
    dirichlet = r_maphosa_dirichlet(
      n = sum(take_h1),
      d = d,
      concentration_multiplier = dirichlet_concentration_multiplier
    )
  )
  x[take_h1, ] <- replacement
  x
}

lg_sigma_shape_design <- function(dimensions, n_values, beta_values, scenarios, M) {
  dimensions <- sort(unique(as.integer(dimensions)))
  n_values <- sort(unique(as.integer(n_values)))
  beta_values <- sort(unique(as.numeric(beta_values)))
  scenarios <- unique(as.character(scenarios))

  if (!all(beta_values %in% lg_sigma_shape_betas)) {
    stop("beta_values must be a subset of c(0, 0.25, 0.5, 1).")
  }
  if (!all(scenarios %in% lg_sigma_shape_scenarios)) {
    stop("scenarios must be a subset of c('t4', 'dirichlet').")
  }

  design <- expand.grid(
    scenario = scenarios,
    d = dimensions,
    n = n_values,
    beta = beta_values,
    stringsAsFactors = FALSE
  )

  # Stable IDs are assigned from the complete two-scenario beta catalogue.
  catalog <- expand.grid(
    scenario = lg_sigma_shape_scenarios,
    d = dimensions,
    n = n_values,
    beta = lg_sigma_shape_betas,
    stringsAsFactors = FALSE
  )
  catalog$design_id <- seq_len(nrow(catalog))
  design_key <- paste(design$scenario, design$d, design$n, design$beta, sep = "|")
  catalog_key <- paste(catalog$scenario, catalog$d, catalog$n, catalog$beta, sep = "|")
  idx <- match(design_key, catalog_key)
  if (anyNA(idx)) stop("Could not assign stable design_id values.")
  design$design_id <- catalog$design_id[idx]

  merge(design, data.frame(replication = seq_len(as.integer(M))), by = NULL)
}

lg_sigma_shape_key <- function(x) {
  paste(
    x$scenario,
    x$d,
    x$n,
    formatC(x$beta, digits = 16L, format = "fg", flag = "#"),
    x$replication,
    sep = "|"
  )
}

empty_results <- function() {
  data.frame(
    scenario = character(), d = integer(), n = integer(), beta = numeric(),
    design_id = integer(), replication = integer(),
    status = character(), error_message = character(), warning_message = character(),
    seed_data = integer(), seed_bootstrap = integer(), seed_derivative = integer(),
    mu_hat = character(), sigma_hat = numeric(), score_mean_norm = numeric(),
    ks_statistic = numeric(), cvm_statistic = numeric(),
    ks_pvalue = numeric(), cvm_pvalue = numeric(),
    ks_reject = logical(), cvm_reject = logical(),
    effective_bootstrap_method = character(), fallback_to_reestimated = logical(),
    fast_backend = character(), fast_kernel = character(), fast_fused = logical(),
    derivative_method = character(), derivative_mc_size = integer(),
    vhat_method = character(), elapsed_seconds = numeric(),
    stringsAsFactors = FALSE
  )
}

write_atomic_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop(sprintf("Could not update '%s'.", path))
}

summarize_results <- function(x) {
  x <- x[x$status == "ok", , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  groups <- split(x, interaction(x$scenario, x$d, x$n, x$beta, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) data.frame(
    scenario = z$scenario[[1L]], d = z$d[[1L]], n = z$n[[1L]],
    beta = z$beta[[1L]], M = nrow(z),
    ks_rejection_percent = 100 * mean(z$ks_reject),
    cvm_rejection_percent = 100 * mean(z$cvm_reject),
    ks_pvalue_mean = mean(z$ks_pvalue),
    cvm_pvalue_mean = mean(z$cvm_pvalue),
    mean_elapsed_seconds = mean(z$elapsed_seconds),
    stringsAsFactors = FALSE
  )))
}

conforming <- function(x, derivative_mc_size) {
  x$status == "ok" &
    x$effective_bootstrap_method == "fast_multiplier" &
    !x$fallback_to_reestimated &
    x$fast_backend == "cpp" &
    x$fast_kernel == "contiguous_double" &
    x$fast_fused &
    x$derivative_method == "score_mc" &
    x$derivative_mc_size == as.integer(derivative_mc_size) &
    x$vhat_method == "logistic_gaussian_sigma_shape_analytic_fisher"
}

write_status <- function(path, total, results, started, cores, derivative_mc_size) {
  ok <- if (nrow(results)) conforming(results, derivative_mc_size) else logical()
  per_job <- if (any(ok)) mean(results$elapsed_seconds[ok]) else NA_real_
  remaining <- total - sum(ok)
  eta <- if (is.finite(per_job)) per_job * remaining / cores else NA_real_
  writeLines(c(
    sprintf("updated: %s", format(Sys.time(), tz = "Europe/Madrid")),
    sprintf("completed: %d/%d", sum(ok), total),
    sprintf("pending: %d", remaining),
    sprintf("elapsed_seconds: %.1f", as.numeric(difftime(Sys.time(), started, units = "secs"))),
    sprintf("mean_seconds_per_conforming_job: %s", if (is.finite(per_job)) per_job else "NA"),
    sprintf("eta_seconds: %s", if (is.finite(eta)) round(eta) else "NA"),
    sprintf("ok: %d", sum(results$status == "ok")),
    sprintf("errors: %d", sum(results$status == "error")),
    sprintf("nonconforming: %d", sum(results$status == "nonconforming"))
  ), path)
}

run_logistic_gaussian_sigma_shape_scenarios <- function(
    output_dir,
    M = 100L,
    B = 499L,
    dimensions = c(2L, 5L),
    n_values = c(50L, 100L, 200L, 400L),
    beta_values = c(0, 0.25, 0.5, 1),
    scenarios = c("t4", "dirichlet"),
    nu = 4,
    t_standardized = FALSE,
    dirichlet_concentration_multiplier = 1,
    derivative_mc_size = 10000L,
    cvm_block_size = 50L,
    cores = 4L,
    checkpoint_results = 100L,
    base_seed = 20260826L,
    show_progress = TRUE) {

  if (nu <= 2 || any(dimensions < 1L) || any(n_values < 2L) ||
      any(beta_values < 0 | beta_values > 1) || M < 1L || B < 1L ||
      derivative_mc_size < 1L || cvm_block_size < 1L || cores < 1L ||
      checkpoint_results < 1L || dirichlet_concentration_multiplier <= 0) {
    stop("Invalid Logistic-Gaussian sigma-shape scenario settings.")
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Outer parallelism requires a Unix platform.")
  }
  if ("dirichlet" %in% scenarios && !requireNamespace("gtools", quietly = TRUE)) {
    stop(paste(
      "The Dirichlet scenario requires package 'gtools'.",
      "Install it locally with renv::install('gtools') before running this script."
    ))
  }

  source(file.path("bootstrap", "logistic_gaussian_sigma_shape_bootstrap.R"), local = environment())
  design <- lg_sigma_shape_design(dimensions, n_values, beta_values, scenarios, M)
  manifest <- unique(design[c("scenario", "d", "n", "beta", "design_id")])
  manifest <- transform(
    manifest,
    M = as.integer(M), B = as.integer(B), nu = as.numeric(nu),
    t_standardized = isTRUE(t_standardized),
    dirichlet_concentration_multiplier = as.numeric(dirichlet_concentration_multiplier),
    base_seed = as.integer(base_seed),
    derivative_mc_size = as.integer(derivative_mc_size),
    cvm_block_size = as.integer(cvm_block_size),
    null_model = "LG(mu,sigma^2 A_d), A_d=diag(1,...,1,1/(d+1))",
    derivative_method = "score_mc",
    statistics = "ks,cvm",
    ks_grid = "sample_points_unique_distances",
    bootstrap_method = "fast_multiplier",
    fast_backend = "cpp",
    fast_kernel = "contiguous_double",
    fast_fused = TRUE,
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(output_dir, ".logistic_gaussian_sigma_shape.lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Output directory is locked by another run.")
  }
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  manifest_path <- file.path(output_dir, "manifest.csv")
  results_path <- file.path(output_dir, "raw_results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  status_path <- file.path(output_dir, "progress_status.txt")

  if (file.exists(manifest_path)) {
    old <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(old, manifest, check.attributes = FALSE))) {
      stop("Existing manifest is incompatible; use a new output directory.")
    }
  } else {
    write_atomic_csv(manifest, manifest_path)
  }

  results <- if (file.exists(results_path)) {
    utils::read.csv(results_path, stringsAsFactors = FALSE)
  } else {
    empty_results()
  }
  if (!all(names(empty_results()) %in% names(results))) {
    stop("Existing raw results have an incompatible schema.")
  }

  done <- if (nrow(results)) {
    lg_sigma_shape_key(results[conforming(results, derivative_mc_size), , drop = FALSE])
  } else {
    character()
  }
  pending <- design[!lg_sigma_shape_key(design) %in% done, , drop = FALSE]
  started <- Sys.time()

  run_one <- function(row) {
    began <- proc.time()[["elapsed"]]
    data_seed <- lg_sigma_shape_seed(base_seed, row$design_id, row$replication, 0L)
    bootstrap_seed <- lg_sigma_shape_seed(base_seed, row$design_id, row$replication, 1L)
    derivative_seed <- lg_sigma_shape_seed(base_seed, row$design_id, row$replication, 2L)

    base <- data.frame(
      scenario = as.character(row$scenario),
      d = as.integer(row$d),
      n = as.integer(row$n),
      beta = as.numeric(row$beta),
      design_id = as.integer(row$design_id),
      replication = as.integer(row$replication),
      status = "ok",
      error_message = NA_character_,
      warning_message = NA_character_,
      seed_data = data_seed,
      seed_bootstrap = bootstrap_seed,
      seed_derivative = derivative_seed,
      stringsAsFactors = FALSE
    )

    tryCatch({
      set.seed(data_seed)
      x <- generate_lg_sigma_shape_sample(
        scenario = base$scenario,
        n = base$n,
        d = base$d,
        beta = base$beta,
        nu = nu,
        t_standardized = t_standardized,
        dirichlet_concentration_multiplier = dirichlet_concentration_multiplier
      )

      warnings <- character()
      fit <- withCallingHandlers(
        multiplier_bootstrap_logistic_gaussian_sigma_shape(
          data = x,
          null = list(type = "composite"),
          statistics = c("ks", "cvm"),
          ks_grid = make_sample_unique_distance_ks_grid(),
          B = B,
          alpha = 0.05,
          n_cores = 1L,
          seed = bootstrap_seed,
          bootstrap_method = "fast_multiplier",
          keep = list(
            observed_process = FALSE,
            bootstrap_statistics = FALSE,
            bootstrap_thetas = FALSE
          ),
          control = list(
            derivative_method = "score_mc",
            derivative_mc_size = derivative_mc_size,
            derivative_mc_seed = derivative_seed,
            fast_multiplier_cvm_block_size = cvm_block_size,
            fast_multiplier_backend = "cpp",
            fast_multiplier_cpp_kernel = "contiguous_double",
            fast_multiplier_fuse_ks_cvm = TRUE,
            fast_multiplier_cache_corrections = "auto",
            fast_multiplier_stream_chunk_size = 100L
          ),
          distance_profile_backend = "r",
          fast_multiplier_backend = "cpp",
          fast_multiplier_cpp_kernel = "contiguous_double",
          fuse_ks_cvm = TRUE,
          cache_block_corrections = "auto"
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      theta_hat <- fit$observed$theta_hat
      score <- logistic_gaussian_sigma_shape_score_matrix(x, theta_hat)
      diagnostics <- fit$diagnostics
      out <- cbind(base, data.frame(
        mu_hat = paste(theta_hat$mu_ilr, collapse = ";"),
        sigma_hat = theta_hat$sigma,
        score_mean_norm = sqrt(sum(colMeans(score)^2)),
        ks_statistic = fit$inference$ks$observed,
        cvm_statistic = fit$inference$cvm$observed,
        ks_pvalue = fit$inference$ks$p_value,
        cvm_pvalue = fit$inference$cvm$p_value,
        ks_reject = fit$inference$ks$reject,
        cvm_reject = fit$inference$cvm$reject,
        effective_bootstrap_method = diagnostics$effective_bootstrap_method %||% NA_character_,
        fallback_to_reestimated = isTRUE(diagnostics$fallback_to_reestimated),
        fast_backend = diagnostics$fast_multiplier_backend_effective %||% NA_character_,
        fast_kernel = diagnostics$fast_multiplier_cpp_kernel_effective %||% NA_character_,
        fast_fused = isTRUE(diagnostics$fast_multiplier_fuse_ks_cvm_effective),
        derivative_method = diagnostics$derivative_method_effective %||% NA_character_,
        derivative_mc_size = diagnostics$derivative_mc_size %||% NA_integer_,
        vhat_method = diagnostics$vhat_method %||% NA_character_,
        elapsed_seconds = proc.time()[["elapsed"]] - began,
        stringsAsFactors = FALSE
      ))
      out$warning_message <- if (length(warnings)) {
        paste(unique(warnings), collapse = " | ")
      } else {
        NA_character_
      }
      if (!conforming(out, derivative_mc_size)) {
        out$status <- "nonconforming"
        out$error_message <- "Requested Logistic-Gaussian sigma-shape fast configuration was not effective."
      }
      out
    }, error = function(e) {
      base$status <- "error"
      base$error_message <- conditionMessage(e)
      cbind(base, data.frame(
        mu_hat = NA_character_, sigma_hat = NA_real_, score_mean_norm = NA_real_,
        ks_statistic = NA_real_, cvm_statistic = NA_real_,
        ks_pvalue = NA_real_, cvm_pvalue = NA_real_,
        ks_reject = NA, cvm_reject = NA,
        effective_bootstrap_method = NA_character_, fallback_to_reestimated = NA,
        fast_backend = NA_character_, fast_kernel = NA_character_, fast_fused = NA,
        derivative_method = NA_character_, derivative_mc_size = NA_integer_,
        vhat_method = NA_character_, elapsed_seconds = proc.time()[["elapsed"]] - began,
        stringsAsFactors = FALSE
      ))
    })
  }

  total <- nrow(design)
  write_status(status_path, total, results, started, cores, derivative_mc_size)

  while (nrow(pending)) {
    batch_n <- min(as.integer(checkpoint_results), nrow(pending))
    batch <- pending[seq_len(batch_n), , drop = FALSE]
    pending <- pending[-seq_len(batch_n), , drop = FALSE]

    rows <- split(batch, seq_len(nrow(batch)))
    batch_results <- if (cores > 1L) {
      parallel::mclapply(rows, run_one, mc.cores = as.integer(cores), mc.preschedule = FALSE)
    } else {
      lapply(rows, run_one)
    }
    results <- rbind(results, do.call(rbind, batch_results))
    write_atomic_csv(results, results_path)
    write_atomic_csv(summarize_results(results), summary_path)
    write_status(status_path, total, results, started, cores, derivative_mc_size)

    if (isTRUE(show_progress)) {
      ok <- sum(conforming(results, derivative_mc_size))
      cat(sprintf("[%s] %d/%d conforming jobs complete\n",
                  format(Sys.time(), tz = "Europe/Madrid"), ok, total))
      flush.console()
    }
  }

  invisible(list(raw_results = results, summary = summarize_results(results), manifest = manifest))
}

if (sys.nframe() == 0L) {
  cli <- parse_args(commandArgs(trailingOnly = TRUE))
  output_dir <- cli$output_dir %||%
    file.path("simulation_results", "logistic_gaussian_sigma_shape_scenarios")

  run_logistic_gaussian_sigma_shape_scenarios(
    output_dir = output_dir,
    M = as.integer(cli$M %||% 100L),
    B = as.integer(cli$B %||% 499L),
    dimensions = parse_csv(cli$d, c(2L, 5L), integer = TRUE),
    n_values = parse_csv(cli$n, c(50L, 100L, 200L, 400L), integer = TRUE),
    beta_values = parse_csv(cli$beta, c(0, 0.25, 0.5, 1)),
    scenarios = parse_character_csv(cli$scenario, c("t4", "dirichlet")),
    nu = as.numeric(cli$nu %||% 4),
    t_standardized = parse_bool(cli$t_standardized, FALSE),
    dirichlet_concentration_multiplier = as.numeric(
      cli$dirichlet_concentration_multiplier %||% 1
    ),
    derivative_mc_size = as.integer(cli$derivative_mc_size %||% 10000L),
    cvm_block_size = as.integer(cli$cvm_block_size %||% 50L),
    cores = as.integer(cli$cores %||% 4L),
    checkpoint_results = as.integer(cli$checkpoint_results %||% 100L),
    base_seed = as.integer(cli$seed %||% 20260826L),
    show_progress = parse_bool(cli$show_progress, TRUE)
  )
}
