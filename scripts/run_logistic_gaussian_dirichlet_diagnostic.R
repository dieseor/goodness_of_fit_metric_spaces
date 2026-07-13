#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

resolve_diag_path <- function(...) {
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

source(resolve_diag_path("bootstrap", "calibration_study.R"))
source(resolve_diag_path("real_data", "logistic_gaussian", "utils_logistic_gaussian_screening.R"))

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

parse_named_args <- function(args) {
  if (length(args) == 0L) {
    return(list())
  }

  args <- args[startsWith(args, "--")]
  if (length(args) == 0L) {
    return(list())
  }

  out <- vector("list", length(args))
  names(out) <- rep("", length(args))

  for (i in seq_along(args)) {
    arg <- substring(args[[i]], 3L)
    pieces <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- if (length(pieces) > 1L) paste(pieces[-1L], collapse = "=") else "TRUE"
    out[[i]] <- value
    names(out)[[i]] <- key
  }

  out
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

rdirichlet_matrix <- function(n, alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 3L || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop("`alpha` must be a positive vector of length 3.")
  }
  draws <- matrix(stats::rgamma(n * 3L, shape = rep(alpha, each = n), rate = 1), nrow = n, ncol = 3L)
  draws / rowSums(draws)
}

make_dirichlet_design <- function(n_values,
                                  symmetric_a_values,
                                  asymmetric_alphas) {
  rows <- list()
  idx <- 1L

  for (a in symmetric_a_values) {
    for (n in n_values) {
      rows[[idx]] <- data.frame(
        family = "symmetric",
        label = sprintf("Dir(%s,%s,%s)", a, a, a),
        alpha1 = a,
        alpha2 = a,
        alpha3 = a,
        n = as.integer(n),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  for (alpha in asymmetric_alphas) {
    for (n in n_values) {
      rows[[idx]] <- data.frame(
        family = "asymmetric",
        label = sprintf("Dir(%s,%s,%s)", alpha[[1L]], alpha[[2L]], alpha[[3L]]),
        alpha1 = alpha[[1L]],
        alpha2 = alpha[[2L]],
        alpha3 = alpha[[3L]],
        n = as.integer(n),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  design <- do.call(rbind, rows)
  design$design_id <- seq_len(nrow(design))
  design
}

seed_from_job <- function(base_seed, design_id, rep) {
  modulus <- 2147483647
  seed <- (as.numeric(base_seed) + 1000003 * as.numeric(design_id) + 1009 * as.numeric(rep)) %% modulus
  seed <- as.integer(seed)
  if (seed <= 0L) {
    seed <- seed + 1L
  }
  seed
}

run_single_job <- function(job,
                           design,
                           B,
                           ks_grid,
                           base_seed) {
  row <- design[design$design_id == job$design_id, , drop = FALSE]
  rep_id <- as.integer(job$rep)
  seed <- seed_from_job(base_seed, design_id = row$design_id, rep = rep_id)
  set.seed(seed)

  x <- rdirichlet_matrix(
    n = as.integer(row$n),
    alpha = c(row$alpha1, row$alpha2, row$alpha3)
  )

  result <- multiplier_bootstrap_logistic_gaussian(
    data = x,
    null = list(type = "composite"),
    statistics = c("ks", "cvm"),
    ks_grid = ks_grid,
    B = as.integer(B),
    alpha = 0.05,
    n_cores = 1L,
    seed = seed + 1L,
    bootstrap_method = "fast_multiplier",
    keep = list(
      observed_process = FALSE,
      bootstrap_statistics = FALSE,
      bootstrap_thetas = FALSE
    ),
    control = list(
      derivative_method = "score_mc",
      derivative_mc_size = 1000L,
      derivative_mc_seed = seed + 2L,
      fast_multiplier_cvm_block_size = 50L,
      logistic_gaussian_quadform_method = "hbe"
    ),
    unknown_param = "both"
  )

  data.frame(
    family = row$family,
    label = row$label,
    alpha1 = row$alpha1,
    alpha2 = row$alpha2,
    alpha3 = row$alpha3,
    n = as.integer(row$n),
    rep = rep_id,
    ks_pvalue = as.numeric(result$inference$ks$p_value),
    cvm_pvalue = as.numeric(result$inference$cvm$p_value),
    ks_stat = as.numeric(result$observed$ks$statistic),
    cvm_stat = as.numeric(result$observed$cvm$statistic),
    effective_bootstrap_method = as.character(result$diagnostics$effective_bootstrap_method %||% NA_character_),
    fallback_to_reestimated = isTRUE(result$diagnostics$fallback_to_reestimated),
    stringsAsFactors = FALSE
  )
}

run_dirichlet_lg_diagnostic <- function(output_dir = file.path("simulation_results", "logistic_gaussian_dirichlet_diagnostic"),
                                        n_values = c(50L, 100L, 200L),
                                        symmetric_a_values = c(0.2, 0.5, 1, 2, 4, 6, 10),
                                        asymmetric_alphas = list(
                                          c(0.5, 0.5, 2),
                                          c(1, 2, 4),
                                          c(1, 1, 8),
                                          c(2, 3, 9)
                                        ),
                                        M = 80L,
                                        B = 200L,
                                        n_cores = 12L,
                                        seed = 20260618L,
                                        show_progress = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  design <- make_dirichlet_design(
    n_values = as.integer(n_values),
    symmetric_a_values = as.numeric(symmetric_a_values),
    asymmetric_alphas = asymmetric_alphas
  )
  jobs <- expand.grid(
    design_id = design$design_id,
    rep = seq_len(as.integer(M)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  ks_grid <- make_sample_unique_distance_ks_grid()

  if (isTRUE(show_progress)) {
    message(sprintf("Running %d jobs with %d core(s).", nrow(jobs), as.integer(n_cores)))
  }

  results <- parallel::mclapply(
    seq_len(nrow(jobs)),
    function(i) run_single_job(
      job = jobs[i, , drop = FALSE],
      design = design,
      B = B,
      ks_grid = ks_grid,
      base_seed = seed
    ),
    mc.cores = max(1L, as.integer(n_cores))
  )

  raw_df <- do.call(rbind, results)
  summary_key <- paste(raw_df$family, raw_df$label, raw_df$n, sep = "|")
  split_rows <- split(raw_df, summary_key, drop = TRUE)

  summary_df <- do.call(rbind, lapply(split_rows, function(df) {
    first <- df[1L, , drop = FALSE]
    data.frame(
      family = first$family,
      label = first$label,
      alpha1 = first$alpha1,
      alpha2 = first$alpha2,
      alpha3 = first$alpha3,
      n = first$n,
      M = nrow(df),
      power_ks_005 = mean(df$ks_pvalue <= 0.05),
      power_cvm_005 = mean(df$cvm_pvalue <= 0.05),
      median_ks_pvalue = stats::median(df$ks_pvalue),
      median_cvm_pvalue = stats::median(df$cvm_pvalue),
      stringsAsFactors = FALSE
    )
  }))

  raw_df <- raw_df[order(raw_df$family, raw_df$alpha1, raw_df$alpha2, raw_df$alpha3, raw_df$n, raw_df$rep), , drop = FALSE]
  summary_df <- summary_df[order(summary_df$family, summary_df$alpha1, summary_df$alpha2, summary_df$alpha3, summary_df$n), , drop = FALSE]

  utils::write.csv(raw_df, file.path(output_dir, "raw_results.csv"), row.names = FALSE)
  utils::write.csv(summary_df, file.path(output_dir, "summary_results.csv"), row.names = FALSE)

  invisible(list(raw = raw_df, summary = summary_df, output_dir = output_dir))
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  result <- run_dirichlet_lg_diagnostic(
    output_dir = as.character(args$output_dir %||% file.path("simulation_results", "logistic_gaussian_dirichlet_diagnostic")),
    n_values = parse_integer_csv(args$n_values, c(50L, 100L, 200L)),
    symmetric_a_values = parse_numeric_csv(args$symmetric_a_values, c(0.2, 0.5, 1, 2, 4, 6, 10)),
    M = as.integer(args$M %||% 80L),
    B = as.integer(args$B %||% 200L),
    n_cores = as.integer(args$n_cores %||% 12L),
    seed = as.integer(args$seed %||% 20260618L),
    show_progress = parse_logical_flag(args$show_progress, default = TRUE)
  )

  message(sprintf("Summary CSV: %s", file.path(result$output_dir, "summary_results.csv")))
}
