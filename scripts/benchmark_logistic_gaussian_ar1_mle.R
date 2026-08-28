#!/usr/bin/env Rscript

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    fields <- strsplit(substring(arg, 3L), "=", fixed = TRUE)[[1L]]
    out[[fields[[1L]]]] <-
      if (length(fields) == 1L) "TRUE" else paste(fields[-1L], collapse = "=")
  }
  out
}

parse_num_csv <- function(x, default, integer = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  ans <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(!is.finite(ans))) stop("Invalid numeric comma-separated argument.")
  if (integer) {
    if (any(ans != as.integer(ans))) stop("Expected integer values.")
    ans <- as.integer(ans)
  }
  ans
}

source(file.path("bootstrap", "logistic_gaussian_ar1_model_spec.R"))

required <- c(
  "logistic_gaussian_ar1_covariance",
  "logistic_gaussian_ar1_profile_moments",
  "logistic_gaussian_ar1_profile_loglik_from_moments",
  "fit_logistic_gaussian_ar1_rho_optimize",
  "fit_logistic_gaussian_ar1_rho_global"
)

missing <- required[!vapply(required, exists, logical(1), mode = "function")]
if (length(missing)) {
  stop(
    "Missing functions from the AR(1) MLE patch: ",
    paste(missing, collapse = ", ")
  )
}

rho_score_from_moments <- function(rho, moments) {
  d <- moments$d
  A <- moments$A
  B <- moments$B
  C <- moments$C

  numerator <-
    -(d - 1L) * rho^3 +
    B * rho^2 +
    ((d - 1L) - (A + C)) * rho +
    B

  numerator / (1 - rho^2)^2
}

dense_reference <- function(moments,
                            rho_bound = 0.995,
                            grid_size = 20001L) {
  grid <- seq(-rho_bound, rho_bound, length.out = grid_size)

  objective <- vapply(
    grid,
    logistic_gaussian_ar1_profile_loglik_from_moments,
    numeric(1),
    moments = moments
  )

  best_grid <- which.max(objective)

  candidates_rho <- grid[[best_grid]]
  candidates_obj <- objective[[best_grid]]

  if (best_grid > 1L && best_grid < length(grid)) {
    lo <- grid[[best_grid - 1L]]
    hi <- grid[[best_grid + 1L]]

    local <- stats::optimize(
      f = function(rho) {
        logistic_gaussian_ar1_profile_loglik_from_moments(rho, moments)
      },
      interval = c(lo, hi),
      maximum = TRUE,
      tol = 1e-12
    )

    candidates_rho <- c(candidates_rho, local$maximum)
    candidates_obj <- c(candidates_obj, local$objective)
  }

  # Explicitly include the admissible boundaries.
  boundary_rho <- c(-rho_bound, rho_bound)
  boundary_obj <- vapply(
    boundary_rho,
    logistic_gaussian_ar1_profile_loglik_from_moments,
    numeric(1),
    moments = moments
  )

  candidates_rho <- c(candidates_rho, boundary_rho)
  candidates_obj <- c(candidates_obj, boundary_obj)

  j <- which.max(candidates_obj)

  list(
    rho = candidates_rho[[j]],
    objective = candidates_obj[[j]]
  )
}

prepare_case <- function(z, weights = NULL) {
  n <- nrow(z)

  probability_weights <- if (is.null(weights)) {
    rep.int(1 / n, n)
  } else {
    weights / sum(weights)
  }

  mu_hat <- colSums(z * probability_weights)
  centered <- sweep(z, 2L, mu_hat, FUN = "-")

  list(
    centered = centered,
    weights = probability_weights,
    moments = logistic_gaussian_ar1_profile_moments(
      centered,
      probability_weights
    )
  )
}

time_method <- function(cases, method, rho_bound, root_imag_tol) {
  # Untimed warm-up to avoid charging lazy initialization to one method.
  if (identical(method, "old")) {
    invisible(fit_logistic_gaussian_ar1_rho_optimize(
      cases[[1L]]$centered,
      cases[[1L]]$weights,
      rho_bound = rho_bound,
      rho_tol = 1e-12
    ))
  } else {
    invisible(fit_logistic_gaussian_ar1_rho_global(
      cases[[1L]]$centered,
      cases[[1L]]$weights,
      rho_bound = rho_bound,
      root_imag_tol = root_imag_tol
    ))
  }

  t0 <- proc.time()[["elapsed"]]

  fits <- if (identical(method, "old")) {
    lapply(cases, function(case) {
      fit_logistic_gaussian_ar1_rho_optimize(
        case$centered,
        case$weights,
        rho_bound = rho_bound,
        rho_tol = 1e-12
      )
    })
  } else {
    lapply(cases, function(case) {
      fit_logistic_gaussian_ar1_rho_global(
        case$centered,
        case$weights,
        rho_bound = rho_bound,
        root_imag_tol = root_imag_tol
      )
    })
  }

  elapsed <- proc.time()[["elapsed"]] - t0

  list(
    fits = fits,
    elapsed = elapsed
  )
}

run_mode <- function(cases,
                     mode,
                     cell_id,
                     d,
                     n,
                     rho0,
                     rho_bound,
                     root_imag_tol,
                     disagreement_tol) {

  # Alternate timing order across cells/modes to avoid a systematic
  # first-method/cache advantage.
  mode_index <- if (identical(mode, "unweighted")) 0L else 1L
  old_first <- ((cell_id + mode_index) %% 2L) == 0L

  if (old_first) {
    old_run <- time_method(
      cases, "old", rho_bound, root_imag_tol
    )
    new_run <- time_method(
      cases, "new", rho_bound, root_imag_tol
    )
  } else {
    new_run <- time_method(
      cases, "new", rho_bound, root_imag_tol
    )
    old_run <- time_method(
      cases, "old", rho_bound, root_imag_tol
    )
  }

  M <- length(cases)

  rows <- vector("list", M)

  for (m in seq_len(M)) {
    case <- cases[[m]]
    old <- old_run$fits[[m]]
    new <- new_run$fits[[m]]

    rho_old <- old$rho
    rho_new <- new$rho

    # Evaluate both estimates with exactly the SAME closed-form
    # profile likelihood. This avoids comparing objective implementations.
    obj_old <- logistic_gaussian_ar1_profile_loglik_from_moments(
      rho_old,
      case$moments
    )
    obj_new <- logistic_gaussian_ar1_profile_loglik_from_moments(
      rho_new,
      case$moments
    )

    score_old <- rho_score_from_moments(
      rho_old,
      case$moments
    )
    score_new <- rho_score_from_moments(
      rho_new,
      case$moments
    )

    rho_diff <- abs(rho_new - rho_old)
    objective_gain <- obj_new - obj_old

    needs_reference <-
      rho_diff > disagreement_tol ||
      abs(objective_gain) > 1e-10

    reference_rho <- NA_real_
    reference_objective <- NA_real_
    new_reference_gap <- NA_real_
    old_reference_gap <- NA_real_

    if (needs_reference) {
      ref <- dense_reference(
        case$moments,
        rho_bound = rho_bound
      )

      reference_rho <- ref$rho
      reference_objective <- ref$objective
      new_reference_gap <- ref$objective - obj_new
      old_reference_gap <- ref$objective - obj_old
    }

    rows[[m]] <- data.frame(
      cell_id = cell_id,
      d = d,
      n = n,
      rho0 = rho0,
      mode = mode,
      replication = m,

      rho_old = rho_old,
      rho_new = rho_new,
      abs_rho_difference = rho_diff,

      abs_error_old = abs(rho_old - rho0),
      abs_error_new = abs(rho_new - rho0),

      objective_old = obj_old,
      objective_new = obj_new,
      objective_gain_new_minus_old = objective_gain,

      old_reported_objective_error = old$objective - obj_old,

      abs_score_old = abs(score_old),
      abs_score_new = abs(score_new),

      old_at_boundary =
        abs(rho_old) >= rho_bound - 1e-10,
      new_at_boundary =
        abs(rho_new) >= rho_bound - 1e-10,

      needs_reference = needs_reference,
      reference_rho = reference_rho,
      reference_objective = reference_objective,
      new_reference_gap = new_reference_gap,
      old_reference_gap = old_reference_gap,

      stringsAsFactors = FALSE
    )
  }

  raw <- do.call(rbind, rows)

  timing <- data.frame(
    cell_id = cell_id,
    d = d,
    n = n,
    rho0 = rho0,
    mode = mode,
    M = M,
    timing_order = if (old_first) "old_then_new" else "new_then_old",
    old_elapsed_seconds = old_run$elapsed,
    new_elapsed_seconds = new_run$elapsed,
    speedup_old_over_new =
      old_run$elapsed / new_run$elapsed,
    old_fits_per_second =
      M / old_run$elapsed,
    new_fits_per_second =
      M / new_run$elapsed,
    stringsAsFactors = FALSE
  )

  list(raw = raw, timing = timing)
}

summarize_accuracy <- function(raw) {
  groups <- split(
    raw,
    interaction(
      raw$d,
      raw$n,
      raw$rho0,
      raw$mode,
      drop = TRUE
    )
  )

  do.call(rbind, lapply(groups, function(x) {
    data.frame(
      d = x$d[[1L]],
      n = x$n[[1L]],
      rho0 = x$rho0[[1L]],
      mode = x$mode[[1L]],
      M = nrow(x),

      mean_abs_rho_difference =
        mean(x$abs_rho_difference),
      q95_abs_rho_difference =
        unname(quantile(x$abs_rho_difference, 0.95)),
      q99_abs_rho_difference =
        unname(quantile(x$abs_rho_difference, 0.99)),
      max_abs_rho_difference =
        max(x$abs_rho_difference),

      n_rho_diff_gt_1e8 =
        sum(x$abs_rho_difference > 1e-8),
      n_rho_diff_gt_1e6 =
        sum(x$abs_rho_difference > 1e-6),
      n_rho_diff_gt_1e4 =
        sum(x$abs_rho_difference > 1e-4),

      mean_objective_gain =
        mean(x$objective_gain_new_minus_old),
      min_objective_gain =
        min(x$objective_gain_new_minus_old),
      max_objective_gain =
        max(x$objective_gain_new_minus_old),

      n_new_strictly_better =
        sum(x$objective_gain_new_minus_old > 1e-10),
      n_new_strictly_worse =
        sum(x$objective_gain_new_minus_old < -1e-10),

      mean_abs_score_old =
        mean(x$abs_score_old),
      mean_abs_score_new =
        mean(x$abs_score_new),

      q99_abs_score_old =
        unname(quantile(x$abs_score_old, 0.99)),
      q99_abs_score_new =
        unname(quantile(x$abs_score_new, 0.99)),

      max_abs_score_old =
        max(x$abs_score_old),
      max_abs_score_new =
        max(x$abs_score_new),

      mae_rho_old =
        mean(x$abs_error_old),
      mae_rho_new =
        mean(x$abs_error_new),

      rmse_rho_old =
        sqrt(mean((x$rho_old - x$rho0)^2)),
      rmse_rho_new =
        sqrt(mean((x$rho_new - x$rho0)^2)),

      old_boundary_percent =
        100 * mean(x$old_at_boundary),
      new_boundary_percent =
        100 * mean(x$new_at_boundary),

      n_reference_checks =
        sum(x$needs_reference),

      max_new_reference_gap =
        if (any(x$needs_reference))
          max(x$new_reference_gap, na.rm = TRUE)
        else 0,

      max_old_reference_gap =
        if (any(x$needs_reference))
          max(x$old_reference_gap, na.rm = TRUE)
        else 0,

      stringsAsFactors = FALSE
    )
  }))
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))

output_dir <- cli$output_dir
if (is.null(output_dir)) {
  output_dir <- file.path(
    "simulation_results",
    "lg_ar1_mle_old_vs_global"
  )
}

M <- as.integer(if (is.null(cli$M)) 1000L else cli$M)
dimensions <- parse_num_csv(
  cli$d,
  c(2L, 5L),
  integer = TRUE
)
n_values <- parse_num_csv(
  cli$n,
  c(50L, 100L, 200L, 400L, 800L),
  integer = TRUE
)
rho_values <- parse_num_csv(
  cli$rho,
  c(-0.9, -0.5, 0, 0.5, 0.9)
)
cores <- as.integer(if (is.null(cli$cores)) 10L else cli$cores)
base_seed <- as.integer(
  if (is.null(cli$seed)) 20260902L else cli$seed
)

rho_bound <- as.numeric(
  if (is.null(cli$rho_bound)) 0.995 else cli$rho_bound
)
root_imag_tol <- as.numeric(
  if (is.null(cli$root_imag_tol)) 1e-9 else cli$root_imag_tol
)
disagreement_tol <- as.numeric(
  if (is.null(cli$disagreement_tol)) 1e-7 else cli$disagreement_tol
)

if (M < 1L || cores < 1L ||
    any(dimensions < 2L) ||
    any(n_values < 2L) ||
    any(abs(rho_values) >= 1) ||
    rho_bound <= 0 || rho_bound >= 1) {
  stop("Invalid benchmark settings.")
}

if (.Platform$OS.type != "unix" && cores > 1L) {
  stop("Parallel benchmark with cores > 1 requires Unix/macOS/Linux.")
}

design <- expand.grid(
  d = dimensions,
  n = n_values,
  rho0 = rho_values,
  stringsAsFactors = FALSE
)
design$cell_id <- seq_len(nrow(design))

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(
  transform(
    design,
    M = M,
    cores = cores,
    base_seed = base_seed,
    rho_bound = rho_bound
  ),
  file.path(output_dir, "manifest.csv"),
  row.names = FALSE
)

run_cell <- function(i) {
  cell <- design[i, ]

  d <- as.integer(cell$d)
  n <- as.integer(cell$n)
  rho0 <- as.numeric(cell$rho0)
  cell_id <- as.integer(cell$cell_id)

  set.seed(
    as.integer(
      (as.numeric(base_seed) +
         1000003 * cell_id) %% 2147483647
    ) + 1L
  )

  Sigma <- logistic_gaussian_ar1_covariance(
    d,
    rho0
  )
  chol_Sigma <- chol(Sigma)

  # Generate all data once. Both methods see exactly the same samples.
  z_samples <- vector("list", M)
  exp_weights <- vector("list", M)

  for (m in seq_len(M)) {
    standard_normal <- matrix(
      stats::rnorm(n * d),
      nrow = n,
      ncol = d
    )

    z_samples[[m]] <-
      standard_normal %*% chol_Sigma

    # Positive mean-one style multipliers. Normalization is performed
    # inside prepare_case(). This stresses the weighted MLE used by
    # bootstrap re-estimation.
    exp_weights[[m]] <- stats::rexp(n, rate = 1)
  }

  unweighted_cases <- lapply(
    z_samples,
    prepare_case,
    weights = NULL
  )

  weighted_cases <- Map(
    prepare_case,
    z_samples,
    exp_weights
  )

  unweighted <- run_mode(
    cases = unweighted_cases,
    mode = "unweighted",
    cell_id = cell_id,
    d = d,
    n = n,
    rho0 = rho0,
    rho_bound = rho_bound,
    root_imag_tol = root_imag_tol,
    disagreement_tol = disagreement_tol
  )

  weighted <- run_mode(
    cases = weighted_cases,
    mode = "exp_weighted",
    cell_id = cell_id,
    d = d,
    n = n,
    rho0 = rho0,
    rho_bound = rho_bound,
    root_imag_tol = root_imag_tol,
    disagreement_tol = disagreement_tol
  )

  message(sprintf(
    "[cell %02d/%02d] d=%d n=%d rho=%+.1f complete",
    cell_id,
    nrow(design),
    d,
    n,
    rho0
  ))

  list(
    raw = rbind(
      unweighted$raw,
      weighted$raw
    ),
    timing = rbind(
      unweighted$timing,
      weighted$timing
    )
  )
}

started <- Sys.time()

results <- parallel::mclapply(
  seq_len(nrow(design)),
  run_cell,
  mc.cores = min(cores, nrow(design)),
  mc.preschedule = FALSE,
  mc.set.seed = FALSE
)

raw <- do.call(
  rbind,
  lapply(results, `[[`, "raw")
)

timing <- do.call(
  rbind,
  lapply(results, `[[`, "timing")
)

accuracy <- summarize_accuracy(raw)

disagreements <- raw[
  raw$needs_reference |
    raw$abs_rho_difference > 1e-7 |
    abs(raw$objective_gain_new_minus_old) > 1e-10,
  ,
  drop = FALSE
]

utils::write.csv(
  raw,
  file.path(output_dir, "raw_results.csv"),
  row.names = FALSE
)

utils::write.csv(
  timing,
  file.path(output_dir, "timing_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  accuracy,
  file.path(output_dir, "accuracy_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  disagreements,
  file.path(output_dir, "disagreement_cases.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Global summary
# ------------------------------------------------------------

overall_accuracy <- data.frame(
  comparisons = nrow(raw),

  mean_abs_rho_difference =
    mean(raw$abs_rho_difference),
  q99_abs_rho_difference =
    unname(quantile(raw$abs_rho_difference, 0.99)),
  max_abs_rho_difference =
    max(raw$abs_rho_difference),

  n_rho_diff_gt_1e8 =
    sum(raw$abs_rho_difference > 1e-8),
  n_rho_diff_gt_1e6 =
    sum(raw$abs_rho_difference > 1e-6),
  n_rho_diff_gt_1e4 =
    sum(raw$abs_rho_difference > 1e-4),

  min_objective_gain =
    min(raw$objective_gain_new_minus_old),
  max_objective_gain =
    max(raw$objective_gain_new_minus_old),

  n_new_strictly_better =
    sum(raw$objective_gain_new_minus_old > 1e-10),
  n_new_strictly_worse =
    sum(raw$objective_gain_new_minus_old < -1e-10),

  max_abs_score_old =
    max(raw$abs_score_old),
  max_abs_score_new =
    max(raw$abs_score_new),

  n_reference_checks =
    sum(raw$needs_reference),

  max_new_reference_gap =
    if (any(raw$needs_reference))
      max(raw$new_reference_gap, na.rm = TRUE)
    else 0,

  max_old_reference_gap =
    if (any(raw$needs_reference))
      max(raw$old_reference_gap, na.rm = TRUE)
    else 0,

  stringsAsFactors = FALSE
)

overall_timing <- data.frame(
  old_total_seconds =
    sum(timing$old_elapsed_seconds),
  new_total_seconds =
    sum(timing$new_elapsed_seconds),
  overall_speedup_old_over_new =
    sum(timing$old_elapsed_seconds) /
      sum(timing$new_elapsed_seconds),
  stringsAsFactors = FALSE
)

utils::write.csv(
  overall_accuracy,
  file.path(output_dir, "overall_accuracy.csv"),
  row.names = FALSE
)

utils::write.csv(
  overall_timing,
  file.path(output_dir, "overall_timing.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Hard validation checks
# ------------------------------------------------------------

if (any(raw$objective_gain_new_minus_old < -1e-9)) {
  stop(
    "VALIDATION FAILURE: new global MLE gave lower profile likelihood ",
    "than the old optimizer by more than tolerance."
  )
}

interior_new <- !raw$new_at_boundary

if (any(raw$abs_score_new[interior_new] > 1e-7)) {
  stop(
    "VALIDATION FAILURE: an interior global MLE has a non-negligible ",
    "rho-score residual."
  )
}

if (any(raw$needs_reference)) {
  reference_rows <- raw$needs_reference

  if (any(
    raw$new_reference_gap[reference_rows] > 1e-8,
    na.rm = TRUE
  )) {
    stop(
      "VALIDATION FAILURE: dense independent reference found a ",
      "meaningfully better objective than the new method."
    )
  }
}

elapsed <- as.numeric(
  difftime(
    Sys.time(),
    started,
    units = "secs"
  )
)

cat("\n============================================\n")
cat("AR(1) MLE OLD vs GLOBAL benchmark complete\n")
cat("============================================\n")
cat(sprintf("Cells:             %d\n", nrow(design)))
cat(sprintf("Samples per cell:  %d\n", M))
cat(sprintf("Comparisons:       %d\n", nrow(raw)))
cat(sprintf("Cores:             %d\n", cores))
cat(sprintf("Wall time:         %.1f sec\n\n", elapsed))

print(overall_accuracy, row.names = FALSE)
cat("\n")
print(overall_timing, row.names = FALSE)

cat("\nOutput directory:\n")
cat(output_dir, "\n")
