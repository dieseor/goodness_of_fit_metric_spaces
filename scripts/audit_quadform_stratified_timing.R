#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): bounded, exploratory timing/accuracy audit only.
# It does not alter the production dispatcher, Section 6, or paper results.
# The available designs are deterministic stress grids; calculation is capped
# at two cores regardless of the caller's request.

Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

`%||%` <- function(x, y) if (is.null(x)) y else x

resolve_path <- function(...) {
  candidates <- c(file.path(...), file.path("..", ...), file.path("..", "..", ...))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Cannot resolve path: ", file.path(...))
  hit[[1L]]
}

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) > 1L) paste(bits[-1L], collapse = "=") else "true"
  }
  out
}

named <- parse_args(args)
output_dir <- named$output_dir %||% file.path(
  "real_data", "bootstrap_audit", "quadform_method_stress_20260724", "stratified_48_cases"
)
n_cores <- min(2L, max(1L, as.integer(named$n_cores %||% 2L)))
points_per_stratum <- max(1L, as.integer(named$points_per_stratum %||% 1L))
design <- tolower(named$design %||% "stratified")
rule_set <- tolower(named$rule_set %||% "historical")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(resolve_path("tests", "benchmark_logistic_gaussian_quadform_utils.R"))

all_cases <- make_stress_quadform_cases()
select_evenly <- function(d, n) {
  if (n > nrow(d)) {
    stop("Requested more points than are available in a stratum.")
  }
  d[unique(round(seq(1L, nrow(d), length.out = n))), , drop = FALSE]
}

if (identical(design, "stratified")) {
  main <- subset(all_cases, q_ratio %in% c(0.1, 1))
  main_groups <- split(main, interaction(main$lambda_family, main$delta_family, main$q_ratio, drop = TRUE))
  main_cases <- do.call(rbind, lapply(main_groups, select_evenly, n = points_per_stratum))

  tail_pool <- subset(
    all_cases,
    lambda_family == "geometric_severe" &
      delta_family %in% c("central", "large_random", "tiny_eig_concentrated", "diffuse_extreme") &
      q_ratio %in% c(0.02, 4)
  )
  tail_groups <- split(tail_pool, interaction(tail_pool$delta_family, tail_pool$q_ratio, drop = TRUE))
  tail_cases <- do.call(rbind, lapply(tail_groups, select_evenly, n = points_per_stratum))
  cases <- rbind(main_cases, tail_cases)
  stopifnot(nrow(cases) == 48L * points_per_stratum)
} else if (identical(design, "preventive_432")) {
  # EN DUDA (2026-07-25): 400 regularly spaced cells plus 32 additional
  # intermediate-q cells in the two adverse spectra. This is three times the
  # previous 144-cell battery, not a probability sample from fitted models.
  regular <- subset(all_cases, q_ratio %in% c(0.02, 0.1, 1, 4))
  targeted <- subset(
    all_cases,
    lambda_family %in% c("geometric_severe", "one_spike") &
      delta_family %in% c("central", "large_random", "tiny_eig_concentrated", "diffuse_extreme") &
      q_ratio %in% c(0.25, 2) & dim %in% c(3L, 10L)
  )
  cases <- rbind(regular, targeted)
  stopifnot(nrow(regular) == 400L, nrow(targeted) == 32L, nrow(cases) == 432L)
} else {
  stop("`design` must be 'stratified' or 'preventive_432'.")
}

analytic_features <- function(row) {
  lambda <- parse_case_vector(row$lambda_json)
  delta <- parse_case_vector(row$delta_json)
  beta <- min(lambda)
  data.frame(
    log_a0_mode1 = 0.5 * sum(log(beta / lambda)) - 0.5 * sum(delta),
    rho_mode1 = 1 - beta / max(lambda),
    log_double_xmin = log(.Machine$double.xmin),
    stringsAsFactors = FALSE
  )
}

features <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) analytic_features(cases[i, , drop = FALSE])))
cases <- cbind(cases, features)

methods <- c("farebrother", "davies", "imhof", "sphunif_hbe")
rows <- parallel::mclapply(
  split(cases, seq_len(nrow(cases))),
  function(case) benchmark_quadform_case(case, methods = methods),
  mc.cores = n_cores,
  mc.preschedule = TRUE,
  mc.set.seed = FALSE
)
results <- do.call(rbind, rows)

fare <- subset(results, method == "farebrother", select = c(case_id, elapsed, ifault, status))
names(fare)[names(fare) == "elapsed"] <- "fare_elapsed"
names(fare)[names(fare) == "ifault"] <- "fare_ifault"
names(fare)[names(fare) == "status"] <- "fare_status"
fare$fare_state <- ifelse(fare$fare_status == "ok" & fare$fare_ifault == 0L, "ifault_0", "ifault_nonzero")

imhof <- subset(results, method == "imhof", select = c(case_id, elapsed, status, abserr))
names(imhof)[names(imhof) == "elapsed"] <- "imhof_elapsed"
names(imhof)[names(imhof) == "status"] <- "imhof_status"

wide <- merge(fare, imhof, by = "case_id", all = TRUE)
summary_timing <- do.call(rbind, lapply(split(wide, wide$fare_state), function(d) {
  data.frame(
    fare_state = d$fare_state[[1L]],
    n = nrow(d),
    fare_median_s = median(d$fare_elapsed),
    fare_p95_s = as.numeric(quantile(d$fare_elapsed, 0.95, names = FALSE)),
    fare_max_s = max(d$fare_elapsed),
    imhof_median_s = median(d$imhof_elapsed),
    imhof_p95_s = as.numeric(quantile(d$imhof_elapsed, 0.95, names = FALSE)),
    imhof_max_s = max(d$imhof_elapsed),
    imhof_valid = sum(d$imhof_status == "ok"),
    stringsAsFactors = FALSE
  )
}))

# EN DUDA (2026-07-25): this classifies numerical operating events, not
# mathematical truth. A Farebrother non-zero ifault is a problem. When both
# methods report their diagnostics, a discrepancy above 1e-6 is separately
# treated as a problem. Cells where Imhof does not meet its own 1e-8 absolute
# error diagnostic remain indeterminate rather than being labelled as false
# positives or false negatives.
fare_detail <- subset(
  results,
  method == "farebrother",
  select = c(case_id, value, ifault, status, n_warnings)
)
names(fare_detail)[names(fare_detail) == "value"] <- "fare_value"
names(fare_detail)[names(fare_detail) == "ifault"] <- "fare_ifault_detail"
names(fare_detail)[names(fare_detail) == "status"] <- "fare_status_detail"
names(fare_detail)[names(fare_detail) == "n_warnings"] <- "fare_n_warnings"
imhof_detail <- subset(
  results,
  method == "imhof",
  select = c(case_id, value, abserr, status, n_warnings)
)
names(imhof_detail)[names(imhof_detail) == "value"] <- "imhof_value"
names(imhof_detail)[names(imhof_detail) == "abserr"] <- "imhof_abserr"
names(imhof_detail)[names(imhof_detail) == "status"] <- "imhof_status_detail"
names(imhof_detail)[names(imhof_detail) == "n_warnings"] <- "imhof_n_warnings"
classification <- Reduce(
  function(x, y) merge(x, y, by = "case_id", all.x = TRUE),
  list(cases, fare_detail, imhof_detail)
)
classification$fare_ok <- classification$fare_status_detail == "ok" &
  classification$fare_ifault_detail == 0L & classification$fare_n_warnings == 0L
classification$imhof_precise <- classification$imhof_status_detail == "ok" &
  is.finite(classification$imhof_value) & is.finite(classification$imhof_abserr) &
  classification$imhof_abserr >= 0 & classification$imhof_abserr <= 1e-8 &
  classification$imhof_n_warnings == 0L
classification$abs_fare_imhof <- abs(classification$fare_value - classification$imhof_value)
classification$fare_ifault_problem <- !classification$fare_ok
classification$reliable_disagreement_problem <- classification$fare_ok &
  classification$imhof_precise & classification$abs_fare_imhof > 1e-6
classification$operating_class <- ifelse(
  classification$fare_ifault_problem | classification$reliable_disagreement_problem,
  "problem", ifelse(
    classification$fare_ok & classification$imhof_precise,
    "good", "indeterminate"
  )
)

if (identical(rule_set, "analytic")) {
  # EN DUDA (2026-07-25): log(a0) is the mode-1 Farebrother leading
  # coefficient on the log scale. The threshold is a conservative floating-
  # point underflow warning; it is not yet a theorem certifying failure.
  rule_matrix <- data.frame(
    case_id = classification$case_id,
    cond_gt_1e4 = classification$condition_number > 1e4,
    log_a0_lt_log_double_xmin = classification$log_a0_mode1 < classification$log_double_xmin,
    stringsAsFactors = FALSE
  )
  rule_matrix$cond_or_log_a0 <- with(
    rule_matrix,
    cond_gt_1e4 | log_a0_lt_log_double_xmin
  )
} else if (identical(rule_set, "historical")) {
  rule_matrix <- data.frame(
    case_id = classification$case_id,
    cond_gt_1e4 = classification$condition_number > 1e4,
    dim5_delta10_qgt_01 = classification$dim >= 5L & classification$delta_max > 10 & classification$q_ratio > 0.1,
    delta1e3_qgt_1 = classification$delta_max > 1e3 & classification$q_ratio > 1,
    stringsAsFactors = FALSE
  )
  rule_matrix$historical_auto <- with(
    rule_matrix,
    cond_gt_1e4 | dim5_delta10_qgt_01 | delta1e3_qgt_1
  )
  rule_matrix$cond_or_dim_delta <- with(rule_matrix, cond_gt_1e4 | dim5_delta10_qgt_01)
  rule_matrix$cond_or_extreme <- with(rule_matrix, cond_gt_1e4 | delta1e3_qgt_1)
} else {
  stop("`rule_set` must be 'analytic' or 'historical'.")
}

score_rule <- function(predicted_imhof, rule_name) {
  d <- classification
  # User-specified definition: a problem is exactly Farebrother ifault != 0.
  # This deliberately makes no appeal to Imhof or to Farebrother--Imhof
  # agreement when calculating the preventive-rule confusion matrix.
  known_problem <- !is.na(d$fare_ifault_detail) & d$fare_ifault_detail != 0L
  known_good <- !is.na(d$fare_ifault_detail) & d$fare_ifault_detail == 0L
  data.frame(
    rule = rule_name,
    flagged = sum(predicted_imhof),
    known_problem = sum(known_problem),
    true_positive = sum(predicted_imhof & known_problem),
    missed_problem = sum(!predicted_imhof & known_problem),
    sensitivity_known_problem = if (any(known_problem)) sum(predicted_imhof & known_problem) / sum(known_problem) else NA_real_,
    known_good = sum(known_good),
    false_positive = sum(predicted_imhof & known_good),
    false_positive_rate_known_good = if (any(known_good)) sum(predicted_imhof & known_good) / sum(known_good) else NA_real_,
    flagged_indeterminate = sum(predicted_imhof & is.na(d$fare_ifault_detail)),
    stringsAsFactors = FALSE
  )
}

rule_scores <- do.call(rbind, lapply(names(rule_matrix)[-1L], function(name) {
  score_rule(rule_matrix[[name]], name)
}))

increment_over_cond <- do.call(rbind, lapply(names(rule_matrix)[-c(1L, 2L)], function(name) {
  extra <- rule_matrix[[name]] & !rule_matrix$cond_gt_1e4
  data.frame(
    rule = name,
    additional_flagged_beyond_cond = sum(extra),
    additional_problem_detected = sum(extra & classification$fare_ifault_detail != 0L, na.rm = TRUE),
    additional_good_false_positive = sum(extra & classification$fare_ifault_detail == 0L, na.rm = TRUE),
    additional_indeterminate = sum(extra & is.na(classification$fare_ifault_detail)),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(cases, file.path(output_dir, "quadform_cases.csv"), row.names = FALSE)
utils::write.csv(results, file.path(output_dir, "quadform_results.csv"), row.names = FALSE)
utils::write.csv(wide, file.path(output_dir, "farebrother_imhof_timing_by_case.csv"), row.names = FALSE)
utils::write.csv(summary_timing, file.path(output_dir, "timing_by_farebrother_status.csv"), row.names = FALSE)
utils::write.csv(classification, file.path(output_dir, "farebrother_operating_classification.csv"), row.names = FALSE)
utils::write.csv(rule_matrix, file.path(output_dir, "preventive_rule_flags.csv"), row.names = FALSE)
utils::write.csv(rule_scores, file.path(output_dir, "preventive_rule_confusion_summary.csv"), row.names = FALSE)
utils::write.csv(increment_over_cond, file.path(output_dir, "preventive_rules_increment_over_cond.csv"), row.names = FALSE)
saveRDS(
  list(cases = cases, results = results, timing = summary_timing,
       classification = classification, rule_scores = rule_scores,
       increment_over_cond = increment_over_cond),
  file.path(output_dir, "stratified_timing_audit.rds")
)

cat("Completed", nrow(cases), "stratified cases with", n_cores, "cores.\n")
