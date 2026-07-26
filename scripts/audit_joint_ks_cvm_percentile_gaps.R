#!/usr/bin/env Rscript

# EN DUDA (2026-07-24): exploratory diagnostic script.  It does not establish
# a production bootstrap or quadrature rule, and its conclusions require
# mathematical review and explicit approval before any use in the paper.

# Examine the *joint* bootstrap reference distribution of the KS and CvM
# percentile ranks.  This is a rank/copula diagnostic only: it reads existing
# bootstrap draws and does not simulate any new samples.

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else file.path(
  "real_data", "bootstrap_audit", "ks_cvm_real_data_20260723"
)
paper_dir <- file.path(
  "real_data", "logistic_gaussian", "screening", "fast",
  "paper_table_B5000_sampleks_fast_rerun_20260718"
)
matched_path <- file.path(output_dir, "matched_fast_vs_reestimated_coxite_wind.rds")

stopifnot(dir.exists(paper_dir), file.exists(matched_path))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

survival_ranks <- function(x) {
  x <- as.numeric(x)
  (1 + vapply(x, function(value) sum(x >= value), integer(1))) / (length(x) + 1)
}

joint_gap_summary <- function(raw_result, case, source, method) {
  ks_star <- as.numeric(raw_result$bootstrap$statistics$ks)
  cvm_star <- as.numeric(raw_result$bootstrap$statistics$cvm)
  stopifnot(length(ks_star) == length(cvm_star), length(ks_star) > 0L)

  ks_rank <- survival_ranks(ks_star)
  cvm_rank <- survival_ranks(cvm_star)
  gap_star <- cvm_rank - ks_rank
  observed_ks_p <- as.numeric(raw_result$inference$ks$p_value)
  observed_cvm_p <- as.numeric(raw_result$inference$cvm$p_value)
  observed_gap <- observed_cvm_p - observed_ks_p
  observed_ks <- as.numeric(raw_result$inference$ks$observed)
  observed_cvm <- as.numeric(raw_result$inference$cvm$observed)
  shape_ratio_star <- cvm_star / ks_star^2

  data.frame(
    case = case,
    source = source,
    method = method,
    B = length(ks_star),
    observed_p_ks = observed_ks_p,
    observed_p_cvm = observed_cvm_p,
    observed_p_cvm_minus_ks = observed_gap,
    observed_ks = observed_ks,
    observed_cvm = observed_cvm,
    observed_shape_ratio_cvm_over_ks_sq = observed_cvm / observed_ks^2,
    bootstrap_rank_correlation = stats::cor(ks_rank, cvm_rank),
    bootstrap_probability_p_cvm_gt_p_ks = mean(cvm_rank > ks_rank),
    bootstrap_probability_gap_at_least_observed = mean(gap_star >= observed_gap),
    bootstrap_gap_q95 = as.numeric(stats::quantile(gap_star, 0.95, names = FALSE)),
    bootstrap_gap_q975 = as.numeric(stats::quantile(gap_star, 0.975, names = FALSE)),
    bootstrap_probability_shape_ratio_at_most_observed = mean(shape_ratio_star <= observed_cvm / observed_ks^2),
    stringsAsFactors = FALSE
  )
}

paper_paths <- sort(Sys.glob(file.path(paper_dir, "*_results.rds")))
paper_rows <- lapply(paper_paths, function(path) {
  result <- readRDS(path)
  joint_gap_summary(
    raw_result = result$bootstrap$raw_result,
    case = sub("_results.rds$", "", basename(path)),
    source = "paper_fast_hbe",
    method = "fast_multiplier"
  )
})

matched <- readRDS(matched_path)
matched_rows <- list(
  joint_gap_summary(
    matched$coxite$fast, "coxite", "matched_same_multipliers", "fast_multiplier"
  ),
  joint_gap_summary(
    matched$coxite$reestimated, "coxite", "matched_same_multipliers", "reestimated_multiplier"
  ),
  joint_gap_summary(
    matched$risoe_nov_dec_125m_start4$fast, "risoe_nov_dec_125m_start4", "matched_same_multipliers", "fast_multiplier"
  ),
  joint_gap_summary(
    matched$risoe_nov_dec_125m_start4$reestimated, "risoe_nov_dec_125m_start4", "matched_same_multipliers", "reestimated_multiplier"
  )
)

out <- do.call(rbind, c(paper_rows, matched_rows))
out <- out[order(out$source, out$bootstrap_probability_gap_at_least_observed), , drop = FALSE]
utils::write.csv(out, file.path(output_dir, "joint_ks_cvm_percentile_gap_audit.csv"), row.names = FALSE)

# When a Farebrother profile audit is available, repeat the same joint-rank
# comparison after replacing only the observed statistic by its exact-profile
# value.  The stored multiplier draws are retained: fast multiplier draws are
# computed from weighted indicators and score derivatives, not by evaluating
# the HBE observed profile.
farebrother_audit_paths <- list.files(
  output_dir,
  pattern = "_hbe_profile_audit[.]csv$",
  full.names = TRUE
)
farebrother_rows <- lapply(farebrother_audit_paths, function(path) {
  audit <- utils::read.csv(path, stringsAsFactors = FALSE)
  audit <- audit[audit$status == "ok", , drop = FALSE]
  lapply(seq_len(nrow(audit)), function(i) {
    row <- audit[i, , drop = FALSE]
    slug <- tolower(gsub("[^a-zA-Z0-9]+", "_", row$dataset))
    result_path <- file.path(paper_dir, paste0(slug, "_results.rds"))
    if (!file.exists(result_path)) {
      warning("No stored paper result for exact-profile audit row: ", row$dataset)
      return(NULL)
    }
    raw_result <- readRDS(result_path)$bootstrap$raw_result
    ks_star <- as.numeric(raw_result$bootstrap$statistics$ks)
    cvm_star <- as.numeric(raw_result$bootstrap$statistics$cvm)
    ks_rank <- survival_ranks(ks_star)
    cvm_rank <- survival_ranks(cvm_star)
    gap_star <- cvm_rank - ks_rank
    exact_ks <- as.numeric(row$farebrother_ks)
    exact_cvm <- as.numeric(row$farebrother_cvm)
    exact_p_ks <- (1 + sum(ks_star >= exact_ks)) / (length(ks_star) + 1)
    exact_p_cvm <- (1 + sum(cvm_star >= exact_cvm)) / (length(cvm_star) + 1)
    exact_gap <- exact_p_cvm - exact_p_ks
    data.frame(
      dataset = row$dataset,
      n = row$n,
      hbe_p_ks = row$hbe_p_ks,
      hbe_p_cvm = row$hbe_p_cvm,
      farebrother_p_ks = exact_p_ks,
      farebrother_p_cvm = exact_p_cvm,
      farebrother_p_cvm_minus_ks = exact_gap,
      bootstrap_rank_correlation = stats::cor(ks_rank, cvm_rank),
      bootstrap_probability_p_cvm_gt_p_ks = mean(cvm_rank > ks_rank),
      bootstrap_probability_gap_at_least_farebrother_gap = mean(gap_star >= exact_gap),
      bootstrap_probability_gap_at_most_farebrother_gap = mean(gap_star <= exact_gap),
      farebrother_shape_ratio_cvm_over_ks_sq = exact_cvm / exact_ks^2,
      bootstrap_probability_shape_ratio_at_most_farebrother = mean(cvm_star / ks_star^2 <= exact_cvm / exact_ks^2),
      max_abs_profile_difference = row$max_abs_profile_difference,
      fraction_hbe_zero_farebrother_gt_001 = row$fraction_hbe_zero_farebrother_gt_001,
      source_audit_file = basename(path),
      stringsAsFactors = FALSE
    )
  })
})
farebrother_rows <- Filter(Negate(is.null), unlist(farebrother_rows, recursive = FALSE))
if (length(farebrother_rows) > 0L) {
  farebrother_out <- do.call(rbind, farebrother_rows)
  farebrother_out <- farebrother_out[order(farebrother_out$bootstrap_probability_gap_at_least_farebrother_gap), , drop = FALSE]
  utils::write.csv(
    farebrother_out,
    file.path(output_dir, "farebrother_joint_ks_cvm_percentile_gap_audit.csv"),
    row.names = FALSE
  )
}

paper_only <- out[out$source == "paper_fast_hbe", , drop = FALSE]
message(sprintf(
  "Paper compositional cases: CvM p > KS p in %d/%d; %d have a conditional bootstrap gap-tail probability <= 0.05.",
  sum(paper_only$observed_p_cvm_minus_ks > 0), nrow(paper_only),
  sum(paper_only$bootstrap_probability_gap_at_least_observed <= 0.05)
))
message("Wrote: ", file.path(output_dir, "joint_ks_cvm_percentile_gap_audit.csv"))
