#!/usr/bin/env Rscript

argument <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  values <- commandArgs(trailingOnly = TRUE)
  hit <- values[startsWith(values, prefix)]
  if (!length(hit)) default else substring(hit[[1L]], nchar(prefix) + 1L)
}

quadrature_dir <- argument("quadrature_dir")
score_dir <- argument("score_dir")
output <- argument("output", file.path(quadrature_dir, "paired_comparison.csv"))
if (is.null(quadrature_dir) || is.null(score_dir)) {
  stop("Supply `--quadrature_dir` and `--score_dir`.")
}

read_results <- function(directory, label) {
  path <- file.path(directory, "raw_results.csv")
  if (!file.exists(path)) stop(sprintf("Missing %s.", path))
  value <- utils::read.csv(path, stringsAsFactors = FALSE)
  value <- value[value$status == "ok" & value$beta == 0, , drop = FALSE]
  value$method <- label
  value
}

quadrature <- read_results(quadrature_dir, "quadrature")
score <- read_results(score_dir, "score_mc")
keys <- c("scenario", "family", "d", "n", "beta", "design_id", "rep")
paired <- merge(
  quadrature, score, by = keys, suffixes = c("_quadrature", "_score_mc"),
  all = FALSE, sort = TRUE
)
if (!nrow(paired)) stop("The two runs have no paired successful jobs.")
if (any(paired$seed_data_quadrature != paired$seed_data_score_mc) ||
    any(paired$seed_bootstrap_quadrature != paired$seed_bootstrap_score_mc)) {
  stop("The runs are not paired: data or bootstrap seeds differ.")
}

summarize_statistic <- function(data, statistic) {
  q_reject <- data[[paste0(statistic, "_reject_quadrature")]]
  s_reject <- data[[paste0(statistic, "_reject_score_mc")]]
  q_pvalue <- data[[paste0(statistic, "_pvalue_quadrature")]]
  s_pvalue <- data[[paste0(statistic, "_pvalue_score_mc")]]
  pvalue_difference <- q_pvalue - s_pvalue
  q_only <- sum(q_reject & !s_reject)
  s_only <- sum(!q_reject & s_reject)
  discordant <- q_only + s_only
  data.frame(
    statistic = toupper(statistic),
    M_paired = length(q_reject),
    rejection_quadrature = mean(q_reject),
    rejection_score_mc = mean(s_reject),
    difference = mean(q_reject) - mean(s_reject),
    distance_to_005_quadrature = abs(mean(q_reject) - 0.05),
    distance_to_005_score_mc = abs(mean(s_reject) - 0.05),
    quadrature_only_rejections = q_only,
    score_mc_only_rejections = s_only,
    mean_pvalue_quadrature = mean(q_pvalue),
    mean_pvalue_score_mc = mean(s_pvalue),
    mean_pvalue_difference = mean(pvalue_difference),
    median_pvalue_quadrature = stats::median(q_pvalue),
    median_pvalue_score_mc = stats::median(s_pvalue),
    pvalue_rmse = sqrt(mean(pvalue_difference^2)),
    pvalue_max_abs_difference = max(abs(pvalue_difference)),
    paired_one_sided_p = if (discordant) {
      stats::binom.test(q_only, discordant, p = 0.5,
                        alternative = "greater")$p.value
    } else 1,
    mean_seconds_quadrature = mean(data$elapsed_seconds_quadrature),
    mean_seconds_score_mc = mean(data$elapsed_seconds_score_mc),
    time_ratio = mean(data$elapsed_seconds_quadrature) /
      mean(data$elapsed_seconds_score_mc)
  )
}

groups <- split(paired, interaction(
  paired$family, paired$scenario, paired$d, paired$n, drop = TRUE
))
summary <- do.call(rbind, lapply(groups, function(group) {
  common <- group[1L, c("family", "scenario", "d", "n")]
  rbind(
    cbind(common, summarize_statistic(group, "ks")),
    cbind(common, summarize_statistic(group, "cvm"))
  )
}))
rownames(summary) <- NULL
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, output, row.names = FALSE)
print(summary)
