#!/usr/bin/env zsh

# Parallel validation launcher only; it does not change production runners.

setopt pipefail

repo_root="${0:A:h:h}"
cd "$repo_root" || exit 1
cores="${PROFILE_LOOKUP_PILOT_CORES:-10}"
if [[ ! "$cores" =~ '^[1-9][0-9]*$' ]]; then
  print -u2 -- "PROFILE_LOOKUP_PILOT_CORES must be a positive integer."
  exit 2
fi

out="${1:-benchmarks/profile_lookup_joint_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$out/cases" "$out/status" || exit 1
out="${out:A}"
status_dir="$out/status"

# One BLAS thread per worker makes the ten-worker limit meaningful.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

# n, data seed, B. The n=200 entry reproduces the archived benchmark exactly.
tasks=(
  "50 20260501 2999" "50 20260502 2999" "50 20260503 2999"
  "100 20261001 2999" "100 20261002 2999" "100 20261003 2999"
  "200 2026080101 4999" "200 20262002 2999" "200 20262003 2999"
  "400 20264001 2999" "400 20264002 2999" "400 20264003 2999"
)

total=${#tasks}
next=1
active=0
completed=0
failed=0

render_progress() {
  local width=32 filled empty done_bar left_bar percent
  filled=$(( completed * width / total ))
  empty=$(( width - filled ))
  done_bar=$(printf '%*s' "$filled" '' | tr ' ' '#')
  left_bar=$(printf '%*s' "$empty" '' | tr ' ' '-')
  percent=$(( completed * 100 / total ))
  printf '\r[%s%s] %3d%%  %2d/%2d completed  %2d active  %d failed' \
    "$done_bar" "$left_bar" "$percent" "$completed" "$total" "$active" "$failed"
}

launch_task() {
  local task="$1" n data_seed B multiplier_seed stem log marker
  read -r n data_seed B <<< "$task"
  if [[ "$data_seed" == "2026080101" ]]; then
    multiplier_seed=2026080102
  else
    multiplier_seed=$(( data_seed + 50000 ))
  fi
  stem="n${n}_B${B}_data${data_seed}"
  log="$out/cases/${stem}.log"
  marker="$status_dir/${stem}"
  (
    if Rscript --vanilla scripts/validate_profile_lookup_end_to_end.R \
      --n="$n" --B="$B" --cores=1 \
      --data-seed="$data_seed" --multiplier-seed="$multiplier_seed" \
      --reference-cpp-kernel=legacy \
      --lookup-cpp-kernel=contiguous_double \
      --output="$out/cases/${stem}.csv" >"$log" 2>&1; then
      print -r -- "$stem" > "${marker}.done"
    else
      print -r -- "$stem" > "${marker}.failed"
    fi
  ) &
}

while (( next <= total && active < cores )); do
  launch_task "${tasks[$next]}"
  (( next++ ))
  (( active++ ))
done

print -- "Profile lookup + kernel joint pilot"
print -- "Output: $out"
print -- "Workers: $cores; expected duration: roughly 10 minutes."

while (( completed < total )); do
  for marker in "$status_dir"/*.done(N) "$status_dir"/*.failed(N); do
    [[ -f "$marker" && ! -f "${marker}.reported" ]] || continue
    touch "${marker}.reported"
    (( completed++ ))
    (( active-- ))
    if [[ "$marker" == *.failed ]]; then
      (( failed++ ))
      stem="${${marker:t}%.failed}"
      print -u2 -- "\nERROR: ${stem} failed. Last 25 log lines:"
      tail -n 25 "$out/cases/${stem}.log" >&2
    fi
  done
  while (( next <= total && active < cores )); do
    launch_task "${tasks[$next]}"
    (( next++ ))
    (( active++ ))
  done
  render_progress
  sleep 1
done
print

if (( failed > 0 )); then
  print -u2 -- "Pilot finished with $failed failed task(s). Inspect $out/cases/*.log."
  exit 1
fi

Rscript --vanilla - "$out" <<'RS'
out <- commandArgs(trailingOnly = TRUE)[1L]
files <- list.files(file.path(out, "cases"), pattern = "\\.csv$", full.names = TRUE)
if (length(files) != 12L) stop("Expected 12 completed case CSVs.")
x <- do.call(rbind, lapply(files, utils::read.csv, check.names = FALSE))
utils::write.csv(x, file.path(out, "all_results.csv"), row.names = FALSE)

precision <- aggregate(
  cbind(max_abs_difference, rmse_difference) ~
    model + n + B + comparison_path + quantity + statistic,
  data = x, FUN = max
)
utils::write.csv(precision, file.path(out, "precision_summary.csv"), row.names = FALSE)

kernel <- subset(x, comparison_path == "kernel_only" & quantity == "bootstrap_replicates")
kernel_max <- max(kernel$max_abs_difference)
if (kernel_max > 1e-12) stop(sprintf("Kernel-only comparison exceeded 1e-12: %.17g", kernel_max))

timing_columns <- c(
  "model", "n", "B", "table_build_seconds", "direct_F_seconds",
  "direct_derivative_prep_seconds", "lookup_F_and_derivative_seconds",
  "direct_legacy_kernel_seconds", "lookup_candidate_kernel_seconds"
)
timing <- unique(x[timing_columns])
timing$direct_total_seconds <- timing$direct_F_seconds +
  timing$direct_derivative_prep_seconds + timing$direct_legacy_kernel_seconds
timing$lookup_reused_table_seconds <- timing$lookup_F_and_derivative_seconds +
  timing$lookup_candidate_kernel_seconds
timing$lookup_one_off_seconds <- timing$table_build_seconds + timing$lookup_reused_table_seconds
timing$reused_table_speedup <- timing$direct_total_seconds / timing$lookup_reused_table_seconds
utils::write.csv(timing, file.path(out, "timing_summary.csv"), row.names = FALSE)

historical <- utils::read.csv("benchmarks/profile_lookup_end_to_end_q10_n200_B4999.csv", check.names = FALSE)
new_historical <- subset(
  x, n == 200L & B == 4999L &
    ((comparison_path == "table_only" & quantity %in% c(
      "F_matrix", "dot_F_matrix", "information_matrix",
      "dot_F_times_information_inverse", "complete_fast_correction_contribution",
      "bootstrap_replicates")) |
     (comparison_path == "combined" & quantity %in% c(
       "observed_statistic", "critical_value", "p_value")))
)
joined <- merge(new_historical, historical,
  by = c("model", "quantity", "statistic"), suffixes = c("_new", "_archived"))
for (field in c("max_abs_difference", "rmse_difference", "mean_difference")) {
  joined[[paste0(field, "_delta")]] <-
    joined[[paste0(field, "_new")]] - joined[[paste0(field, "_archived")]]
}
utils::write.csv(joined, file.path(out, "historical_comparison.csv"), row.names = FALSE)

p_value_max <- max(subset(x, comparison_path == "combined" & quantity == "p_value")$max_abs_difference)
cat(sprintf("Kernel-only maximum difference: %.5g\n", kernel_max))
cat(sprintf("Combined maximum p-value difference: %.5g\n", p_value_max))
cat(sprintf("Wrote summaries to %s\n", out))
RS

print -- "Pilot completed successfully. Review: $out/precision_summary.csv and $out/timing_summary.csv"
