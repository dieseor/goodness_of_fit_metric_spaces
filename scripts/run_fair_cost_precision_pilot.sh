#!/usr/bin/env zsh

# Fair production-path timing and strict table-precision pilot.  This script
# is validation-only and never changes a production method selector.

setopt pipefail

repo_root="${0:A:h:h}"
cd "$repo_root" || exit 1
cores="${FAIR_PILOT_CORES:-10}"
if [[ ! "$cores" =~ '^[1-9][0-9]*$' ]]; then
  print -u2 -- "FAIR_PILOT_CORES must be a positive integer."
  exit 2
fi

out="${1:-benchmarks/fair_cost_precision_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$out/kernel" "$out/profile_precision" "$out/status" || exit 1
out="${out:A}"
status_dir="$out/status"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

tasks=()
for model in vmf hvmf; do
  for q in 2 10; do
    for n in 50 100 200 400; do
      tasks+=("kernel $model $q $n")
    done
  done
done
for model in vmf hvmf; do
  for q in 2 10; do
    if [[ "$model" == vmf ]]; then
      tasks+=("profile $model $q 1.5 2.5 1 3.141592653589793 kappa2")
      tasks+=("profile $model $q 7 14 1 3.141592653589793 kappa10")
      tasks+=("profile $model $q 40 60 1 3.141592653589793 kappa50")
    else
      tasks+=("profile $model $q 1.5 2.5 3 4 kappa2")
      tasks+=("profile $model $q 5 12 3 4 kappa10")
      tasks+=("profile $model $q 40 60 3 4 kappa50")
    fi
  done
done

total=${#tasks}
completed=0
if [[ "${FAIR_PILOT_RESUME:-0}" == "1" ]]; then
  pending_tasks=()
  for task in "${tasks[@]}"; do
    kind="${task%% *}"
    if [[ "$kind" == kernel ]]; then
      read -r ignored model q n <<< "$task"
      stem="kernel_${model}_q${q}_n${n}"
      result="$out/kernel/$stem/summary.csv"
    else
      read -r ignored model q kmin kmax gmax tmax label <<< "$task"
      stem="profile_${model}_q${q}_${label}"
      result="$out/profile_precision/$stem.csv"
    fi
    if [[ -f "$result" ]]; then
      (( completed++ ))
    else
      pending_tasks+=("$task")
    fi
  done
  tasks=("${pending_tasks[@]}")
fi
pending_total=${#tasks}
next=1 active=0 failed=0

render_progress() {
  local width=36 filled empty left right percent
  filled=$(( completed * width / total ))
  empty=$(( width - filled ))
  left=$(printf '%*s' "$filled" '' | tr ' ' '#')
  right=$(printf '%*s' "$empty" '' | tr ' ' '-')
  percent=$(( completed * 100 / total ))
  printf '\r[%s%s] %3d%%  %2d/%2d complete  %2d active  %d failed' \
    "$left" "$right" "$percent" "$completed" "$total" "$active" "$failed"
}

launch_task() {
  local task="$1" kind model q n kmin kmax gmax tmax label stem log log_dir marker
  kind="${task%% *}"
  if [[ "$kind" == kernel ]]; then
    read -r kind model q n <<< "$task"
    stem="kernel_${model}_q${q}_n${n}"
    log_dir="$out/kernel"
  else
    read -r kind model q kmin kmax gmax tmax label <<< "$task"
    stem="profile_${model}_q${q}_${label}"
    log_dir="$out/profile_precision"
  fi
  log="$log_dir/${stem}.log"
  marker="$status_dir/${stem}"

  (
    if [[ "$kind" == kernel ]]; then
      Rscript --vanilla scripts/validate_fast_multiplier_cpp_kernel_end_to_end.R \
        --models="$model" --dimensions="$q" --sizes="$n" --B=4999 --cores=1 \
        --derivative-method=quadrature --output-dir="$out/kernel/$stem" >"$log" 2>&1
    else
      Rscript --vanilla scripts/benchmark_profile_lookup_interpolation.R \
        --model="$model" --q="$q" --kappa-min="$kmin" --kappa-max="$kmax" \
        --geometry-max="$gmax" --t-max="$tmax" \
        --n-kappa=41 --n-geometry=65 --n-t=513 --queries=200 --stencil-size=6 \
        --integration-grid-size=16385 --reference-grid-size=32769 --cores=1 \
        --output="$out/profile_precision/${stem}.csv" >"$log" 2>&1
    fi
    if (( $? == 0 )); then
      print -r -- "$stem" > "${marker}.done"
    else
      print -r -- "$stem" > "${marker}.failed"
    fi
  ) &
}

while (( next <= pending_total && active < cores )); do
  launch_task "${tasks[$next]}"
  (( next++ ))
  (( active++ ))
done

print -- "Fair timing + precision pilot"
print -- "Output: $out"
print -- "Workers: $cores; 16 production-kernel and 12 strict-precision tasks."
if [[ "${FAIR_PILOT_RESUME:-0}" == "1" ]]; then
  print -- "Resume mode: $completed completed task(s) were retained."
fi

while (( completed < total )); do
  for marker in "$status_dir"/*.done(N) "$status_dir"/*.failed(N); do
    [[ -f "$marker" && ! -f "${marker}.reported" ]] || continue
    touch "${marker}.reported"
    (( completed++ ))
    (( active-- ))
    if [[ "$marker" == *.failed ]]; then
      (( failed++ ))
      stem="${${marker:t}%.failed}"
      log_kernel="$out/kernel/${stem}.log"
      log_profile="$out/profile_precision/${stem}.log"
      log="$log_kernel"
      [[ -f "$log" ]] || log="$log_profile"
      print -u2 -- "\nERROR: ${stem} failed. Last 25 log lines:"
      tail -n 25 "$log" >&2
    fi
  done
  while (( next <= pending_total && active < cores )); do
    launch_task "${tasks[$next]}"
    (( next++ ))
    (( active++ ))
  done
  render_progress
  sleep 1
done
print

if (( failed > 0 )); then
  print -u2 -- "Pilot finished with $failed failed task(s). Inspect $out/kernel/*.log and $out/profile_precision/*.log."
  exit 1
fi

Rscript --vanilla - "$out" <<'RS'
out <- commandArgs(trailingOnly = TRUE)[1L]
kernel_files <- list.files(file.path(out, "kernel"), pattern = "summary\\.csv$", recursive = TRUE, full.names = TRUE)
profile_files <- list.files(file.path(out, "profile_precision"), pattern = "\\.csv$", full.names = TRUE)
if (length(kernel_files) != 16L || length(profile_files) != 12L) stop("Missing completed result files.")

kernel <- do.call(rbind, lapply(kernel_files, utils::read.csv, check.names = FALSE))
kernel$legacy_preparation_seconds <- kernel$legacy_common_observed_seconds + kernel$legacy_fast_prep_seconds
kernel$candidate_preparation_seconds <- kernel$candidate_common_observed_seconds + kernel$candidate_fast_prep_seconds
kernel$loop_speedup <- kernel$legacy_fast_loop_seconds / kernel$candidate_fast_loop_seconds
kernel$total_speedup <- kernel$legacy_fast_total_seconds / kernel$candidate_fast_total_seconds
utils::write.csv(kernel, file.path(out, "production_kernel_timing_decomposition.csv"), row.names = FALSE)

profiles <- do.call(rbind, lapply(profile_files, utils::read.csv, check.names = FALSE))
utils::write.csv(profiles, file.path(out, "profile_table_precision_summary.csv"), row.names = FALSE)

cat(sprintf("Kernel-only maximum relative discrepancy: %.5g\n", max(c(kernel$max_rel_ks_difference, kernel$max_rel_cvm_difference))))
cat(sprintf("Maximum table F error: %.5g\n", max(profiles$max_abs_F_error)))
cat(sprintf("Maximum table derivative-component error: %.5g\n", max(profiles$max_component_derivative_error)))
RS

print -- "Pilot completed successfully. Results: $out"
