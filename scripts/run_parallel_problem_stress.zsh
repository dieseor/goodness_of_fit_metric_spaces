#!/bin/zsh

# Interactive concurrency stress test for the numerical/profile paths that
# previously produced failures or expected compatibility warnings.  Every
# worker has a separate log; failures are printed immediately to this terminal.

setopt pipefail

cores=10
waves=7
output=""

for argument in "$@"; do
  case "$argument" in
    --cores=*) cores="${argument#--cores=}" ;;
    --waves=*) waves="${argument#--waves=}" ;;
    --output=*) output="${argument#--output=}" ;;
    --worker) worker=1 ;;
    *) ;;
  esac
done

if [[ -n "${worker:-}" ]]; then
  output=$2
  wave=$3
  label=$4
  filter=$5
  log="$output/wave${wave}_${label}.log"
  status_file="$output/wave${wave}_${label}.status"
  XDG_CACHE_HOME="$output/cache" \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    Rscript --vanilla -e \
      "testthat::test_dir(\"tests/testthat\", filter = \"${filter}\", reporter = \"summary\", stop_on_failure = TRUE, stop_on_warning = FALSE)" \
      > "$log" 2>&1
  exit_code=$?
  print -r -- "$exit_code" > "$status_file"
  if (( exit_code != 0 )); then
    print -u2 -r -- "\n[ERROR] wave ${wave}, ${label} (exit ${exit_code})"
    cat "$log" >&2
    print -u2 -r -- "[ERROR] log retained at ${log}\n"
  fi
  exit "$exit_code"
fi

if [[ ! "$cores" =~ '^[0-9]+$' ]] || (( cores < 1 )); then
  print -u2 -- "--cores must be a positive integer."
  exit 2
fi
if [[ ! "$waves" =~ '^[0-9]+$' ]] || (( waves < 1 )); then
  print -u2 -- "--waves must be a positive integer."
  exit 2
fi

if [[ -z "$output" ]]; then
  output="tests/stress_problem_cases_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$output/cache" || exit 1
git status --short > "$output/git_status_before.txt"
Rscript --vanilla -e 'sessionInfo()' > "$output/sessionInfo.txt" 2>&1

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MKL_NUM_THREADS=1

tasks=(
  "gaussian_1:gaussian_quadrature_derivatives"
  "gaussian_2:gaussian_quadrature_derivatives"
  "gaussian_3:gaussian_quadrature_derivatives"
  "c3_1:c3_section6_production_workflow"
  "c3_2:c3_section6_production_workflow"
  "multiplier:multiplier_bootstrap"
  "mvnormal:mvnormal_bootstrap"
  "deterministic:deterministic_profile_derivatives"
  "lookup:profile_lookup_interpolation"
  "comets:comets_spherical_cauchy_runner"
)
total=$(( waves * ${#tasks} ))
started=$(date +%s)

render_progress() {
  completed=$(find "$output" -maxdepth 1 -name '*.status' | wc -l | tr -d ' ')
  failed=0
  for status_file in "$output"/*.status(N); do
    [[ "$(cat "$status_file")" == 0 ]] || (( failed += 1 ))
  done
  filled=$(( completed * 32 / total ))
  empty=$(( 32 - filled ))
  bar="$(printf '%*s' "$filled" '' | tr ' ' '#')$(printf '%*s' "$empty" '' | tr ' ' '-')"
  elapsed=$(( $(date +%s) - started ))
  printf '\r[%s] %3d/%3d  failures=%d  elapsed=%02d:%02d' \
    "$bar" "$completed" "$total" "$failed" \
    $(( elapsed / 60 )) $(( elapsed % 60 ))
}

print -- "Results: $output"
print -- "Running $total test jobs in $waves wave(s), with at most $cores concurrent workers."

for wave in $(seq 1 "$waves"); do
  task_file="$output/wave${wave}_tasks.nul"
  : > "$task_file"
  for task in "${tasks[@]}"; do
    label="${task%%:*}"
    filter="${task#*:}"
    printf '%s\0%s\0%s\0%s\0' "$output" "$wave" "$label" "$filter" >> "$task_file"
  done
  wave_done="$output/wave${wave}.done"
  (
    xargs -0 -n 4 -P "$cores" zsh "$0" --worker < "$task_file"
    print -r -- "$?" > "$wave_done"
  ) &
  xargs_pid=$!
  while [[ ! -f "$wave_done" ]]; do
    render_progress
    sleep 1
  done
  wait "$xargs_pid"
  render_progress
done

print
print -- "\nCompleted. Non-zero statuses:"
for status_file in "$output"/*.status(N); do
  [[ "$(cat "$status_file")" == 0 ]] || print -- "${status_file}: $(cat "$status_file")"
done
print -- "Logs and statuses: $output"
