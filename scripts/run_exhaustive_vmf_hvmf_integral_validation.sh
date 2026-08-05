#!/bin/zsh

# Reproducible validation of the deterministic ("quadrature"/integral)
# distance-profile derivative for vMF and HvMF.  This runner does not modify
# production results, the Monte Carlo implementation, or any paper source.

set -u
set -o pipefail

cd "${0:A:h}/.." || exit 1

validation_stamp=$(date +%Y%m%d_%H%M%S)
validation_dir="benchmarks/vmf_hvmf_integral_validation_${validation_stamp}"
mkdir -p "$validation_dir"

export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

status_file="$validation_dir/status.txt"
unit_log="$validation_dir/unit_tests.log"
pilot_log="$validation_dir/console.log"

{
  echo "status: running"
  echo "started_at: $(date -Iseconds)"
  echo "git_head: $(git rev-parse HEAD)"
  echo "cores: 15"
  echo "production_results_modified: FALSE"
  echo "tex_modified: FALSE"
  echo "score_mc_modified: FALSE"
} > "$status_file"

echo "[1/2] Tests matematicos y numericos especificos de dot F"
if ! Rscript -e 'testthat::test_file("tests/testthat/test_deterministic_profile_derivatives.R", reporter = "progress")' \
    2>&1 | tee "$unit_log"; then
  {
    echo "status: failed_unit_tests"
    echo "finished_at: $(date -Iseconds)"
  } >> "$status_file"
  echo "ERROR: fallaron los tests especificos. Log: $unit_log" >&2
  exit 1
fi

echo "[2/2] Barrido emparejado AoS: integral, refinamiento y score_mc"
echo "      15 procesos; B=1999; limite total del piloto=15 minutos"

if ! Rscript scripts/run_aos_vmf_hvmf_integral_precision_pilot.R \
    --output-dir="$validation_dir/paired_pilot" \
    --main-M=60 \
    --refine-M=20 \
    --ultra-M=10 \
    --mc-M=10 \
    --mc-all-designs=true \
    --kernel-M=0 \
    --B=1999 \
    --cores=15 \
    --seed=20260803 \
    --cvm-block-size=25 \
    --max-wall-minutes=15 \
    --checkpoint-tasks=5 \
    2>&1 | tee "$pilot_log"; then
  {
    echo "status: failed_paired_pilot"
    echo "finished_at: $(date -Iseconds)"
  } >> "$status_file"
  echo "ERROR: fallo el piloto emparejado. Log: $pilot_log" >&2
  exit 1
fi

pilot_status="$validation_dir/paired_pilot/progress_status.txt"
if [[ ! -s "$pilot_status" ]] || \
    rg -q 'failed_mode_runs: [1-9]|nonconforming_mode_runs: [1-9]|wall_time_cap_reached: TRUE' "$pilot_status"; then
  {
    echo "status: failed_or_nonconforming_results"
    echo "finished_at: $(date -Iseconds)"
  } >> "$status_file"
  echo "ERROR: hay tareas fallidas o no conformes. Revisa: $pilot_status" >&2
  exit 1
fi

{
  echo "status: completed"
  echo "finished_at: $(date -Iseconds)"
  echo "results: $validation_dir"
} >> "$status_file"

echo
echo "VALIDACION COMPLETADA SIN ERRORES"
echo "Resultados: $validation_dir"
echo "Resumen de precision: $validation_dir/paired_pilot/paired_comparison_summary.csv"
echo "Estado: $status_file"
