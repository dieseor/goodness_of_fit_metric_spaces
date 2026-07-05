#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

N_CORES="${N_CORES:-12}"
VALIDATION_B="${VALIDATION_B:-1000}"
VALIDATION_M_OUTER="${VALIDATION_M_OUTER:-5}"
DERIVATIVE_MC_SIZE="${DERIVATIVE_MC_SIZE:-1000}"
DERIVATIVE_MC_SEED="${DERIVATIVE_MC_SEED:-20260613}"
OUTPUT_ROOT="${OUTPUT_ROOT:-output}"
LOG_DIR="$OUTPUT_ROOT/logs"

mkdir -p "$LOG_DIR"

run_step() {
  local step_idx="$1"
  local n_steps="$2"
  local label="$3"
  local log_path="$4"
  shift 4

  printf '\n========================================================================\n'
  printf '[step %s/%s] %s\n' "$step_idx" "$n_steps" "$label"
  printf '[step %s/%s] log: %s\n' "$step_idx" "$n_steps" "$log_path"
  printf '========================================================================\n'

  "$@" 2>&1 | tee "$log_path"

  printf '[step %s/%s] completed: %s\n' "$step_idx" "$n_steps" "$label"
}

printf 'Starting overnight fast-multiplier pipeline at %s\n' "$(date)"
printf 'ROOT_DIR=%s\n' "$ROOT_DIR"
printf 'N_CORES=%s\n' "$N_CORES"
printf 'VALIDATION_B=%s\n' "$VALIDATION_B"
printf 'VALIDATION_M_OUTER=%s\n' "$VALIDATION_M_OUTER"
printf 'DERIVATIVE_MC_SIZE=%s\n' "$DERIVATIVE_MC_SIZE"
printf 'DERIVATIVE_MC_SEED=%s\n' "$DERIVATIVE_MC_SEED"
printf 'OUTPUT_ROOT=%s\n' "$OUTPUT_ROOT"

run_step 1 2 "Fast multiplier regression tests" \
  "$LOG_DIR/01_fast_multiplier_regression_tests.log" \
  env RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_AUTO_SNAPSHOT=FALSE \
  Rscript scripts/run_fast_multiplier_tests.R

run_step 2 2 "Sequential fast validation and paper experiments" \
  "$LOG_DIR/02_fast_validation_and_paper.log" \
  env RENV_CONFIG_SANDBOX_ENABLED=FALSE RENV_CONFIG_AUTO_SNAPSHOT=FALSE \
  Rscript scripts/run_all_fast_multiplier_sequential.R \
  --n_cores="$N_CORES" \
  --validation_B="$VALIDATION_B" \
  --validation_M_outer="$VALIDATION_M_OUTER" \
  --derivative_mc_size="$DERIVATIVE_MC_SIZE" \
  --derivative_mc_seed="$DERIVATIVE_MC_SEED" \
  --output_root="$OUTPUT_ROOT"

printf '\nOvernight fast-multiplier pipeline completed at %s\n' "$(date)"
