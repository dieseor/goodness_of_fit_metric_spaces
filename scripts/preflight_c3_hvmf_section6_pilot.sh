#!/bin/bash
# Lightweight preflight for the HvMF Section 6 pilot.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash scripts/preflight_c3_hvmf_section6_pilot.sh {hvmf_1_mixture|hvmf_2_angular} [2,5]" >&2
  exit 2
fi
scenario="$1"
dimensions="${2:-2,5}"
case "$scenario" in hvmf_1_mixture|hvmf_2_angular) ;; *) echo "Unsupported scenario '$scenario'." >&2; exit 2 ;; esac
case "$dimensions" in 2|5|2,5) ;; *) echo "Unsupported dimensions '$dimensions'." >&2; exit 2 ;; esac
[[ -f scripts/run_section6_new_scenarios.R ]] || { echo "Run from repository root." >&2; exit 2; }
module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load R/4.2.1
export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"
Rscript --vanilla -e 'source("scripts/run_section6_new_scenarios.R", encoding="UTF-8"); d <- make_section6_design("hvmf", dimensions=as.integer(strsplit(commandArgs(TRUE)[1], ",")[[1]]), n_values=c(50L,100L,200L,400L), beta_values=c(0,0.5), scenarios=commandArgs(TRUE)[2]); print(d[,c("scenario","d","n","beta")]); stopifnot(all(d$beta %in% c(0,0.5)))' "$dimensions" "$scenario"
