#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
out="${1:-paper/rendered}"
mkdir -p "$out"/vmf "$out"/sunspots "$out"/wind "$out"/simplex

Rscript scripts/plot_paper_vmf_convergence_from_rds.R \
  --input_dir=paper/artifacts/vmf_convergence --output_dir="$out/vmf" \
  --h0=simple --adjust_50=5 --adjust_100=1 --adjust_500=0.85
Rscript scripts/plot_paper_vmf_convergence_from_rds.R \
  --input_dir=paper/artifacts/vmf_convergence --output_dir="$out/vmf" \
  --h0=composite --adjust_50=4.5 --adjust_100=1 --adjust_500=1
Rscript -e 'source("scripts/plot_spherical_parametric_hdr_bands.R"); run_sunspot_parametric_hdr_bands(source_dir="paper/artifacts/sunspots", output_dir=commandArgs(TRUE)[1])' "$out/sunspots"
Rscript real_data/wind/plot_risoe_extended_r_density_panels.R \
  --cases_dir=paper/artifacts/wind --output_dir="$out/wind"
for name in ArcticLake Sediments PogoJump WhiteCells_microscopic Yatquat_panel; do
  slug=$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')
  Rscript scripts/run_sediments_simplex_contours.R \
    --dataset="$name" --input_csv="paper/artifacts/simplex/$slug.csv" \
    --output_dir="$out/simplex"
done
