#!/bin/bash

#SBATCH --job-name=rspike_prod
#SBATCH --account=uc3m-influentia-15
#SBATCH --output=logs/rspike_finalpilot_%j.out
#SBATCH --error=logs/rspike_finalpilot_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:40:00
#SBATCH --mem=7G
#SBATCH --partition=batch

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch script CONFIG BETA"
    exit 2
fi

mean_config="$1"
beta_value="$2"
beta_label=$(echo "$beta_value" | sed "s/\\./p/")

cd "$SLURM_SUBMIT_DIR"

output_dir="simulation_results/restricted_spiked_normal_finalpilot/${mean_config}_d2_d5_n50_100_200_400_M1000_B5000_beta${beta_label}_Nderiv10000"

mkdir -p logs

module purge
module load gnu12/12.2.0
module load gsl/2.7.1
module load openblas/0.3.21
module load R/4.2.1

export R_LIBS_USER="$HOME/R/c3-R-4.2-library"
export RENV_CONFIG_AUTOLOADER_ENABLED="FALSE"

export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "===== Restricted-spiked final pilot ====="
echo "output_dir=${output_dir}"
echo "hostname=$(hostname)"
echo "date=$(date --iso-8601=seconds)"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"

Rscript --vanilla scripts/check_c3_restricted_spiked_normal.R \
  --mode=preflight \
  --mean_config="${mean_config}" \
  --output_dir="${output_dir}" \
  --M=1000 \
  --B=5000 \
  --dimensions=2,5 \
  --n_values=50,100,200,400 \
  --beta_values=${beta_value} \
  --lambda=2 \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --allow_new=true

start=$(date +%s)

Rscript --vanilla scripts/run_restricted_spiked_normal_covariance_alternatives.R \
  --mean_config="${mean_config}" \
  --output_dir="${output_dir}" \
  --M=1000 \
  --B=5000 \
  --dimensions=2,5 \
  --n_values=50,100,200,400 \
  --beta_values=${beta_value} \
  --lambda=2 \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --cores="${SLURM_CPUS_PER_TASK}" \
  --checkpoint_results=2000 \
  --show_progress=true

end=$(date +%s)
wall=$((end-start))

Rscript --vanilla scripts/check_c3_restricted_spiked_normal.R \
  --mode=final \
  --mean_config="${mean_config}" \
  --output_dir="${output_dir}" \
  --M=1000 \
  --B=5000 \
  --dimensions=2,5 \
  --n_values=50,100,200,400 \
  --beta_values=${beta_value} \
  --lambda=2 \
  --derivative_mc_size=10000 \
  --cvm_block_size=50 \
  --runner_wall_seconds="${wall}" \
  --allocated_cores="${SLURM_CPUS_PER_TASK}"

echo "===== END ====="
echo "wall_seconds=${wall}"
echo "output_dir=${output_dir}"
