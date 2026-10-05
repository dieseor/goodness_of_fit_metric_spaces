#!/usr/bin/env bash
# Monte Carlo campaigns behind the rejection-rate table of the paper, with the
# settings and seeds of the runs used in the paper (one line per campaign).
# Run all of them, or only some scenarios:
#
#   bash paper/run_simulations.sh          # scenarios 1-8
#   bash paper/run_simulations.sh 5 6      # only scenarios 5 and 6
#
# CORES sets the cores per campaign (default 8); results do not depend on it,
# because every replication has its own seed. The runners save checkpoints
# and resume an interrupted campaign. These runs were made on a cluster and
# take a long time. Rscript is called without --vanilla so that the renv
# library (the package versions in renv.lock) is used.
set -euo pipefail

cd "$(dirname "$0")/.."
CORES="${CORES:-8}"
SCENARIOS="${*:-1 2 3 4 5 6 7 8}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

campaign() {
  local scenario="$1" output_dir="$2" script="$3"
  shift 3
  case " $SCENARIOS " in *" $scenario "*) ;; *) return 0 ;; esac
  printf '\n===== scenario %s: %s =====\n' "$scenario" "$output_dir"
  Rscript "$script" "$@" --cores="$CORES" --output_dir="$output_dir"
}


# Scenario 1: restricted spiked normal.
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n50_100_200_400_M1000_B1000_beta0_Nderiv10000 \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n50_100_200_400_M1000_B1000_beta0p5_Nderiv10000 \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n50_100_200_400_M1000_B1000_beta1_Nderiv10000 \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n800_M1000_B1000_beta0_Nderiv10000_cachetrue \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n800_M1000_B1000_beta0p5_Nderiv10000_cachetrue \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833
campaign 1 simulation_results/restricted_spiked_normal_B1000/diagonal_100_d2_d5_n800_M1000_B1000_beta1_Nderiv10000_cachetrue \
  scripts/run_restricted_spiked_normal_covariance_alternatives.R \
    --mean_config=diagonal_100 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --lambda=2 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260833

# Scenario 2: normal versus multivariate t.
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B1000_beta0_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B1000_beta0p5_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=50,100,200,400 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n50_100_200_400_M1000_B1000_beta1_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=50,100,200,400 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n800_M1000_B1000_beta0_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n800_M1000_B1000_beta0p5_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=800 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d2_nu3_n800_M1000_B1000_beta1_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=2 \
    --nu=3 \
    --n_values=800 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B1000_beta0_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B1000_beta0p5_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=50,100,200,400 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n50_100_200_400_M1000_B1000_beta1_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=50,100,200,400 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n800_M1000_B1000_beta0_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n800_M1000_B1000_beta0p5_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=800 \
    --beta_values=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728
campaign 2 simulation_results/section6_new_scenarios/final_normal_sigma_Id_t_d5_nu6_n800_M1000_B1000_beta1_Nderiv10000 \
  scripts/run_normal_sigma_Id_t_pilot.R \
    --dimensions=5 \
    --nu=6 \
    --n_values=800 \
    --beta_values=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --cache_corrections=auto \
    --seed=20260728

# Scenario 3: logistic Gaussian versus Dirichlet(1.5).
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta0_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta0p5_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta1_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta0_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=800 \
    --beta=0 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta0p5_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=800 \
    --beta=0.5 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904
campaign 3 simulation_results/lg_dirichlet15_mu_only_final_beta1_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_dirichlet15_mu_only.R \
    --scenario=dirichlet15 \
    --d=2,5 \
    --n=800 \
    --beta=1 \
    --M=1000 \
    --B=1000 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260904

# Scenario 4: logistic Gaussian AR(1) versus t4.
campaign 4 simulation_results/lg_t4_ar1_final_beta0_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=0 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903
campaign 4 simulation_results/lg_t4_ar1_final_beta0p5_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=0.5 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903
campaign 4 simulation_results/lg_t4_ar1_final_beta1_d2_d5_n50_100_200_400_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=50,100,200,400 \
    --beta=1 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903
campaign 4 simulation_results/lg_t4_ar1_final_beta0_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=800 \
    --beta=0 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903
campaign 4 simulation_results/lg_t4_ar1_final_beta0p5_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=800 \
    --beta=0.5 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903
campaign 4 simulation_results/lg_t4_ar1_final_beta1_d2_d5_n800_M1000_B1000_Nderiv10000 \
  scripts/run_logistic_gaussian_sigma_shape_scenarios.R \
    --scenario=t4 \
    --d=2,5 \
    --n=800 \
    --beta=1 \
    --M=1000 \
    --B=1000 \
    --nu=4 \
    --t_standardized=true \
    --t_ar1_rho=0.5 \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260903

# Scenario 5: vMF with fixed kappa, orthogonal mixture.
campaign 5 simulation_results/section6_new_scenarios/final_calibration_vmf_31_orthogonal_mu_only_kappa2_beta0_d2_5_n50_100_200_400_M1000_B1000 \
  scripts/run_vmf_mu_only_fixed_kappa_pilot.R \
    --kappa_values=2 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260908
campaign 5 simulation_results/section6_new_scenarios/final_power_vmf_31_orthogonal_mu_only_kappa2_beta0p5_1_d2_5_n50_100_200_400_M1000_B1000 \
  scripts/run_vmf_mu_only_fixed_kappa_pilot.R \
    --kappa_values=2 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260915
campaign 5 simulation_results/section6_new_scenarios/final_calibration_vmf_31_orthogonal_mu_only_kappa2_beta0_d2_5_n800_M1000_B1000 \
  scripts/run_vmf_mu_only_fixed_kappa_pilot.R \
    --kappa_values=2 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260908
campaign 5 simulation_results/section6_new_scenarios/final_power_vmf_31_orthogonal_mu_only_kappa2_beta0p5_1_d2_5_n800_M1000_B1000 \
  scripts/run_vmf_mu_only_fixed_kappa_pilot.R \
    --kappa_values=2 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260916

# Scenario 6: vMF versus projected normal.
campaign 6 simulation_results/section6_new_scenarios/final_calibration_vmf_32_mean_d_beta0_d2_5_n50_100_200_400_M1000_B1000 \
  scripts/run_vmf_antipodal_fixed_kappa_pilot.R \
    --scenario_type=projected_normal_mean_d \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260902
campaign 6 simulation_results/section6_new_scenarios/final_power_vmf_32_projected_normal_2sqrt_d_kappa2d_beta_half_d2_5_n50_100_200_400_M1000_B1000 \
  scripts/run_vmf_antipodal_fixed_kappa_pilot.R \
    --scenario_type=projected_normal_2sqrt_d_kappa_2d_beta_half \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260913
campaign 6 simulation_results/section6_new_scenarios/final_calibration_vmf_32_mean_d_beta0_d2_5_n800_M1000_B1000 \
  scripts/run_vmf_antipodal_fixed_kappa_pilot.R \
    --scenario_type=projected_normal_mean_d \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260902
campaign 6 simulation_results/section6_new_scenarios/final_power_vmf_32_projected_normal_2sqrt_d_kappa2d_beta_half_d2_5_n800_M1000_B1000 \
  scripts/run_vmf_antipodal_fixed_kappa_pilot.R \
    --scenario_type=projected_normal_2sqrt_d_kappa_2d_beta_half \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --cvm_block_size=50 \
    --seed=20260914

# Scenario 7: HvMF radial mixture.
campaign 7 simulation_results/section6_new_scenarios/final_calibration_hvmf_1_mixture_beta0_d2_5_n50_100_200_400_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_1_mixture \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 7 simulation_results/section6_new_scenarios/final_hvmf_1_radial_c_sqrt2_beta0p5_1_d2_5_n50_100_200_400_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_1_radial_c_sqrt2 \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 7 simulation_results/section6_new_scenarios/final_calibration_hvmf_1_mixture_beta0_d2_5_n800_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_1_mixture \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 7 simulation_results/section6_new_scenarios/final_hvmf_1_radial_c_sqrt2_beta0p5_1_d2_5_n800_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_1_radial_c_sqrt2 \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830

# Scenario 8: HvMF angular concentration.
campaign 8 simulation_results/section6_new_scenarios/final_calibration_hvmf_2_angular_sqrt_d_concentration_beta0_d2_5_n50_100_200_400_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_2_angular_sqrt_d_concentration \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 8 simulation_results/section6_new_scenarios/final_power_hvmf_2_angular_sqrt_d_concentration_beta0p5_1_d2_5_n50_100_200_400_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_2_angular_sqrt_d_concentration \
    --dimensions=2,5 \
    --n_values=50,100,200,400 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 8 simulation_results/section6_new_scenarios/final_calibration_hvmf_2_angular_sqrt_d_concentration_beta0_d2_5_n800_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_2_angular_sqrt_d_concentration \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830
campaign 8 simulation_results/section6_new_scenarios/final_power_hvmf_2_angular_sqrt_d_concentration_beta0p5_1_d2_5_n800_M1000_B1000_quadrature_F_integral_pmax0999 \
  scripts/run_section6_new_scenarios.R \
    --family=hvmf \
    --scenarios=hvmf_2_angular_sqrt_d_concentration \
    --dimensions=2,5 \
    --n_values=800 \
    --beta_values=0.5,1 \
    --M=1000 \
    --B=1000 \
    --derivative_method=quadrature \
    --derivative_mc_size=10000 \
    --cvm_block_size=50 \
    --seed=20260830

printf '\nDone.\n'
