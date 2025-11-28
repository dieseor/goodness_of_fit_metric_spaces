# Gaussian Process Analysis - Code Structure

The goal is to show that the empirical process converges to the correct Gaussian process under the null hypotesis. This is tested for simple and composite hypotheses, for the normal and thw vMF distributions.

## Files

### Main files and functions
- `gaussian_process_normal.R` — Specific functions for the normal distribution (vectorized & parallelized)
  - Functions: `cov_normal()`, `row_cov_normal()`, `simulate_limit_gaussian_normal()`, `simulate_empirical_process_normal()`, `visualize_convergence_to_limit_vmf()`.
- `gaussian_process_vmf.R` — Specific functions for the vMF distribution (vectorized & parallelized)
  - Functions: `cov_vmf()`, `row_cov_vmf()`, `simulate_limit_gaussian_vmf()`, `simulate_empirical_process_vmf()`, `visualize_convergence_to_limit_vmf()`.

### Performance Testing
- **`gaussian_process_tests.R`** (7.9KB) - TESTING/BENCHMARKING
  - Performance comparison between original and vectorized implementations
  - Main function: `test_vectorization_performance()`
  - Baseline implementations for comparison, used to calculate the covariance matrix computing one element at a time: `covariance_gaussian_process()` and `create_covariance_matrix_loops()`.

## Usage

### Obtain results
- `run_gaussian_normal.R` - Produce results for the normal dsitribution
### Run Performance Tests
```r
source("gaussian_process_tests.R")
performance_results <- test_vectorization_performance()
```

## Configuration

Parameters in `gaussian_process_analysis.R` (lines 19-37):
- `OMEGA_MIN`, `OMEGA_MAX`, `OMEGA_POINTS` - Omega grid
- `T_MIN`, `T_MAX`, `T_POINTS` - Distance threshold grid
- `N_SAMPLE` - Sample size for empirical process
- `M_SIMULATIONS` - Number of Monte Carlo simulations
- `MU_TRUE`, `SIGMA_TRUE` - True distribution parameters
- `N_CORES` - Number of parallel cores

## Output

- Histogram plots: `output/gaussian_process/comparison_histogram_NxM_grid.png`
- KS test results printed to console
- Returns list with `limit_values`, `empirical_values`, `plot`, `ks_test`

## Function Naming Convention

- `compute_*` - Auxiliary/internal functions
- `simulate_*` - Main simulation functions
- `create_*` - Matrix/data structure creation
- `compare_*` - High-level analysis functions