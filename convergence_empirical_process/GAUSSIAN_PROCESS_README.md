# Gaussian Process Analysis - Code Structure

## 📁 Files

### Main Analysis
- **`gaussian_process_analysis.R`** (16KB) - PRODUCTION CODE
  - Clean, optimized implementation for generating comparison plots
  - Uses vectorized + parallelized computation
  - Main function: `compare_processes()`
  - Functions included:
    - `compute_distance_profile()` - Theoretical distance profile
    - `compute_joint_probability()` - Joint probability for covariance
    - `compute_covariance_row()` - Vectorized row computation (internal)
    - `create_covariance_matrix()` - Full matrix creation (vectorized + parallel)
    - `simulate_limit_gaussian()` - Simulate Gaussian process limit
    - `simulate_empirical_process()` - Simulate empirical process
    - `compare_processes()` - Main comparison function

### Performance Testing
- **`gaussian_process_tests.R`** (7.9KB) - TESTING/BENCHMARKING
  - Performance comparison between original and vectorized implementations
  - Main function: `test_vectorization_performance()`
  - Baseline implementations for comparison:
    - `covariance_gaussian_process()` - Scalar version (slow)
    - `create_covariance_matrix_loops()` - Nested loops version (slow)

## 🎯 Usage

### Run Main Analysis
```r
source("gaussian_process_analysis.R")
# Adjust parameters in the CONFIGURATION section (lines 15-35)
# Then run the script or call:
result <- compare_processes()
```

### Run Performance Tests
```r
source("gaussian_process_tests.R")
performance_results <- test_vectorization_performance()
```

## 🔧 Configuration

Edit parameters in `gaussian_process_analysis.R` (lines 19-37):
- `OMEGA_MIN`, `OMEGA_MAX`, `OMEGA_POINTS` - Omega grid
- `T_MIN`, `T_MAX`, `T_POINTS` - Distance threshold grid
- `N_SAMPLE` - Sample size for empirical process
- `M_SIMULATIONS` - Number of Monte Carlo simulations
- `MU_TRUE`, `SIGMA_TRUE` - True distribution parameters
- `N_CORES` - Number of parallel cores

## 📊 Output

- Histogram plots: `output/gaussian_process/comparison_histogram_NxM_grid.png`
- KS test results printed to console
- Returns list with `limit_values`, `empirical_values`, `plot`, `ks_test`

### Plot smoothing and missing limit curves
 - `density_adjust` (numeric) — parameter used by `visualize_convergence_to_limit_*` functions to scale kernel density bandwidth; increase (e.g. 1.25, 1.5, or 2.0) to reduce 'wiggly' density curves when using large `M_SIMULATIONS`.
 - If the covariance matrix is not numerically positive semi-definite, `simulate_limit_gaussian()` will error: the function does **not** apply any PSD correction and requires the user to provide a corrected covariance matrix upstream.
The PSD diagnostic CSVs contain the following fields of interest:
- `*_n_neg` – number of eigenvalues < -1e-10 (default, strict)
- `*_n_neg_1e8` – number of eigenvalues < -1e-8 (looser)
- `*_n_neg_1e6` – number of eigenvalues < -1e-6 (loose)

This gives you the ability to see how many eigenvalues are slightly negative at different numerical thresholds.
If you want to intentionally produce a corrected limit curve for visualization, compute a PSD-corrected covariance matrix externally using your preferred method and pass the corrected matrix to `simulate_limit_gaussian()`.

## 🚀 Performance Notes

- Small grids (5x5): Parallelization overhead may dominate
- Large grids (50x50+): Vectorization provides exponential speedup
- See `gaussian_process_tests.R` for detailed benchmarks

## 📝 Function Naming Convention

- `compute_*` - Auxiliary/internal functions
- `simulate_*` - Main simulation functions
- `create_*` - Matrix/data structure creation
- `compare_*` - High-level analysis functions