# Gaussian Process Analysis - Code Structure

## 📁 Files

### Main Analysis (Canonical implementations)
- `gaussian_process_normal.R` — Normal distribution canonical implementation (vectorized & parallelized)
  - Functions: `cov_normal`, `row_cov_normal`, `simulate_limit_gaussian_normal`, `simulate_empirical_process_normal`.
- `gaussian_process_vmf.R` — vMF canonical implementation (vectorized & parallelized)
  - Functions: `cov_vmf`, `row_cov_vmf`, `simulate_limit_gaussian_vmf`, `simulate_empirical_process_vmf`.
  - `gaussian_process_analysis.R` was a legacy analysis script that duplicated functionality; it has been deprecated and removed in favor of the canonical files above.

### Performance Testing
- **`gaussian_process_tests.R`** (7.9KB) - TESTING/BENCHMARKING
  - Performance comparison between original and vectorized implementations
  - Main function: `test_vectorization_performance()`
  - Baseline implementations for comparison:
    - `covariance_gaussian_process()` - Scalar version (slow)
    - `create_covariance_matrix_loops()` - Nested loops version (slow)

## 🎯 Usage

### Run Main Analysis
Use the canonical functions in `convergence_empirical_process/gaussian_process_normal.R` and `convergence_empirical_process/gaussian_process_vmf.R`.
Example: source the canonical file and call the wrapper:
```r
source(file.path("convergence_empirical_process", "gaussian_process_normal.R"))
result <- simulate_limit_gaussian_normal(omega_grid, t_grid, mu, sigma)
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