# Goodness-of-fit in metric spaces

This project explores goodness of fit tests and methods in metric space contexts.

The multiplier-bootstrap pipeline covers Euclidean normal and multivariate normal models, logistic-Gaussian models on the simplex, and directional models including vMF, HvMF, Jones--Pewsey, spherical cardioid, spherical Cauchy / Poisson-kernel, small-circle, and Watson-type models.

## How a goodness-of-fit test is run

The family-specific wrappers have a common execution path:

```
data and null hypothesis
        │
        ▼
multiplier_bootstrap_<family>()
bootstrap/multiplier_bootstrap.R       # public wrappers and common engine
        │
        ▼
make_<family>_spec()
bootstrap/model_specs.R                # common spec interface and adapter loader
bootstrap/*_model_spec.R               # one file per model family
        │
        ├── fit the null parameter and compute distances
        ├── evaluate theoretical distance profiles
        │   utils.R                     # shared profile formulae and numerical helpers
        │   bootstrap/*_model_spec.R    # family-specific profile and fitting code
        │
        ▼
multiplier_bootstrap_gof()
bootstrap/multiplier_bootstrap.R       # KS/CvM statistics and multiplier resampling
        │
        └── optional compiled kernels
            distance_profile_backend.R # opt-in backend selection and loading
            cpp/distance_profile_backend.cpp
```

`utils.R` is therefore the main location for generic mathematical and numerical utilities, including most shared theoretical distance-profile evaluations.  `bootstrap/model_specs.R` defines the model-specification interface and loads all family adapters. Each `bootstrap/*_model_spec.R` file owns the normalization, fitting, distance, profile, and fast-bootstrap components of one model family.

### Restricted mean-aligned spiked normal model

`bootstrap/restricted_spiked_normal_model_spec.R` implements the separate
restricted family
\(N_d(\theta, I_d + \lambda u(\theta)u(\theta)^\top)\), where
\(u(\theta)=\theta/\|\theta\|\), \(\theta\ne0\), and \(\lambda>0\).
It is not PPCA and does not estimate the mean and spike direction separately.
Its production MLE profiles over \(\lambda\) and obtains the leading
eigenvector at each value; `bootstrap/restricted_spiked_normal_openmx.R`
contains an independent OpenMx reference used only by the validation suite.
The custom estimator stops with an explicit error if its profiled radius is
numerically indistinguishable from zero; it does not regularise that case.

The paper's double-spike power experiment is implemented by
`scripts/run_restricted_spiked_normal_covariance_alternatives.R`. A single
runner selects one of four mean configurations with `--mean_config`:
`axis_075`, `axis_100`, `diagonal_150`, or `diagonal_100`. Its paper defaults
are `d=2,5`, `n=50,100,200,400`, `beta=0,0.25,0.5,1`, `M=1000`, `B=5000`,
`N_deriv=10000`, and `lambda=2`. Both KS and CvM use the observed sample
points and their unique observed distance thresholds, sharing the same fast
sample cache. The historical base seeds are retained per mean configuration
(`20260831` for the non-unit means and `20260833` for the unit means). Each
configuration has a separate output directory and a strict
manifest, atomic checkpoints, deterministic seeds, and an output lock. The C3
intermediate benchmark is submitted through
`scripts/run_c3_restricted_spiked_normal_benchmark.sbatch`; production Slurm
resources are intentionally chosen only after those benchmark measurements.

### Fast multivariate-normal loop

For `multiplier_bootstrap_mvnormal(..., bootstrap_method = "fast_multiplier")`,
the sample-based KS/CvM loop uses C++ by default.  Its independent controls are
`fast_multiplier_backend = c("cpp", "r")`, `fuse_ks_cvm = TRUE`, and
`cache_block_corrections = c("auto", "true", "false")`.  Both backends fuse KS
and CvM by default.  In `"auto"` mode, correction matrices are cached only when
the sample size is at most 500 and the estimated cache is at most 128 MiB;
the two explicit values override this policy.

The compiled loop is loaded lazily and is separate from
`distance_profile_backend`, whose default remains R.  The result diagnostics
record the requested and effective loop backend, fusion, and cache choices.

## Repository structure

```
goodness_of_fit_metric_spaces/
├── bootstrap/                # Test engine, generic model interface, and model adapters
│   ├── multiplier_bootstrap.R # Main bootstrap engine and public test wrappers
│   ├── model_specs.R          # Shared model-spec interface and adapter loader
│   ├── calibration_study.R    # Simulation/calibration study machinery
│   └── *_model_spec.R         # One family-specific adapter per file
├── utils.R                   # Shared distance profiles, distributions, and numerical tools
├── distance_profile_backend.R # Optional R/C++ backend dispatcher
├── cpp/                      # C++ kernels used only when the compiled backend is selected
├── distance_profiles/        # Distance-profile analyses and diagnostics
├── convergence_empirical_process/ # Empirical-process convergence experiments
├── bahadur/                  # Bahadur-efficiency analyses
├── real_data/                # Data-specific pipelines: comets, wind, sunspots, etc.
├── simulation_results/       # Outputs from simulation and power studies
├── output/
│   └── ...                   # General diagnostics, calibration, and validation outputs
├── scripts/                  # Reproducible runners, benchmarks, and plotting entry points
├── tests/
│   ├── testthat/             # Automated regression and integration tests
│   └── *.R                   # Focused numerical checks and experiments
└── benchmarks/               # Benchmark reports, including the C++ backend assessment
```

Generated results are kept either in the general `output/` tree, in `simulation_results/`, or next to the relevant real-data pipeline when they are specific to a dataset.

## Dependencies

The project automatically installs these packages if needed:
- ggplot2 (plotting)
- dplyr (data manipulation)
- readr (data reading)
- here (path management)
- ncdf4 (NetCDF reading for `wind/preprocess_risoe_modern_hvmf.R`)

## Reproducible development environment (renv)

This project uses `renv` to capture package dependencies so others can reproduce the environment.

To restore the project environment locally:

```bash
# From project root
# Install renv if needed and restore the environment
Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"
Rscript -e "remotes::install_github('rstudio/renv')"
Rscript -e "renv::restore()"
```
