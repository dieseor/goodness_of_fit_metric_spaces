library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "multiplier_bootstrap.R"))

make_cpp_kernel_stream_prep <- function(B, n, p, n_centers, seed,
                                        extreme_weights = FALSE) {
  set.seed(seed)
  centered_weights <- matrix(stats::rexp(B * n) - 1, nrow = B, ncol = n)
  if (extreme_weights) {
    centered_weights[1L, ] <- seq(-1e8, 1e8, length.out = n)
    centered_weights[2L, ] <- rep(c(1e-8, -1e-8), length.out = n)
  }
  S_obs <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  order <- t(replicate(n_centers, sample.int(n)))
  tie_end <- matrix(0L, nrow = n_centers, ncol = n)
  tie_groups <- c(rep.int(3L, n %/% 3L), n %% 3L)
  tie_groups <- tie_groups[tie_groups > 0L]
  tie_end_row <- rep.int(cumsum(tie_groups), tie_groups)
  for (center in seq_len(n_centers)) {
    tie_end[center, ] <- tie_end_row
  }
  list(
    centered_weights = centered_weights,
    extreme_weights = extreme_weights,
    stream_prep = list(
      S_obs = S_obs,
      obs_order_matrix = order,
      tie_end_matrix = tie_end,
      correction_cache = list(
        values = matrix(
          stats::rnorm(n_centers * n * p),
          nrow = n_centers * n,
          ncol = p
        )
      ),
      n = n
    )
  )
}

test_that("the cache-friendly C++ kernel is the explicit default", {
  expect_identical(normalize_fast_multiplier_cpp_kernel(), "contiguous_double")
  expect_identical(
    normalize_fast_multiplier_cpp_kernel("contiguous_double"),
    "contiguous_double"
  )
  expect_error(
    normalize_fast_multiplier_cpp_kernel("structured"),
    "fast_multiplier_cpp_kernel"
  )
})

test_that("contiguous-double C++ kernel matches legacy with ties and extremes", {
  ensure_distance_profile_cpp_loaded()
  cases <- list(
    make_cpp_kernel_stream_prep(
      B = 13L, n = 17L, p = 4L, n_centers = 17L,
      seed = 20260811L, extreme_weights = FALSE
    ),
    make_cpp_kernel_stream_prep(
      B = 13L, n = 17L, p = 4L, n_centers = 17L,
      seed = 20260812L, extreme_weights = TRUE
    )
  )

  for (case in cases) {
    for (requested in list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))) {
      legacy <- compute_fast_sample_ks_cvm_stats_cpp(
        centered_weight_block = case$centered_weights,
        stream_prep = case$stream_prep,
        scale_factor = 1.25,
        compute_ks = requested[[1L]],
        compute_cvm = requested[[2L]],
        cpp_kernel = "legacy"
      )
      candidate <- compute_fast_sample_ks_cvm_stats_cpp(
        centered_weight_block = case$centered_weights,
        stream_prep = case$stream_prep,
        scale_factor = 1.25,
        compute_ks = requested[[1L]],
        compute_cvm = requested[[2L]],
        cpp_kernel = "contiguous_double"
      )
      if (requested[[1L]]) {
        ks_difference <- max(abs(legacy$ks - candidate$ks))
        expect_lte(ks_difference, 1e-9)
        expect_lte(
          ks_difference / max(1, max(abs(legacy$ks))),
          1e-12
        )
      } else {
        expect_null(legacy$ks)
        expect_null(candidate$ks)
      }
      if (requested[[2L]]) {
        cvm_difference <- max(abs(legacy$cvm - candidate$cvm))
        expect_lte(
          cvm_difference / max(1, max(abs(legacy$cvm))),
          1e-12
        )
        if (!isTRUE(case$extreme_weights)) {
          expect_lte(cvm_difference, 1e-9)
        }
      } else {
        expect_null(legacy$cvm)
        expect_null(candidate$cvm)
      }
    }
  }
})

test_that("production C++ kernel preserves paired vMF and HvMF bootstraps", {
  compare_runs <- function(run) {
    legacy <- run("legacy")
    candidate <- run("contiguous_double")
    expect_identical(
      legacy$diagnostics$fast_multiplier_cpp_kernel_effective,
      "legacy"
    )
    expect_identical(
      candidate$diagnostics$fast_multiplier_cpp_kernel_effective,
      "contiguous_double"
    )
    expect_equal(
      candidate$bootstrap$statistics$ks,
      legacy$bootstrap$statistics$ks,
      tolerance = 1e-12
    )
    expect_equal(
      candidate$bootstrap$statistics$cvm,
      legacy$bootstrap$statistics$cvm,
      tolerance = 1e-12
    )
    expect_equal(candidate$inference$ks, legacy$inference$ks, tolerance = 1e-12)
    expect_equal(candidate$inference$cvm, legacy$inference$cvm, tolerance = 1e-12)
  }

  set.seed(20260813L)
  vmf_data <- normalize_vmf_data(
    rotasym::r_vMF(14L, mu = c(1, 0, 0), kappa = 4)
  )
  compare_runs(function(kernel) {
    multiplier_bootstrap_vmf(
      data = vmf_data,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = 13L,
      seed = 20260814L,
      n_cores = 1L,
      distance_type = "geodesic",
      unknown_param = "xi",
      bootstrap_method = "fast_multiplier",
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = list(
        derivative_method = "quadrature",
        fast_multiplier_backend = "cpp",
        fast_bootstrap_chunk_size = 13L,
        fast_multiplier_cpp_kernel = kernel
      )
    )
  })

  set.seed(20260815L)
  hvmf_data <- rhvmf_h2_polar(
    14L, mu = c(cosh(0.25), sinh(0.25), 0), kappa = 4
  )
  compare_runs(function(kernel) {
    multiplier_bootstrap_hvmf(
      data = hvmf_data,
      null = list(type = "composite"),
      statistics = c("ks", "cvm"),
      ks_grid = make_sample_unique_distance_ks_grid(),
      B = 13L,
      seed = 20260816L,
      n_cores = 1L,
      unknown_param = "both",
      bootstrap_method = "fast_multiplier",
      fast_multiplier_backend = "cpp",
      keep = list(
        observed_process = FALSE,
        bootstrap_statistics = TRUE,
        bootstrap_thetas = FALSE
      ),
      control = list(
        derivative_method = "quadrature",
        fast_bootstrap_chunk_size = 13L,
        fast_multiplier_cpp_kernel = kernel
      )
    )
  })
})
