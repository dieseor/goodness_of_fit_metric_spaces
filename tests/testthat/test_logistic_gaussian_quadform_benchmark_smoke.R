library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

suppressWarnings(source(file.path("tests", "benchmark_logistic_gaussian_quadform_utils.R")))

test_that("quadratic-form benchmark smoke run produces expected artifacts", {
  output_dir <- file.path(tempdir(), "lg_quadform_benchmark_smoke")
  result <- run_logistic_gaussian_quadform_benchmark(
    output_dir = output_dir,
    n_cores = 1L,
    random_n = 8L,
    max_cases = 24L,
    random_seed = 321L,
    save_plots = FALSE
  )

  expect_true(nrow(result$cases) > 0)
  expect_true(nrow(result$results) > 0)
  expect_true(all(c("farebrother", "davies", "imhof", "sphunif_sw", "sphunif_hbe") %in% unique(result$results$method)))
  expect_true(file.exists(file.path(output_dir, "quadform_cases.csv")))
  expect_true(file.exists(file.path(output_dir, "quadform_results.csv")))
  expect_true(file.exists(file.path(output_dir, "quadform_method_summary.csv")))
  expect_true(file.exists(file.path(output_dir, "farebrother_slow_cases.csv")))
  expect_true(file.exists(file.path(output_dir, "recommendation.md")))
  expect_true(file.exists(file.path(output_dir, "quadform_benchmark.rds")))
})
