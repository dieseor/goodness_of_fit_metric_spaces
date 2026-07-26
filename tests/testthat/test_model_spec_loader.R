library(testthat)

oldwd <- setwd(normalizePath(file.path("..", "..")))
on.exit(setwd(oldwd), add = TRUE)

source(file.path("bootstrap", "model_specs.R"))

test_that("the model-spec entry point loads every family adapter", {
  constructors <- c(
    "make_normal_spec",
    "make_mvnormal_spec",
    "make_logistic_gaussian_spec",
    "make_vmf_spec",
    "make_jp_spec",
    "make_hvmf_spec",
    "make_spherical_cauchy_spec",
    "make_beta_mixture2_spec",
    "make_uniform_beta_mixture_spec",
    "make_logitnormal_mixture2_spec",
    "make_axial_truncnorm_mixture2_spec",
    "make_cardioid_spec",
    "make_small_circle_spec",
    "make_watson_spec",
    "make_small_circle_symmetric_mixture2_spec",
    "make_small_circle_weighted_mixture2_spec"
  )

  expect_true(all(vapply(constructors, exists, logical(1), mode = "function")))
})
