# Run after installing dpgof, from any working directory.
stopifnot(requireNamespace("dpgof", quietly = TRUE))
expected <- c("multiplier_bootstrap_normal", "multiplier_bootstrap_mvnormal",
              "multiplier_bootstrap_vmf", "multiplier_bootstrap_hvmf",
              "multiplier_bootstrap_jp", "multiplier_bootstrap_cardioid",
              "multiplier_bootstrap_logistic_gaussian",
              "multiplier_bootstrap_small_circle", "multiplier_bootstrap_watson")
stopifnot(all(expected %in% getNamespaceExports("dpgof")))
set.seed(11)
x <- rnorm(12)
null <- list(type = "simple", theta = list(mu = 0, sigma = 1))
r <- dpgof::multiplier_bootstrap_normal(x, null, statistics = "ks", B = 19,
                                         seed = 22, distance_profile_backend = "r")
c <- dpgof::multiplier_bootstrap_normal(x, null, statistics = "ks", B = 19,
                                         seed = 22, distance_profile_backend = "cpp")
stopifnot(isTRUE(all.equal(r$inference$ks$p_value, c$inference$ks$p_value,
                           tolerance = 1e-10)))
set.seed(3)
y <- dpgof:::r_sph_car(12, c(0, 0, 1), 0.4, 2)
cardioid <- dpgof::multiplier_bootstrap_cardioid(
  y, null = list(type = "simple", theta = list(mu = c(0, 0, 1), rho = 0.4, k = 2)),
  k = 2, statistics = "ks", B = 9, seed = 7)
stopifnot(is.finite(cardioid$inference$ks$p_value))
stopifnot(inherits(dpgof:::make_sunspots_joint_time_space_spec("shared"), "gof_model_spec"))
cat("Package smoke checks passed.\n")
