# Mathematical and numerical helpers: simplex.

small_circle_logistic_bounded <- function(x, upper = 1 - 1e-6) {
  upper / (1 + exp(-x))
}
small_circle_inverse_logistic_bounded <- function(y, upper = 1 - 1e-6) {
  y <- min(max(y, 1e-8), upper - 1e-8)
  stats::qlogis(y / upper)
}
