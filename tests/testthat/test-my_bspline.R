# tests/testthat/test-bspline.R

# Tests are based on the fda package

test_that("B-spline coefficients match fda package for Runge function", {
  # Test case 1: Runge function
  x0 <- seq(-1, 1, length.out = 15)
  f <- function(x) 1/(1 + 25*x^2)
  y0 <- f(x0)
  my_coef <- my_bspline(x0, y0)

  # From fda package
  library(fda)
  degree <- 3
  n_basis <- length(x0)
  basis <- create.bspline.basis(
    rangeval = c(min(x0), max(x0)),
    nbasis = n_basis,
    norder = degree + 1
  )
  basis_matrix <- eval.basis(x0, basis)
  fit_fda <- lm(y0 ~ basis_matrix - 1)
  fda_coef <- as.numeric(fit_fda$coefficients)

  # Test equality with some tolerance
  expect_equal(my_coef, fda_coef, tolerance = 1e-6)
})


test_that("B-spline coefficients match fda package for Random generated values", {
  set.seed(625)
  x0 <- seq(-1, 1, length.out = 15)
  y0 <- rnorm(length(x0), 0, 0.1)
  my_coef <- my_bspline(x0, y0)

  # From fda package
  library(fda)
  degree <- 3
  n_basis <- length(x0)
  basis <- create.bspline.basis(
    rangeval = c(min(x0), max(x0)),
    nbasis = n_basis,
    norder = degree + 1
  )
  basis_matrix <- eval.basis(x0, basis)
  fit_fda <- lm(y0 ~ basis_matrix - 1)
  fda_coef <- as.numeric(fit_fda$coefficients)

  expect_equal(my_coef, fda_coef, tolerance = 1e-6)
})

test_that("B-spline coefficients match fda package for trigonometric function", {
  x0 <- seq(0, 10, length.out = 20)
  y0 <- sin(x0) + cos(x0)
  my_coef <- my_bspline(x0, y0)

  # Get coefficients from fda package
  library(fda)
  degree <- 3
  n_basis <- length(x0)
  basis <- create.bspline.basis(
    rangeval = c(min(x0), max(x0)),
    nbasis = n_basis,
    norder = degree + 1
  )
  basis_matrix <- eval.basis(x0, basis)
  fit_fda <- lm(y0 ~ basis_matrix - 1)
  fda_coef <- as.numeric(fit_fda$coefficients)
  expect_equal(my_coef, fda_coef, tolerance = 1e-6)
})
