#' Evaluate single B-spline basis function using Cox-de Boor recursion
#'
#' @param x Numeric value at which to evaluate the basis function
#' @param i Integer index of the basis function
#' @param order Integer order of the B-spline
#' @param knots Numeric vector of knot positions
#' @return Numeric value of the basis function at x
#' @keywords internal
#'
bspline_basis <- function(x, i, order, knots) {
  if (order == 0) {
    # For order 0, return indicator function
    return(ifelse(x >= knots[i] & x < knots[i + 1], 1, 0))
  }
  # Handle division by zero cases
  coef1 <- ifelse(knots[i + order] == knots[i], 0, (x - knots[i]) / (knots[i + order] - knots[i]))
  coef2 <- ifelse(knots[i + order + 1] == knots[i + 1], 0, (knots[i + order + 1] - x) / (knots[i + order + 1] - knots[i + 1]))
  # Recursive calculation
  term1 <- coef1 * bspline_basis(x, i, order - 1, knots)
  term2 <- coef2 * bspline_basis(x, i + 1, order - 1, knots)
  return(term1 + term2)
}


#' Evaluate a b-spline basis function
#'
#' This function calculates b-spline basis function using Cox-de Boor recursion.
#'
#' @param x A numeric value where the b-spline is evaluated.
#' @param i The index of the basis function.
#' @param order The order of the B-spline.
#' @param knots A numeric vector of knots for the b-spline.
#' @return A numeric value representing the evaluated b-spline basis function used for constructing the resulted whole splines
#' @name bspline_basis
#' @export

my_bspline <- function(x, y) {
  m <- length(x)
  order <- 4  # as the bspline uses one additional degree compare to cubic spline

  # Creating knot sequence:
  internal_knots <- seq(min(x), max(x), length.out = m - order + 2)
  knots <- c(rep(min(x), order-1), internal_knots, rep(max(x), order-1))

  # Initialize a basis matrix
  mat_b <- matrix(0, nrow = length(x), ncol = m)

  # Fill the basis matrix mat_b
  for (i in 1:length(x)) {
    for (j in 1:(m)) {
      mat_b[i,j] <- bspline_basis(x[i], j, order-1, knots)
    }
  }

  mat_b[length(x), m] <- 1  #Set rightmost point to cope with boundary issue

  coef <- solve(mat_b, y)
  return(coef)
}


# Example use of the function
#x <- seq(0, 10, length.out = 20)
#y <- sin(x) + rnorm(length(x), 0, 0.1)
#coef <- my_bspline(x, y)





