#' B-spline Basis Function
#'
#' @export
#' @param x Numeric vector at which to evaluate the basis function
#' @param i Integer index of the basis function
#' @param order Integer order of the B-spline, default 4
#' @param knots Numeric vector of knot positions
#' @return A numeric value representing how much this particular piece of the B-spline curve contributes at point x.
#' @examples
#' # Create knots and evaluation points
#' knots <- c(0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4)
#' x <- seq(0, 4, length.out = 100)
#' basis1 <- sapply(x, bspline_basis, i = 1, order = 3, knots = knots)
#' print(basis1)

bspline_basis <- function(x, i, order = 4, knots) {
  # For order 0, return indicator function
  if (order == 0) {
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

