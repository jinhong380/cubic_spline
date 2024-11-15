#' Evaluate a B-spline basis function
#'
#' This function calculates B-spline basis function using Cox-de Boor recursion.
#'
#' @param x A numeric value where the B-spline is evaluated.
#' @param i The index of the basis function.
#' @param order The order of the B-spline.
#' @param knots A numeric vector of knots for the B-spline.
#' @return A numeric value representing the evaluated B-spline basis function.
#' @name bspline_basis
#' @export


my_bspline <- function(x0, y0) {

  # Function to evaluate single B-spline basis function using Cox-de Boor recursion
  bspline_basis <- function(x, i, order, knots) {
    if (order == 0) {
      # For order 0, return indicator function
      return(ifelse(x >= knots[i] & x < knots[i + 1], 1, 0))
    }

    # Handle division by zero cases
    coef1 <- ifelse(knots[i + order] == knots[i], 0, (x - knots[i]) / (knots[i + order] - knots[i]))

    coef2 <- ifelse (knots[i + order + 1] == knots[i + 1], 0,
      (knots[i + order + 1] - x) / (knots[i + order + 1] - knots[i + 1]))

    # Recursive calculation
    term1 <- coef1 * bspline_basis(x, i, order - 1, knots)
    term2 <- coef2 * bspline_basis(x, i + 1, order - 1, knots)

    return(term1 + term2)
  }

  m <- length(x0)
  order <- 4       # as the bspline use one additional degree compare to cubic spline

  # Creating knot sequence:
  internal_knots <- seq(min(x0), max(x0), length.out = m - order + 2)
  knots <- c(rep(min(x0), order-1), internal_knots, rep(max(x0), order-1))

  # Initialize a basis matrix
  mat_b <- matrix(0, nrow = length(x0), ncol = m)

  # Fill the basis matrix mat_b
  for (i in 1:length(x0)) {
    for (j in 1:(m)) {
      mat_b[i,j] <- bspline_basis(x0[i], j, order-1, knots) # Evaluate basis function
    }
  }

  mat_b[length(x0), m] <- 1   #Set rightmost point to cope with boundary issue

  #Solve for coefficients using QR decomposition for better numerical stability
  coef <- qr.solve(mat_b, y0)

  return(coef) # Return the coefficients
}


# Example use of the function
#x0 <- seq(0, 10, length.out = 20)
#y0 <- sin(x0) + rnorm(length(x0), 0, 0.1)
#coef <- my_bspline(x0, y0)





