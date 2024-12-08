#' Fit B-spline curves using Cox-de Boor recursion to get the coefficients that generate weight of corresponding
#'
#' This function calculates b-spline basis function using Cox-de Boor recursion.
#' @param x A numeric vector where the b-spline is evaluated.
#' @param y A numeric vector of the corresponding response value against x
#' @return my_bspline returns a numeric vector of coefficients for fitting the B-spline curve. These coefficients
#'   are the weights applied to each basis function to construct the final B-spline curve
#'   that best fits the input data points (x, y)
#' @examples
#' # Example 1: Simple trigonometric function
#' set.seed(625) # For reproducibility
#' # Generate sample variable x and its response y
#' x <- seq(0, 10, length.out = 20)
#' y <- sin(x) + 2*cos(x)
#' # Calculate the coefficients (for weight calculation) for b-spline basis functions
#' coef1 <- my_bspline(x, y)
#'
#' # Example 2: Fitting and prediction
#' # Generate example data
#' set.seed(625)
#' x <- seq(0, 1, length.out = 10)
#' y <- sin(2*pi*x) + rnorm(10, 0, 0.1)
#'
#' # Calculate and print the coefficients
#' coef <- my_bspline(x, y)
#' cat("B-spline coefficients: ", round(coef, 4))
#'
#' # Predict at new points
#' x_new <- seq(0, 1, length.out = 100)
#' y_pred <- numeric(length(x_new))
#'
#' # Calculate predictions using existing basis functions
#' order <- 4
#' # Creates equally spaced knots between min(x) and max(x)
#' internal_knots <- seq(min(x), max(x),
#'                      length.out = length(x) - order + 2)
#' # Adds repeated knots at boundaries to avoid poor behavior at boundaries
#' knots <- c(rep(min(x), order-1), internal_knots, rep(max(x), order-1))
#'
#' for(i in 1:length(x_new)) {
#'   y_pred[i] <- sum(coef * sapply(1:length(coef), function(j)
#'     bspline_basis(x_new[i], j, order-1, knots)))
#' }
#'
#' # Plot the results
#' plot(x, y, col = "blue", pch = 16)
#' lines(x_new, y_pred, col = "red")
#' legend("topright", c("Data", "B-spline fit"),
#'        col = c("blue", "red"),
#'        pch = c(16, NA),
#'        lty = c(NA, 1))
#' @importFrom graphics plot lines legend
#' @importFrom stats rnorm
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





