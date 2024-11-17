# cubicSpline

<!-- badges: start -->

[![R-CMD-check](https://github.com/jinhong380/cubic_spline/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jinhong380/cubic_spline/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/jinhong380/cubic_spline/graph/badge.svg)](https://app.codecov.io/gh/jinhong380/cubic_spline)

<!-- badges: end -->

The goal of cubicSpline is to implement the b-spline with the Cox-de Boor recursion formula. The package provides tools for both the coefficient calculated for each b-spline basis functions and their corresponding basis function (at x) numeric values.

## Installation

You can install the cubicSpline package from [GitHub](https://github.com/) with:

``` r
devtools::install_github("jinhong380/cubic_spline") #, build_vignettes = T) add this part to include vignettes
```

## Example

This is a basic example of using my_bspline to obtain coefficients of the b-spline:

``` r
library(cubicSpline)
set.seed(625) # For reproducibility

# Generate sample variable x and its response y
x <- seq(0, 10, length.out = 20)
y <- sin(x) + 2*cos(x)

# Calculate the coefficients (for weight calculation) for b-spline basis functions
my_bspline(x, y)
```

This is an application for point prediction with my_bsline and bspline_basis

``` r
# Generate example data
set.seed(625)
x <- seq(0, 1, length.out = 10)
y <- sin(2*pi*x) + rnorm(10, 0, 0.1)

# Get coefficients
coef <- my_bspline(x, y)

# Predict at new points
x_new <- seq(0, 1, length.out = 100)
y_pred <- numeric(length(x_new))

# Calculate predictions using existing basis functions
order <- 4
# Creates equally spaced knots between min(x) and max(x)
internal_knots <- seq(min(x), max(x), length.out = length(x) - order + 2)
# Adds repeated knots at boundaries to avoid poor behavior at boundaries
knots <- c(rep(min(x), order-1), internal_knots, rep(max(x), order-1))
for(i in 1:length(x_new)) {
  y_pred[i] <- sum(coef * sapply(1:length(coef), function(j) 
    bspline_basis(x_new[i], j, order-1, knots)))
}
```

## Development

This package uses:

-   testthat for unit testing

-   roxygen2 for documentation

## License

GNU-3

## Contact

[Jintong\@umich.edu](mailto:Jintong@umich.edu){.email}
