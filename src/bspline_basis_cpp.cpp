#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double bspline_basis_cpp(double x, int i, int order, NumericVector knots) {
  if (order == 0) {
    return (x >= knots[i] && x < knots[i + 1]) ? 1.0 : 0.0;
  }

  double coef1 = 0.0;
  if (knots[i + order] != knots[i]) {
    coef1 = (x - knots[i]) / (knots[i + order] - knots[i]);
  }

  double coef2 = 0.0;
  if (knots[i + order + 1] != knots[i + 1]) {
    coef2 = (knots[i + order + 1] - x) / (knots[i + order + 1] - knots[i + 1]);
  }

  return coef1 * bspline_basis_cpp(x, i, order - 1, knots) +
    coef2 * bspline_basis_cpp(x, i + 1, order - 1, knots);
}
