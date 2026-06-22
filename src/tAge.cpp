#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
NumericVector tage_elastic_net_predict_cpp(
    NumericMatrix expr,
    IntegerVector feature_match,
    NumericVector imputer,
    NumericVector center,
    NumericVector scale,
    NumericVector coef,
    double intercept) {
  const int n_features = feature_match.size();
  const int n_samples = expr.ncol();
  if (imputer.size() != n_features || center.size() != n_features || coef.size() != n_features) {
    stop("feature_match, imputer, center, and coef must have the same length");
  }
  const bool use_scale = scale.size() > 0;
  if (use_scale && scale.size() != n_features) {
    stop("scale must be empty or have the same length as coef");
  }

  NumericVector prediction(n_samples, intercept);
  for (int j = 0; j < n_features; ++j) {
    const int row = feature_match[j];
    const bool present = row != NA_INTEGER && row > 0 && row <= expr.nrow();
    const double c = center[j];
    const double s = use_scale ? scale[j] : 1.0;
    const double beta = coef[j];
    if (beta == 0.0) {
      continue;
    }
    for (int i = 0; i < n_samples; ++i) {
      double value = present ? expr(row - 1, i) : imputer[j];
      if (NumericVector::is_na(value)) {
        value = imputer[j];
      }
      prediction[i] += ((value - c) / s) * beta;
    }
  }
  return prediction;
}
