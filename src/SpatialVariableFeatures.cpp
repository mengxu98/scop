#include <Rcpp.h>

#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector spatial_variable_score_dense_cpp(
    NumericMatrix expr,
    IntegerVector from,
    IntegerVector to,
    std::string method) {
  const int n_features = expr.nrow();
  const int n_spots = expr.ncol();
  if (from.size() != to.size() || n_spots < 1 ||
      (method != "moran" && method != "geary")) {
    stop("invalid spatial variable feature inputs");
  }
  for (int edge = 0; edge < from.size(); ++edge) {
    if (from[edge] < 1 || from[edge] > n_spots ||
        to[edge] < 1 || to[edge] > n_spots) {
      stop("spatial graph edge index is out of bounds");
    }
  }

  // R matrices are column-major. Accumulating all features while each spot
  // or graph-edge column is resident avoids striding by n_features for every
  // individual feature, while retaining the reference iteration order.
  const double* values = expr.begin();
  std::vector<int> n_finite(n_features, 0);
  std::vector<int> n_edges(n_features, 0);
  std::vector<double> sums(n_features, 0.0);
  std::vector<double> means(n_features, NA_REAL);
  std::vector<double> denom(n_features, 0.0);
  std::vector<double> numerator(n_features, 0.0);

  for (int spot = 0; spot < n_spots; ++spot) {
    const double* column = values + static_cast<size_t>(spot) * n_features;
    for (int feature = 0; feature < n_features; ++feature) {
      const double value = column[feature];
      if (std::isfinite(value)) {
        ++n_finite[feature];
        sums[feature] += value;
      }
    }
  }
  for (int feature = 0; feature < n_features; ++feature) {
    if (n_finite[feature] >= 3) {
      means[feature] = sums[feature] / static_cast<double>(n_finite[feature]);
    }
  }
  for (int spot = 0; spot < n_spots; ++spot) {
    const double* column = values + static_cast<size_t>(spot) * n_features;
    for (int feature = 0; feature < n_features; ++feature) {
      const double value = column[feature];
      if (std::isfinite(value) && n_finite[feature] >= 3) {
        const double centered = value - means[feature];
        denom[feature] += centered * centered;
      }
    }
  }
  for (int edge = 0; edge < from.size(); ++edge) {
    const double* left_column = values +
      static_cast<size_t>(from[edge] - 1) * n_features;
    const double* right_column = values +
      static_cast<size_t>(to[edge] - 1) * n_features;
    for (int feature = 0; feature < n_features; ++feature) {
      if (n_finite[feature] < 3 || !std::isfinite(denom[feature]) ||
          denom[feature] <= 0.0) continue;
      const double left = left_column[feature];
      const double right = right_column[feature];
      if (!std::isfinite(left) || !std::isfinite(right)) continue;
      ++n_edges[feature];
      if (method == "moran") {
        numerator[feature] += (left - means[feature]) *
          (right - means[feature]);
      } else {
        const double difference = left - right;
        numerator[feature] += difference * difference;
      }
    }
  }

  NumericVector statistic(n_features, NA_REAL);
  for (int feature = 0; feature < n_features; ++feature) {
    if (n_finite[feature] < 3 || !std::isfinite(denom[feature]) ||
        denom[feature] <= 0.0 || n_edges[feature] == 0) continue;
    if (method == "moran") {
      statistic[feature] = static_cast<double>(n_finite[feature]) /
        static_cast<double>(n_edges[feature]) * numerator[feature] /
        denom[feature];
    } else {
      statistic[feature] = static_cast<double>(n_finite[feature] - 1) /
        (2.0 * static_cast<double>(n_edges[feature])) * numerator[feature] /
        denom[feature];
    }
  }
  return statistic;
}
