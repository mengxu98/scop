// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>

using namespace Rcpp;

// Calculate the per-cell geometric mean of a non-negative dgCMatrix without
// extracting every sparse column into R.  The caller retains the legacy R
// path for negative or non-finite inputs, where base R's NA/NaN semantics are
// part of the public behaviour.
// [[Rcpp::export]]
NumericVector feature_cor_geometric_mean_sparse_cpp(S4 x, bool log_normalized) {
  IntegerVector dims = x.slot("Dim");
  if (dims.size() != 2) {
    stop("x must have a two-dimensional Dim slot");
  }
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  IntegerVector pointers = x.slot("p");
  NumericVector values = x.slot("x");
  if (pointers.size() != n_cols + 1) {
    stop("sparse matrix p slot is incompatible with Dim");
  }

  NumericVector out(n_cols);
  for (int col_i = 0; col_i < n_cols; ++col_i) {
    const int begin = pointers[col_i];
    const int end = pointers[col_i + 1];
    if (begin < 0 || end < begin || end > values.size()) {
      stop("sparse matrix contains invalid column pointers");
    }
    if (end - begin < n_rows) {
      out[col_i] = 0.0;
      continue;
    }

    double log_sum = 0.0;
    bool has_zero = false;
    for (int value_i = begin; value_i < end; ++value_i) {
      const double value = values[value_i];
      if (value == 0.0) {
        has_zero = true;
        break;
      }
      log_sum += log_normalized ? std::log(std::expm1(value)) : std::log(value);
    }
    if (has_zero) {
      out[col_i] = 0.0;
      continue;
    }

    const double mean_log = log_sum / static_cast<double>(n_rows);
    out[col_i] = log_normalized ? std::log1p(std::exp(mean_log)) : std::exp(mean_log);
  }
  return out;
}
