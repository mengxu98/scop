// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>

#include <unordered_set>
#include <vector>

using namespace Rcpp;

namespace {

struct RowUniqueValues {
  std::unordered_set<double> values;
  bool has_na = false;
  bool has_nan = false;
};

inline void add_value(RowUniqueValues& row, const double value) {
  if (R_IsNA(value)) {
    row.has_na = true;
  } else if (ISNAN(value)) {
    row.has_nan = true;
  } else {
    row.values.insert(value);
  }
}

inline int unique_count(const RowUniqueValues& row) {
  return static_cast<int>(row.values.size()) +
    static_cast<int>(row.has_na) + static_cast<int>(row.has_nan);
}

}  // namespace

// Count distinct values in every row of a dense numeric matrix.  NA and NaN
// remain distinct, matching `length(unique(row))` in base R.
// [[Rcpp::export]]
IntegerVector dynamic_row_unique_counts_dense_cpp(NumericMatrix x) {
  const int n_rows = x.nrow();
  const int n_cols = x.ncol();
  IntegerVector out(n_rows);

  for (int row_i = 0; row_i < n_rows; ++row_i) {
    RowUniqueValues values;
    for (int col_i = 0; col_i < n_cols; ++col_i) {
      add_value(values, x(row_i, col_i));
    }
    out[row_i] = unique_count(values);
  }
  return out;
}

// Count distinct values in each row of a dgCMatrix without materialising it.
// The final implicit-zero increment intentionally follows the existing R
// implementation: it is based on the number of stored entries, including an
// explicitly stored zero should one be present.
// [[Rcpp::export]]
IntegerVector dynamic_row_unique_counts_sparse_cpp(S4 x) {
  IntegerVector dims = x.slot("Dim");
  if (dims.size() != 2) {
    stop("x must have a two-dimensional Dim slot");
  }
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  IntegerVector row_index = x.slot("i");
  NumericVector values = x.slot("x");
  if (row_index.size() != values.size()) {
    stop("sparse matrix i and x slots must have the same length");
  }

  std::vector<RowUniqueValues> row_values(n_rows);
  std::vector<int> row_stored(n_rows, 0);
  for (R_xlen_t value_i = 0; value_i < values.size(); ++value_i) {
    const int row_i = row_index[value_i];
    if (row_i < 0 || row_i >= n_rows) {
      stop("sparse matrix contains an out-of-bounds row index");
    }
    add_value(row_values[row_i], values[value_i]);
    ++row_stored[row_i];
  }

  IntegerVector out(n_rows);
  for (int row_i = 0; row_i < n_rows; ++row_i) {
    out[row_i] = unique_count(row_values[row_i]) +
      static_cast<int>(row_stored[row_i] < n_cols);
  }
  return out;
}
