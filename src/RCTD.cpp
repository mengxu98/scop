// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <thisutils/log_message.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct DgCSlots {
  int n_rows;
  int n_cols;
  IntegerVector row_idx;
  IntegerVector col_ptr;
  NumericVector values;
};

DgCSlots get_dgc_slots(S4 mat, const char* arg_name) {
  IntegerVector dims = mat.slot("Dim");
  IntegerVector row_idx = mat.slot("i");
  IntegerVector col_ptr = mat.slot("p");
  NumericVector values = mat.slot("x");
  if (dims.size() != 2) {
    thisutils::log_message_fmt("error", "%s must have two dimensions.",  arg_name);
  }
  if (col_ptr.size() != dims[1] + 1) {
    thisutils::log_message_fmt("error", "%s must be a valid dgCMatrix.",  arg_name);
  }
  if (row_idx.size() != values.size()) {
    thisutils::log_message_fmt("error", "%s has inconsistent sparse matrix slots.",  arg_name);
  }
  return DgCSlots{
    dims[0],
    dims[1],
    row_idx,
    col_ptr,
    values
  };
}

inline bool is_positive_count(double value) {
  return std::isfinite(value) && value > 0.0;
}

void mark_positive_rows(const DgCSlots& mat, std::vector<unsigned char>& has_value) {
  for (int col = 0; col < mat.n_cols; ++col) {
    if ((col & 1023) == 0) {
      Rcpp::checkUserInterrupt();
    }
    for (int ptr = mat.col_ptr[col]; ptr < mat.col_ptr[col + 1]; ++ptr) {
      const int row = mat.row_idx[ptr];
      if (row < 0 || row >= mat.n_rows) {
        thisutils::log_message("dgCMatrix row index is out of bounds.", "error");
      }
      if (is_positive_count(mat.values[ptr])) {
        has_value[static_cast<std::size_t>(row)] = 1;
      }
    }
  }
}

NumericVector kept_column_sums(const DgCSlots& mat, const LogicalVector& keep_rows) {
  NumericVector sums(mat.n_cols);
  for (int col = 0; col < mat.n_cols; ++col) {
    if ((col & 1023) == 0) {
      Rcpp::checkUserInterrupt();
    }
    double total = 0.0;
    for (int ptr = mat.col_ptr[col]; ptr < mat.col_ptr[col + 1]; ++ptr) {
      const int row = mat.row_idx[ptr];
      if (keep_rows[row] && is_positive_count(mat.values[ptr])) {
        total += mat.values[ptr];
      }
    }
    sums[col] = total;
  }
  return sums;
}

}

// [[Rcpp::export]]
List rctd_sparse_quality_cpp(S4 st_counts, S4 ref_counts) {
  DgCSlots st = get_dgc_slots(st_counts, "st_counts");
  DgCSlots ref = get_dgc_slots(ref_counts, "ref_counts");
  if (st.n_rows != ref.n_rows) {
    thisutils::log_message("st_counts and ref_counts must have the same number of rows.", "error");
  }

  std::vector<unsigned char> st_has(static_cast<std::size_t>(st.n_rows), 0);
  std::vector<unsigned char> ref_has(static_cast<std::size_t>(ref.n_rows), 0);
  mark_positive_rows(st, st_has);
  mark_positive_rows(ref, ref_has);

  LogicalVector keep_features(st.n_rows);
  for (int row = 0; row < st.n_rows; ++row) {
    keep_features[row] = st_has[static_cast<std::size_t>(row)] &&
      ref_has[static_cast<std::size_t>(row)];
  }

  return List::create(
    _["keep_features"] = keep_features,
    _["st_numi"] = kept_column_sums(st, keep_features),
    _["ref_numi"] = kept_column_sums(ref, keep_features)
  );
}
