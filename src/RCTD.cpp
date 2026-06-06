// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
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
    stop("%s must have two dimensions.", arg_name);
  }
  if (col_ptr.size() != dims[1] + 1) {
    stop("%s must be a valid dgCMatrix.", arg_name);
  }
  if (row_idx.size() != values.size()) {
    stop("%s has inconsistent sparse matrix slots.", arg_name);
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
        stop("dgCMatrix row index is out of bounds.");
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

NumericMatrix normalize_weights_copy(NumericMatrix weights) {
  NumericMatrix out = clone(weights);
  const int n_rows = out.nrow();
  const int n_cols = out.ncol();

  for (int row = 0; row < n_rows; ++row) {
    double total = 0.0;
    for (int col = 0; col < n_cols; ++col) {
      double value = out(row, col);
      if (!std::isfinite(value) || value < 0.0) {
        value = 0.0;
        out(row, col) = 0.0;
      }
      total += value;
    }
    if (std::isfinite(total) && total > 0.0) {
      for (int col = 0; col < n_cols; ++col) {
        out(row, col) /= total;
      }
    } else {
      for (int col = 0; col < n_cols; ++col) {
        out(row, col) = 0.0;
      }
    }
  }

  out.attr("dimnames") = weights.attr("dimnames");
  return out;
}

List metadata_from_weights(NumericMatrix weights, CharacterVector all_spots) {
  List dimnames = weights.attr("dimnames");
  if (dimnames.size() < 2 || dimnames[0] == R_NilValue || dimnames[1] == R_NilValue) {
    stop("weights must have row and column names.");
  }
  CharacterVector weight_spots = dimnames[0];
  CharacterVector cell_types = dimnames[1];
  const int n_spots = all_spots.size();
  const int n_types = weights.ncol();
  if (n_types == 0) {
    stop("weights must contain at least one cell type.");
  }

  std::unordered_map<std::string, int> weight_index;
  weight_index.reserve(static_cast<std::size_t>(weight_spots.size()));
  for (int row = 0; row < weight_spots.size(); ++row) {
    if (CharacterVector::is_na(weight_spots[row])) {
      continue;
    }
    std::string key = as<std::string>(weight_spots[row]);
    if (weight_index.find(key) == weight_index.end()) {
      weight_index[key] = row;
    }
  }

  NumericMatrix full_weights(n_spots, n_types);
  std::fill(full_weights.begin(), full_weights.end(), NA_REAL);
  CharacterVector dominant(n_spots);
  NumericVector max_prop(n_spots, NA_REAL);
  for (int spot = 0; spot < n_spots; ++spot) {
    dominant[spot] = NA_STRING;
  }

  for (int spot = 0; spot < n_spots; ++spot) {
    if (CharacterVector::is_na(all_spots[spot])) {
      continue;
    }
    const std::string key = as<std::string>(all_spots[spot]);
    const std::unordered_map<std::string, int>::const_iterator hit =
      weight_index.find(key);
    if (hit == weight_index.end()) {
      continue;
    }

    const int source_row = hit->second;
    double row_max = -std::numeric_limits<double>::infinity();
    int max_index = 0;
    bool row_has_value = false;
    for (int col = 0; col < n_types; ++col) {
      const double value = weights(source_row, col);
      full_weights(spot, col) = value;
      const bool is_missing = NumericVector::is_na(value);
      const double comparable = is_missing ? 0.0 : value;
      row_has_value = row_has_value || !is_missing;
      if (col == 0 || comparable > row_max) {
        row_max = comparable;
        max_index = col;
      }
    }
    if (row_has_value) {
      max_prop[spot] = row_max;
      if (row_max > 0.0) {
        dominant[spot] = cell_types[max_index];
      }
    }
  }

  full_weights.attr("dimnames") = List::create(all_spots, cell_types);
  return List::create(
    _["weights"] = full_weights,
    _["dominant"] = dominant,
    _["max_prop"] = max_prop
  );
}

} // namespace

// [[Rcpp::export]]
List rctd_sparse_quality_cpp(S4 st_counts, S4 ref_counts) {
  DgCSlots st = get_dgc_slots(st_counts, "st_counts");
  DgCSlots ref = get_dgc_slots(ref_counts, "ref_counts");
  if (st.n_rows != ref.n_rows) {
    stop("st_counts and ref_counts must have the same number of rows.");
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

// [[Rcpp::export]]
NumericMatrix rctd_normalize_weights_cpp(NumericMatrix weights) {
  return normalize_weights_copy(weights);
}

// [[Rcpp::export]]
List rctd_metadata_cpp(NumericMatrix weights, CharacterVector all_spots) {
  return metadata_from_weights(weights, all_spots);
}

// [[Rcpp::export]]
List rctd_finalize_weights_cpp(NumericMatrix weights, CharacterVector all_spots) {
  NumericMatrix normalized = normalize_weights_copy(weights);
  List metadata = metadata_from_weights(normalized, all_spots);
  NumericMatrix full_weights = metadata["weights"];
  CharacterVector dominant = metadata["dominant"];
  NumericVector max_prop = metadata["max_prop"];
  return List::create(
    _["weights"] = normalized,
    _["full_weights"] = full_weights,
    _["dominant"] = dominant,
    _["max_prop"] = max_prop
  );
}
