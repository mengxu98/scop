#include <Rcpp.h>

#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix augur_subsample_cpp(S4 mat, IntegerVector cols) {
  IntegerVector dims = mat.slot("Dim");
  int n_features = dims[0];
  int n_selected = cols.size();
  if (n_selected < 1) {
    stop("At least one column must be selected");
  }

  IntegerVector p = mat.slot("p");
  IntegerVector i = mat.slot("i");
  NumericVector x = mat.slot("x");
  List dimnames = mat.slot("Dimnames");
  CharacterVector feature_names = dimnames[0];
  CharacterVector cell_names = dimnames[1];

  std::vector<int> selected_cols(n_selected);
  for (int k = 0; k < n_selected; ++k) {
    int col = cols[k] - 1;
    if (col < 0 || col >= dims[1]) {
      stop("Column index out of bounds");
    }
    selected_cols[k] = col;
  }

  // Identify variable genes directly from the sparse columns. This avoids
  // materialising the full selected-cell-by-feature matrix before allocating
  // the final dense classifier input.
  std::vector<int> stored(n_features, 0);
  std::vector<bool> seen(n_features, false);
  std::vector<bool> variable(n_features, false);
  std::vector<double> first(n_features, 0.0);
  for (int out_col = 0; out_col < n_selected; ++out_col) {
    int col = selected_cols[out_col];
    for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
      int gene = i[ptr];
      double value = x[ptr];
      ++stored[gene];
      if (!seen[gene]) {
        seen[gene] = true;
        first[gene] = value;
      } else if (value != first[gene]) {
        variable[gene] = true;
      }
    }
  }
  std::vector<int> kept;
  kept.reserve(n_features);
  for (int gene = 0; gene < n_features; ++gene) {
    if (stored[gene] < n_selected && seen[gene] && first[gene] != 0.0) {
      variable[gene] = true;
    }
    if (variable[gene]) {
      kept.push_back(gene);
    }
  }

  int n_kept = kept.size();
  NumericMatrix out(n_selected, n_kept);
  std::vector<int> output_col(n_features, -1);
  for (int col = 0; col < n_kept; ++col) {
    output_col[kept[col]] = col;
  }
  for (int out_col = 0; out_col < n_selected; ++out_col) {
    int col = selected_cols[out_col];
    for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
      int result_col = output_col[i[ptr]];
      if (result_col >= 0) {
        out(out_col, result_col) = x[ptr];
      }
    }
  }

  CharacterVector row_names(n_selected);
  for (int row = 0; row < n_selected; ++row) {
    row_names[row] = cell_names[selected_cols[row]];
  }
  CharacterVector col_names(n_kept);
  for (int col = 0; col < n_kept; ++col) {
    col_names[col] = feature_names[kept[col]];
  }
  out.attr("dimnames") = List::create(row_names, col_names);
  return out;
}
