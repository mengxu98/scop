#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix augur_subsample_cpp(S4 mat, IntegerVector cols) {
  IntegerVector dims = mat.slot("Dim");
  int n_features = dims[0];
  int n_selected = cols.size();

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

  NumericMatrix dense(n_selected, n_features);
  for (int out_col = 0; out_col < n_selected; ++out_col) {
    int col = selected_cols[out_col];
    for (int ptr = p[col]; ptr < p[col + 1]; ++ptr) {
      dense(out_col, i[ptr]) = x[ptr];
    }
  }

  std::vector<int> kept;
  kept.reserve(n_features);
  for (int gene = 0; gene < n_features; ++gene) {
    bool variable = false;
    double first = dense(0, gene);
    for (int row = 1; row < n_selected; ++row) {
      if (dense(row, gene) != first) {
        variable = true;
        break;
      }
    }
    if (variable) {
      kept.push_back(gene);
    }
  }

  int n_kept = kept.size();
  NumericMatrix out(n_selected, n_kept);
  for (int col = 0; col < n_kept; ++col) {
    int gene = kept[col];
    for (int row = 0; row < n_selected; ++row) {
      out(row, col) = dense(row, gene);
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
