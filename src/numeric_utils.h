#ifndef SCOP_NUMERIC_UTILS_H
#define SCOP_NUMERIC_UTILS_H

#include <Rcpp.h>
#include <vector>

// Convert an R numeric matrix to a row-major float buffer (used by the
// kNN/UMAP C++ kernels that consume flat float arrays).
inline std::vector<float> matrix_as_row_float(Rcpp::NumericMatrix x) {
  const int rows = x.nrow();
  const int cols = x.ncol();
  std::vector<float> out(static_cast<size_t>(rows) * cols);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      out[static_cast<size_t>(row) * cols + col] =
        static_cast<float>(x(row, col));
    }
  }
  return out;
}

#endif // SCOP_NUMERIC_UTILS_H
