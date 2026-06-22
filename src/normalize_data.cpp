#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static inline double column_total(const double* values, int begin, int end) {
  double total = 0.0;
  for (int offset = begin; offset < end; ++offset) {
    total += values[offset];
  }
  return total;
}

static inline void rescale_column(double* values, int begin, int end,
                                  double multiplier) {
  for (int offset = begin; offset < end; ++offset) {
    values[offset] = std::log1p(values[offset] * multiplier);
  }
}

// [[Rcpp::export]]
void log_normalize_dgc(S4 mat, double scale_factor, int grain_size = 100) {
  NumericVector x = mat.slot("x");
  IntegerVector p = mat.slot("p");
  const int columns = p.size() - 1;
  double* values = REAL(x);
  const int* colptr = INTEGER(p);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, grain_size)
  #endif
  for (int col = 0; col < columns; ++col) {
    const int first = colptr[col];
    const int last = colptr[col + 1];
    const double total = column_total(values, first, last);
    if (total > 0.0) {
      rescale_column(values, first, last, scale_factor / total);
    }
  }
}
