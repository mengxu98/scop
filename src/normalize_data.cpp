#include <Rcpp.h>
#include <cmath>
#include "thread_utils.h"
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

// grain_size is retained for API compatibility but is no longer used;
// column-level parallelism is now handled by OpenMP (each column is
// independent: it normalises only its own non-zero entries in @x).
// [[Rcpp::export]]
void log_normalize_dgc(S4 mat, double scale_factor, int grain_size = 100, int n_threads = 0) {
  NumericVector x = mat.slot("x");
  IntegerVector p = mat.slot("p");
  const int columns = p.size() - 1;
  double* values = REAL(x);
  const int* colptr = INTEGER(p);
  const int threads = omp_thread_count(n_threads, columns);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
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
