#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static double clipped(double value, double limit) {
  return value > limit ? limit : value;
}

// [[Rcpp::export]]
NumericMatrix scale_sparse_full(S4 sparse_mat,
                                       IntegerVector gene_indices, double scale_max) {
  IntegerVector i_vec = sparse_mat.slot("i");
  IntegerVector p_vec = sparse_mat.slot("p");
  NumericVector x_vec = sparse_mat.slot("x");
  IntegerVector dim_vec = sparse_mat.slot("Dim");
  const int* ip = i_vec.begin();
  const int* pp = p_vec.begin();
  const double* xp = x_vec.begin();
  const int n_total_genes = dim_vec[0];
  const int n_cells = dim_vec[1];
  const int n_sel = gene_indices.size();

  std::vector<int> selected_lookup(n_total_genes, -1);
  for (int out_row = 0; out_row < n_sel; ++out_row) {
    selected_lookup[gene_indices[out_row]] = out_row;
  }

  std::vector<double> center(n_sel, 0.0);
  std::vector<double> inv_sd(n_sel, 1.0);
  std::vector<double> zero_score(n_sel, 0.0);
  std::vector<double> sum(n_sel, 0.0);
  std::vector<double> sumsq(n_sel, 0.0);

  #ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  std::vector<double> thread_sum(static_cast<size_t>(max_threads) * n_sel, 0.0);
  std::vector<double> thread_sumsq(static_cast<size_t>(max_threads) * n_sel, 0.0);
  #pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    double* local_sum = thread_sum.data() + static_cast<size_t>(tid) * n_sel;
    double* local_sumsq = thread_sumsq.data() + static_cast<size_t>(tid) * n_sel;
    #pragma omp for schedule(static)
    for (int col = 0; col < n_cells; ++col) {
      for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
        const int mapped = selected_lookup[ip[pos]];
        if (mapped >= 0) {
          const double value = xp[pos];
          local_sum[mapped] += value;
          local_sumsq[mapped] += value * value;
        }
      }
    }
  }
  for (int tid = 0; tid < max_threads; ++tid) {
    const double* local_sum = thread_sum.data() + static_cast<size_t>(tid) * n_sel;
    const double* local_sumsq = thread_sumsq.data() + static_cast<size_t>(tid) * n_sel;
    for (int row = 0; row < n_sel; ++row) {
      sum[row] += local_sum[row];
      sumsq[row] += local_sumsq[row];
    }
  }
  #else
  for (int col = 0; col < n_cells; ++col) {
    for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
      const int mapped = selected_lookup[ip[pos]];
      if (mapped >= 0) {
        const double value = xp[pos];
        sum[mapped] += value;
        sumsq[mapped] += value * value;
      }
    }
  }
  #endif
  for (int row = 0; row < n_sel; ++row) {
    center[row] = sum[row] / n_cells;
    double variance = (sumsq[row] / n_cells - center[row] * center[row]) *
      n_cells / (n_cells - 1.0);
    if (variance > 0.0) {
      inv_sd[row] = 1.0 / std::sqrt(variance);
    }
    zero_score[row] = clipped(-center[row] * inv_sd[row], scale_max);
  }

  NumericMatrix result = Rcpp::no_init_matrix(n_sel, n_cells);
  double* out = REAL(result);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int col = 0; col < n_cells; ++col) {
    const size_t col_offset = static_cast<size_t>(col) * n_sel;
    std::memcpy(out + col_offset, zero_score.data(), n_sel * sizeof(double));
    for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
      const int row = selected_lookup[ip[pos]];
      if (row >= 0) {
        double value = (xp[pos] - center[row]) * inv_sd[row];
        out[col_offset + row] = clipped(value, scale_max);
      }
    }
  }

  return result;
}
