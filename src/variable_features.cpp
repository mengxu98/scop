#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

struct RowAccumulation {
  std::vector<double> sum;
  std::vector<double> sumsq;
  std::vector<int> count;

  explicit RowAccumulation(int rows) : sum(rows, 0.0), sumsq(rows, 0.0),
    count(rows, 0) {}
};

static RowAccumulation accumulate_rows(const int* colptr, const int* row_index,
                                       const double* values, int rows,
                                       int columns) {
  RowAccumulation out(rows);
#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  if (max_threads > 1 && columns >= max_threads * 64) {
    std::vector<RowAccumulation> local;
    local.reserve(max_threads);
    for (int t = 0; t < max_threads; ++t) {
      local.push_back(RowAccumulation(rows));
    }
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      RowAccumulation& acc = local[tid];
#pragma omp for schedule(static)
      for (int col = 0; col < columns; ++col) {
        for (int pos = colptr[col]; pos < colptr[col + 1]; ++pos) {
          const int row = row_index[pos];
          const double value = values[pos];
          acc.sum[row] += value;
          acc.sumsq[row] += value * value;
          acc.count[row] += 1;
        }
      }
    }
    for (int t = 0; t < max_threads; ++t) {
      for (int row = 0; row < rows; ++row) {
        out.sum[row] += local[t].sum[row];
        out.sumsq[row] += local[t].sumsq[row];
        out.count[row] += local[t].count[row];
      }
    }
    return out;
  }
#endif
  for (int col = 0; col < columns; ++col) {
    for (int pos = colptr[col]; pos < colptr[col + 1]; ++pos) {
      const int row = row_index[pos];
      const double value = values[pos];
      out.sum[row] += value;
      out.sumsq[row] += value * value;
      out.count[row] += 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
List sparse_row_mean_var(IntegerVector p, IntegerVector i, NumericVector x,
                                int nrow, int ncol) {
  const int* pp = INTEGER(p);
  const int* ip = INTEGER(i);
  const double* xp = REAL(x);
  RowAccumulation rows = accumulate_rows(pp, ip, xp, nrow, ncol);

  NumericVector mu(nrow);
  NumericVector variance(nrow);
  IntegerVector nnz(nrow);
  const double n = static_cast<double>(ncol);
  const double denom = n - 1.0;

  for (int row = 0; row < nrow; ++row) {
    mu[row] = rows.sum[row] / n;
    variance[row] = (rows.sumsq[row] - n * mu[row] * mu[row]) / denom;
    if (variance[row] < 0.0) variance[row] = 0.0;
    nnz[row] = rows.count[row];
  }

  return List::create(Named("mean") = mu, Named("variance") = variance,
                      Named("nnz") = nnz);
}

// [[Rcpp::export]]
List sparse_row_mean_var_dgc_list(List mats, int nrow) {
  const int n_layers = mats.size();
  if (n_layers < 1) {
    stop("mats must contain at least one dgCMatrix");
  }

  std::vector<double> sum(nrow, 0.0);
  std::vector<double> sumsq(nrow, 0.0);
  std::vector<int> nnz(nrow, 0);
  int ncol_total = 0;

  for (int layer = 0; layer < n_layers; ++layer) {
    S4 mat = mats[layer];
    IntegerVector p = mat.slot("p");
    IntegerVector i = mat.slot("i");
    NumericVector x = mat.slot("x");
    IntegerVector dim = mat.slot("Dim");
    if (dim.size() != 2 || dim[0] != nrow) {
      stop("All matrices must have the requested row count");
    }
    const int ncol = dim[1];
    ncol_total += ncol;
    const int* pp = INTEGER(p);
    const int* ip = INTEGER(i);
    const double* xp = REAL(x);
    for (int col = 0; col < ncol; ++col) {
      for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
        const int row = ip[pos];
        const double value = xp[pos];
        sum[row] += value;
        sumsq[row] += value * value;
        nnz[row] += 1;
      }
    }
  }

  NumericVector mu(nrow);
  NumericVector variance(nrow);
  IntegerVector nnz_out(nrow);
  const double n = static_cast<double>(ncol_total);
  const double denom = n - 1.0;
  for (int row = 0; row < nrow; ++row) {
    mu[row] = sum[row] / n;
    variance[row] = (sumsq[row] - n * mu[row] * mu[row]) / denom;
    if (variance[row] < 0.0) variance[row] = 0.0;
    nnz_out[row] = nnz[row];
  }

  return List::create(Named("mean") = mu,
                      Named("variance") = variance,
                      Named("nnz") = nnz_out,
                      Named("ncol") = ncol_total);
}

static void prepare_standardization(NumericVector mu, NumericVector sd,
                                    std::vector<double>& inv_sd,
                                    std::vector<double>& zero_offset) {
  const int rows = mu.size();
  for (int row = 0; row < rows; ++row) {
    if (sd[row] == 0.0) {
      inv_sd[row] = 0.0;
      zero_offset[row] = 0.0;
    } else {
      inv_sd[row] = 1.0 / sd[row];
      zero_offset[row] = mu[row] * inv_sd[row];
    }
  }
}

// [[Rcpp::export]]
NumericVector sparse_row_var_std(IntegerVector p, IntegerVector i, NumericVector x,
                                         int nrow, int ncol,
                                         NumericVector mu, NumericVector sd,
                                         double vmax, IntegerVector nnzPerRow) {
  std::vector<double> inv_sd_vec(nrow);
  std::vector<double> mu_isd_vec(nrow);
  prepare_standardization(mu, sd, inv_sd_vec, mu_isd_vec);

  std::vector<double> sumSq(nrow, 0.0);
  const int* pp = INTEGER(p);
  const int* ip = INTEGER(i);
  const double* xp = REAL(x);

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  if (max_threads > 1 && ncol >= max_threads * 64) {
    std::vector<std::vector<double> > local(max_threads, std::vector<double>(nrow, 0.0));
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      std::vector<double>& acc = local[tid];
#pragma omp for schedule(static)
      for (int col = 0; col < ncol; ++col) {
        for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
          const int row = ip[pos];
          const double isd = inv_sd_vec[row];
          if (isd == 0.0) continue;
          double z = xp[pos] * isd - mu_isd_vec[row];
          if (z > vmax) z = vmax;
          acc[row] += z * z;
        }
      }
    }
    for (int t = 0; t < max_threads; ++t) {
      for (int row = 0; row < nrow; ++row) {
        sumSq[row] += local[t][row];
      }
    }
  } else {
#endif
  for (int col = 0; col < ncol; ++col) {
    for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
      const int row = ip[pos];
      const double isd = inv_sd_vec[row];
      if (isd == 0.0) continue;
      double z = xp[pos] * isd - mu_isd_vec[row];
      if (z > vmax) z = vmax;
      sumSq[row] += z * z;
    }
  }
#ifdef _OPENMP
  }
#endif

  NumericVector result(nrow);
  const double denom = ncol - 1.0;
  for (int row = 0; row < nrow; ++row) {
    if (inv_sd_vec[row] == 0.0) {
      result[row] = 0.0;
      continue;
    }
    const int nZero = ncol - nnzPerRow[row];
    const double zeroVal = mu_isd_vec[row];
    const double total = sumSq[row] + zeroVal * zeroVal * nZero;
    result[row] = total / denom;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector sparse_row_var_std_dgc_list(List mats, int nrow,
                                          NumericVector mu, NumericVector sd,
                                          double vmax) {
  std::vector<double> inv_sd_vec(nrow);
  std::vector<double> mu_isd_vec(nrow);
  prepare_standardization(mu, sd, inv_sd_vec, mu_isd_vec);

  std::vector<double> sumSq(nrow, 0.0);
  std::vector<int> nnz(nrow, 0);
  int ncol_total = 0;
  const int n_layers = mats.size();

  for (int layer = 0; layer < n_layers; ++layer) {
    S4 mat = mats[layer];
    IntegerVector p = mat.slot("p");
    IntegerVector i = mat.slot("i");
    NumericVector x = mat.slot("x");
    IntegerVector dim = mat.slot("Dim");
    if (dim.size() != 2 || dim[0] != nrow) {
      stop("All matrices must have the requested row count");
    }
    const int ncol = dim[1];
    ncol_total += ncol;
    const int* pp = INTEGER(p);
    const int* ip = INTEGER(i);
    const double* xp = REAL(x);
    for (int col = 0; col < ncol; ++col) {
      for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
        const int row = ip[pos];
        const double isd = inv_sd_vec[row];
        nnz[row] += 1;
        if (isd == 0.0) continue;
        double z = xp[pos] * isd - mu_isd_vec[row];
        if (z > vmax) z = vmax;
        sumSq[row] += z * z;
      }
    }
  }

  NumericVector result(nrow);
  const double denom = ncol_total - 1.0;
  for (int row = 0; row < nrow; ++row) {
    if (inv_sd_vec[row] == 0.0) {
      result[row] = 0.0;
      continue;
    }
    const int nZero = ncol_total - nnz[row];
    const double zeroVal = mu_isd_vec[row];
    const double total = sumSq[row] + zeroVal * zeroVal * nZero;
    result[row] = total / denom;
  }
  return result;
}
