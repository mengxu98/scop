#include <Rcpp.h>
#include <thisutils/log_message.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;



// [[Rcpp::export]]
List sparse_row_mean_var(IntegerVector p, IntegerVector i, NumericVector x,
                                int nrow, int ncol) {
  const int* pp = INTEGER(p);
  const int* ip = INTEGER(i);
  const double* xp = REAL(x);

  std::vector<double> sum(nrow, 0.0);
  const int nnz_total = x.size();
  for (int k = 0; k < nnz_total; ++k) {
    sum[ip[k]] += xp[k];
  }

  NumericVector mu(nrow);
  NumericVector variance(nrow);
  IntegerVector nnz_out(nrow);
  const double n = static_cast<double>(ncol);
  const double denom = n - 1.0;
  for (int row = 0; row < nrow; ++row) {
    mu[row] = sum[row] / n;
  }
  std::vector<double> rowvar(nrow, 0.0);
  std::vector<int> nzero(nrow, ncol);
  for (int k = 0; k < nnz_total; ++k) {
    const int row = ip[k];
    const double diff = xp[k] - mu[row];
    rowvar[row] += diff * diff;
    nzero[row] -= 1;
  }
  for (int row = 0; row < nrow; ++row) {
    variance[row] = (rowvar[row] + mu[row] * mu[row] * nzero[row]) / denom;
    if (variance[row] < 0.0) variance[row] = 0.0;
    nnz_out[row] = ncol - nzero[row];
  }

  return List::create(Named("mean") = mu, Named("variance") = variance,
                      Named("nnz") = nnz_out);
}

// [[Rcpp::export]]
List sparse_row_mean_var_dgc_list(List mats, int nrow) {
  const int n_layers = mats.size();
  if (n_layers < 1) {
    thisutils::log_message("mats must contain at least one dgCMatrix", "error");
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
      thisutils::log_message("All matrices must have the requested row count", "error");
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


// [[Rcpp::export]]
NumericVector sparse_row_var_std(IntegerVector p, IntegerVector i, NumericVector x,
                                         int nrow, int ncol,
                                         NumericVector mu, NumericVector sd,
                                         double vmax, IntegerVector nnzPerRow) {
  const double* mup = REAL(mu);
  const double* sdp = REAL(sd);
  const int* nnz = INTEGER(nnzPerRow);

  std::vector<double> sumSq(nrow, 0.0);
  const int* pp = INTEGER(p);
  const int* ip = INTEGER(i);
  const double* xp = REAL(x);

  for (int col = 0; col < ncol; ++col) {
    for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
      const int row = ip[pos];
      if (sdp[row] == 0.0) continue;
      double z = (xp[pos] - mup[row]) / sdp[row];
      if (z > vmax) z = vmax;
      sumSq[row] += z * z;
    }
  }

  NumericVector result(nrow);
  const double denom = ncol - 1.0;
  for (int row = 0; row < nrow; ++row) {
    if (sdp[row] == 0.0) {
      result[row] = 0.0;
      continue;
    }
    const int nZero = ncol - nnz[row];
    const double zeroVal = (0.0 - mup[row]) / sdp[row];
    const double total = sumSq[row] + zeroVal * zeroVal * nZero;
    result[row] = total / denom;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector sparse_row_var_std_dgc_list(List mats, int nrow,
                                          NumericVector mu, NumericVector sd,
                                          double vmax) {
  const double* mup = REAL(mu);
  const double* sdp = REAL(sd);

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
      thisutils::log_message("All matrices must have the requested row count", "error");
    }
    const int ncol = dim[1];
    ncol_total += ncol;
    const int* pp = INTEGER(p);
    const int* ip = INTEGER(i);
    const double* xp = REAL(x);
    for (int col = 0; col < ncol; ++col) {
      for (int pos = pp[col]; pos < pp[col + 1]; ++pos) {
        const int row = ip[pos];
        nnz[row] += 1;
        if (sdp[row] == 0.0) continue;
        double z = (xp[pos] - mup[row]) / sdp[row];
        if (z > vmax) z = vmax;
        sumSq[row] += z * z;
      }
    }
  }

  NumericVector result(nrow);
  const double denom = ncol_total - 1.0;
  for (int row = 0; row < nrow; ++row) {
    if (sdp[row] == 0.0) {
      result[row] = 0.0;
      continue;
    }
    const int nZero = ncol_total - nnz[row];
    const double zeroVal = (0.0 - mup[row]) / sdp[row];
    const double total = sumSq[row] + zeroVal * zeroVal * nZero;
    result[row] = total / denom;
  }
  return result;
}
