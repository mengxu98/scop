#include <Rcpp.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

struct CsrMatrixView {
  const int* row_start;
  const int* column;
  const double* value;
};

static double clamp_value(double value, double lower, double upper) {
  if (value < lower) return lower;
  if (value > upper) return upper;
  return value;
}

static double nb_variance(double mean, double theta, double lower_bound) {
  double variance = mean + mean * mean / theta;
  return variance < lower_bound ? lower_bound : variance;
}

static double pearson_residual(double observed, double mean,
                               double theta, double min_var) {
  return (observed - mean) / std::sqrt(nb_variance(mean, theta, min_var));
}

// [[Rcpp::export]]
List csc_to_csr(IntegerVector csc_i,
                IntegerVector csc_p,
                NumericVector csc_x,
                int nrow,
                int ncol) {
  const int nnz = csc_x.size();
  IntegerVector row_ptr(nrow + 1);
  for (int pos = 0; pos < nnz; ++pos) {
    ++row_ptr[csc_i[pos] + 1];
  }
  for (int row = 0; row < nrow; ++row) {
    row_ptr[row + 1] += row_ptr[row];
  }

  IntegerVector columns(nnz);
  NumericVector values(nnz);
  std::vector<int> cursor(row_ptr.begin(), row_ptr.end());
  for (int col = 0; col < ncol; ++col) {
    for (int pos = csc_p[col]; pos < csc_p[col + 1]; ++pos) {
      const int row = csc_i[pos];
      const int target = cursor[row]++;
      columns[target] = col;
      values[target] = csc_x[pos];
    }
  }

  return List::create(
    Named("row_ptr") = row_ptr,
    Named("col_idx") = columns,
    Named("vals") = values
  );
}

struct CorrectedRow {
  std::vector<int> columns;
  std::vector<double> values;
};

static int csc_size(const std::vector<CorrectedRow>& rows) {
  int total = 0;
  for (std::vector<CorrectedRow>::const_iterator it = rows.begin();
       it != rows.end(); ++it) {
    total += static_cast<int>(it->columns.size());
  }
  return total;
}

static List corrected_rows_to_csc(const std::vector<CorrectedRow>& rows,
                                  int cells) {
  std::vector<int> counts(cells, 0);
  for (size_t row = 0; row < rows.size(); ++row) {
    for (size_t k = 0; k < rows[row].columns.size(); ++k) {
      ++counts[rows[row].columns[k]];
    }
  }

  IntegerVector p(cells + 1);
  for (int col = 0; col < cells; ++col) {
    p[col + 1] = p[col] + counts[col];
  }

  const int total = csc_size(rows);
  IntegerVector i(total);
  NumericVector x(total);
  std::vector<int> cursor(p.begin(), p.end());
  for (int row = 0; row < static_cast<int>(rows.size()); ++row) {
    for (size_t k = 0; k < rows[row].columns.size(); ++k) {
      const int col = rows[row].columns[k];
      const int target = cursor[col]++;
      i[target] = row;
      x[target] = rows[row].values[k];
    }
  }

  return List::create(Named("csc_i") = i, Named("csc_p") = p,
                      Named("csc_x") = x);
}

// [[Rcpp::export]]
List sct_stats_correct_sparse(NumericVector intercepts,
                              NumericVector cell_mu_base,
                              IntegerVector csr_row_ptr,
                              IntegerVector csr_col_idx,
                              NumericVector csr_vals,
                              IntegerVector gene_idx,
                              NumericVector theta,
                              NumericVector corr_factor,
                              double min_var,
                              double clip_lo,
                              double clip_hi,
                              bool do_correct) {
  const int genes = gene_idx.size();
  const int cells = cell_mu_base.size();
  NumericVector residual_var(genes);
  NumericVector residual_mean(genes);
  std::vector<CorrectedRow> corrected(do_correct ? genes : 0);
  CsrMatrixView csr{INTEGER(csr_row_ptr), INTEGER(csr_col_idx), REAL(csr_vals)};
  const int* source_gene = INTEGER(gene_idx);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 50)
  #endif
  for (int out_row = 0; out_row < genes; ++out_row) {
    const int input_row = source_gene[out_row];
    const double theta_g = theta[out_row];
    const double intercept = std::exp(intercepts[out_row]);
    int cursor = csr.row_start[input_row];
    const int end = csr.row_start[input_row + 1];
    double sum = 0.0;
    double sumsq = 0.0;
    CorrectedRow local;

    for (int cell = 0; cell < cells; ++cell) {
      double observed = 0.0;
      if (cursor < end && csr.column[cursor] == cell) {
        observed = csr.value[cursor++];
      }
      const double mean = intercept * cell_mu_base[cell];
      double residual = pearson_residual(observed, mean, theta_g, min_var);
      residual = clamp_value(residual, clip_lo, clip_hi);
      sum += residual;
      sumsq += residual * residual;

      if (do_correct) {
        const double corrected_mean = mean * corr_factor[cell];
        const double corrected_var = nb_variance(corrected_mean, theta_g, min_var);
        const double count = std::round(corrected_mean +
          residual * std::sqrt(corrected_var));
        if (count > 0.0) {
          local.columns.push_back(cell);
          local.values.push_back(count);
        }
      }
    }

    const double mean_residual = sum / cells;
    residual_mean[out_row] = mean_residual;
    residual_var[out_row] =
      (sumsq / cells - mean_residual * mean_residual) *
      cells / (cells - 1.0);
    if (do_correct) {
      corrected[out_row].columns.swap(local.columns);
      corrected[out_row].values.swap(local.values);
    }
  }

  if (!do_correct) {
    return List::create(Named("res_var") = residual_var,
                        Named("res_mean") = residual_mean);
  }

  List csc = corrected_rows_to_csc(corrected, cells);
  csc["res_var"] = residual_var;
  csc["res_mean"] = residual_mean;
  return csc;
}

// [[Rcpp::export]]
NumericMatrix sct_fused_resid_center_sparse(
    NumericVector intercepts,
    NumericVector cell_mu_base,
    IntegerVector csr_row_ptr,
    IntegerVector csr_col_idx,
    NumericVector csr_vals,
    IntegerVector gene_idx,
    NumericVector theta,
    double min_var,
    double wide_clip_lo,
    double wide_clip_hi,
    double narrow_clip_lo,
    double narrow_clip_hi) {
  const int genes = gene_idx.size();
  const int cells = cell_mu_base.size();
  NumericMatrix output(genes, cells);
  CsrMatrixView csr{INTEGER(csr_row_ptr), INTEGER(csr_col_idx), REAL(csr_vals)};
  const int* source_gene = INTEGER(gene_idx);
  const bool skip_wide_clip =
    narrow_clip_lo >= wide_clip_lo && narrow_clip_hi <= wide_clip_hi;

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 50)
  #endif
  for (int out_row = 0; out_row < genes; ++out_row) {
    const int input_row = source_gene[out_row];
    const double theta_g = theta[out_row];
    const double intercept = std::exp(intercepts[out_row]);
    int cursor = csr.row_start[input_row];
    const int end = csr.row_start[input_row + 1];
    double total = 0.0;

    for (int cell = 0; cell < cells; ++cell) {
      double observed = 0.0;
      if (cursor < end && csr.column[cursor] == cell) {
        observed = csr.value[cursor++];
      }
      const double mean = intercept * cell_mu_base[cell];
      double residual = pearson_residual(observed, mean, theta_g, min_var);
      if (!skip_wide_clip) {
        residual = clamp_value(residual, wide_clip_lo, wide_clip_hi);
      }
      residual = clamp_value(residual, narrow_clip_lo, narrow_clip_hi);
      output(out_row, cell) = residual;
      total += residual;
    }

    const double center = total / cells;
    for (int cell = 0; cell < cells; ++cell) {
      output(out_row, cell) -= center;
    }
  }
  return output;
}
