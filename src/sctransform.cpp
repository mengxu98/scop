#include <Rcpp.h>
#include <algorithm>
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

  for (int out_row = 0; out_row < genes; ++out_row) {
    const int input_row = source_gene[out_row];
    const double theta_g = theta[out_row];
    const double intercept = std::exp(intercepts[out_row]);
    int cursor = csr.row_start[input_row];
    const int end = csr.row_start[input_row + 1];
    double sum = 0.0;
    CorrectedRow local;

    for (int cell = 0; cell < cells; ++cell) {
      double observed = 0.0;
      if (cursor < end && csr.column[cursor] == cell) {
        observed = csr.value[cursor++];
      }
      const double mean = intercept * cell_mu_base[cell];
      const double raw_residual = pearson_residual(
        observed, mean, theta_g, min_var
      );
      const double residual = clamp_value(raw_residual, clip_lo, clip_hi);
      sum += residual;

      if (do_correct) {
        const double corrected_mean = mean * corr_factor[cell];
        const double corrected_var = nb_variance(
          corrected_mean, theta_g, 0.0
        );
        const double count = std::nearbyint(corrected_mean +
          raw_residual * std::sqrt(corrected_var));
        if (count > 0.0) {
          local.columns.push_back(cell);
          local.values.push_back(count);
        }
      }
    }

    const double mean_residual = sum / cells;
    cursor = csr.row_start[input_row];
    double squared_deviation = 0.0;
    for (int cell = 0; cell < cells; ++cell) {
      double observed = 0.0;
      if (cursor < end && csr.column[cursor] == cell) {
        observed = csr.value[cursor++];
      }
      const double mean = intercept * cell_mu_base[cell];
      const double residual = clamp_value(
        pearson_residual(observed, mean, theta_g, min_var),
        clip_lo,
        clip_hi
      );
      const double delta = residual - mean_residual;
      squared_deviation += delta * delta;
    }
    residual_mean[out_row] = mean_residual;
    residual_var[out_row] = squared_deviation / (cells - 1.0);
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
    double narrow_clip_hi,
    bool do_center,
    bool do_scale) {
  const int genes = gene_idx.size();
  const int cells = cell_mu_base.size();
  NumericMatrix output(genes, cells);
  CsrMatrixView csr{INTEGER(csr_row_ptr), INTEGER(csr_col_idx), REAL(csr_vals)};
  const int* source_gene = INTEGER(gene_idx);
  const bool skip_wide_clip =
    narrow_clip_lo >= wide_clip_lo && narrow_clip_hi <= wide_clip_hi;

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

    const double center = do_center ? total / cells : 0.0;
    double sumsq = 0.0;
    for (int cell = 0; cell < cells; ++cell) {
      output(out_row, cell) -= center;
      sumsq += output(out_row, cell) * output(out_row, cell);
    }
    if (do_scale) {
      double scale = std::sqrt(sumsq / std::max(cells - 1, 1));
      if (!std::isfinite(scale) || scale == 0.0) scale = 1.0;
      for (int cell = 0; cell < cells; ++cell) {
        output(out_row, cell) /= scale;
      }
    }
  }
  return output;
}

struct SctDispersionStats {
  std::vector<double> values;
  std::vector<double> frequencies;
  double sum;
  double sumsq;
  double max_value;
  bool all_zero;
};

static SctDispersionStats sct_dispersion_stats(const double* y, int n) {
  SctDispersionStats stats;
  stats.sum = 0.0;
  stats.sumsq = 0.0;
  stats.max_value = 0.0;
  stats.all_zero = true;
  bool nonnegative_integer = true;
  for (int i = 0; i < n; ++i) {
    const double value = y[i];
    stats.sum += value;
    stats.sumsq += value * value;
    stats.max_value = std::max(stats.max_value, value);
    stats.all_zero = stats.all_zero && value == 0.0;
    nonnegative_integer = nonnegative_integer && value >= 0.0 &&
      value == std::floor(value);
  }
  const bool use_table = nonnegative_integer && stats.max_value <= 1000000.0 &&
    stats.max_value <= std::max(256.0, 4.0 * n);
  if (use_table) {
    std::vector<int> table(static_cast<size_t>(stats.max_value) + 1, 0);
    for (int i = 0; i < n; ++i) {
      ++table[static_cast<size_t>(y[i])];
    }
    for (size_t value = 0; value < table.size(); ++value) {
      if (table[value] > 0) {
        stats.values.push_back(static_cast<double>(value));
        stats.frequencies.push_back(static_cast<double>(table[value]));
      }
    }
  } else {
    std::vector<double> sorted(y, y + n);
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < n; ++i) {
      const double value = sorted[i];
      if (stats.values.empty() || value != stats.values.back()) {
        stats.values.push_back(value);
        stats.frequencies.push_back(1.0);
      } else {
        stats.frequencies.back() += 1.0;
      }
    }
  }
  return stats;
}

// Intercept-only specialization of glmGamPoi's fitBeta_one_group(). It keeps
// the same Newton update and convergence rule, but uses a shared cell offset
// vector and omits the unused per-gene deviance calculation.
// [[Rcpp::export]]
List sct_fit_beta_intercept_offset(
    NumericMatrix y,
    NumericVector log_offset,
    NumericVector dispersions,
    Nullable<NumericVector> beta_start = R_NilValue,
    bool return_mu = false,
    double tolerance = 1e-8,
    int max_iter = 100) {
  const int genes = y.nrow();
  const int cells = y.ncol();
  if (log_offset.size() != cells || dispersions.size() != genes) {
    stop("offset and dispersion dimensions do not match y");
  }
  NumericVector starts;
  const bool supplied_start = beta_start.isNotNull();
  if (supplied_start) {
    starts = NumericVector(beta_start);
    if (starts.size() != genes) stop("beta_start length does not match y");
  }
  NumericVector beta(genes);
  LogicalVector converged(genes, true);
  NumericMatrix mu(return_mu ? genes : 0, return_mu ? cells : 0);

  for (int gene = 0; gene < genes; ++gene) {
    if (gene % 100 == 0) checkUserInterrupt();
    double current;
    if (supplied_start) {
      current = starts[gene];
    } else {
      double normalized_sum = 0.0;
      for (int cell = 0; cell < cells; ++cell) {
        normalized_sum += y(gene, cell) / std::exp(log_offset[cell]);
      }
      current = std::log(normalized_sum / cells);
    }
    const double dispersion = dispersions[gene];
    if (NumericVector::is_na(current) || NumericVector::is_na(dispersion)) {
      beta[gene] = NA_REAL;
      converged[gene] = false;
      continue;
    }
    bool all_zero = true;
    int iter = 0;
    for (; iter < max_iter; ++iter) {
      double dl = 0.0;
      double ddl = 0.0;
      all_zero = true;
      for (int cell = 0; cell < cells; ++cell) {
        const double count = y(gene, cell);
        all_zero = all_zero && count == 0.0;
        const double mean = std::exp(current + log_offset[cell]);
        const double denom = 1.0 + mean * dispersion;
        dl += (count - mean) / denom;
        ddl += mean * (1.0 + count * dispersion) / (denom * denom);
      }
      if (all_zero) {
        current = R_NegInf;
        break;
      }
      const double step = dl / ddl;
      current += step;
      if (std::abs(step) < tolerance) break;
      if (std::isnan(current)) break;
    }
    if (iter == max_iter || std::isnan(current)) {
      converged[gene] = false;
    }
    beta[gene] = std::max(current, -1e8);
    if (return_mu && converged[gene]) {
      for (int cell = 0; cell < cells; ++cell) {
        mu(gene, cell) = std::exp(current + log_offset[cell]);
      }
    }
  }
  return List::create(
    Named("beta") = beta,
    Named("converged") = converged,
    Named("mu") = mu
  );
}

static double sct_dispersion_score(
    const double* y,
    const double* mu,
    int n,
    double log_theta,
    const SctDispersionStats& stats) {
  const double theta = std::exp(log_theta);
  const double theta_inv = 1.0 / theta;
  double sum_w = 0.0;
  double sum_dw = 0.0;
  double raw_digamma = 0.0;
  double sum_prod_y = 0.0;
  for (size_t i = 0; i < stats.values.size(); ++i) {
    const double count = stats.values[i];
    const double frequency = stats.frequencies[i];
    raw_digamma += frequency * R::digamma(count + theta_inv);
    sum_prod_y += frequency * (count - 1.0) * count;
  }
  double score_digamma = raw_digamma;
  const double correction = theta_inv > 1e5 ?
    sum_prod_y / (2.0 * theta_inv) : 0.0;
  if (stats.max_value * 1e6 < theta_inv) {
    score_digamma = stats.sum - correction;
  } else {
    score_digamma -= n * R::digamma(theta_inv);
    score_digamma *= theta_inv;
    score_digamma = std::min(score_digamma, stats.sum - correction);
  }
  double likelihood_part = 0.0;
  for (int i = 0; i < n; ++i) {
    const double mean = mu[i] == 0.0 ? 1e-6 : mu[i];
    const double mean_theta = mean * theta;
    if (mean_theta < 1e-10) {
      likelihood_part += mean_theta * mean_theta *
        (1.0 / (1.0 + mean_theta) - 0.5);
    } else if (mean_theta < 1e-4) {
      const double inverse = 1.0 / (1.0 + mean_theta);
      const double upper = mean_theta * mean_theta * inverse;
      const double lower = mean_theta * mean_theta * (inverse - 0.5);
      const double suggested = std::log(1.0 + mean_theta) -
        mean / (mean + theta_inv);
      likelihood_part += std::max(std::min(suggested, upper), lower);
    } else {
      likelihood_part += std::log(1.0 + mean_theta) -
        mean / (mean + theta_inv);
    }
    likelihood_part += y[i] / (mean + theta_inv);
    const double weight = mean / (1.0 + mean_theta);
    sum_w += weight;
    sum_dw -= weight * weight;
  }
  likelihood_part *= theta_inv;
  const double cr_score = -0.5 * sum_dw /
    (sum_w + 1e-6) * 0.99 * theta;
  return likelihood_part - score_digamma + cr_score;
}

// Specialized single-core dispersion MLE for the intercept-only design used
// by SCTransform v2. This removes glmGamPoi's per-gene callback into R while
// preserving its score function and Cox-Reid adjustment.
// [[Rcpp::export]]
NumericVector sct_overdispersion_mle_intercept(
    NumericMatrix y,
    NumericMatrix mu,
    int max_iter = 200,
    double tolerance = 1e-8) {
  if (y.nrow() != mu.nrow() || y.ncol() != mu.ncol()) {
    stop("y and mu must have identical dimensions");
  }
  const int genes = y.nrow();
  const int cells = y.ncol();
  NumericVector estimates(genes);
  const double lower = std::log(1e-16);
  const double upper = std::log(1e16);
  const double far_left = std::log(1e-8);
  std::vector<double> y_row(cells);
  std::vector<double> mu_row(cells);

  for (int gene = 0; gene < genes; ++gene) {
    if (gene % 100 == 0) checkUserInterrupt();
    for (int cell = 0; cell < cells; ++cell) {
      y_row[cell] = y(gene, cell);
      mu_row[cell] = mu(gene, cell);
    }
    const SctDispersionStats stats = sct_dispersion_stats(y_row.data(), cells);
    if (stats.all_zero) {
      estimates[gene] = 0.0;
      continue;
    }
    const double left_score = sct_dispersion_score(
      y_row.data(), mu_row.data(), cells, far_left, stats
    );
    if (left_score < 0.0) {
      estimates[gene] = 0.0;
      continue;
    }
    double mean = stats.sum / cells;
    double variance = cells > 1 ?
      (stats.sumsq - cells * mean * mean) / (cells - 1.0) : 0.0;
    double start = (variance - mean) / (mean * mean);
    if (!std::isfinite(start) || start <= 0.0) start = 0.5;
    double lo = lower;
    double hi = upper;
    double current = std::max(lo, std::min(hi, std::log(start)));
    for (int iter = 0; iter < max_iter; ++iter) {
      const double score = sct_dispersion_score(
        y_row.data(), mu_row.data(), cells, current, stats
      );
      if (score > 0.0) lo = current; else hi = current;
      const double next = 0.5 * (lo + hi);
      if (std::abs(next - current) < tolerance) {
        current = next;
        break;
      }
      current = next;
    }
    estimates[gene] = std::exp(current);
  }
  return estimates;
}
