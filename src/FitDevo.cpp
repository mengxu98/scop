#include <Rcpp.h>
#include <thisutils/log_message.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector fitdevo_spearman_weights_cpp(
    Rcpp::NumericMatrix scaled,
    Rcpp::NumericVector target_centered) {
  const int n_features = scaled.nrow();
  const int n_cells = scaled.ncol();
  if (target_centered.size() != n_cells) {
    thisutils::log_message("Target rank length must match the number of cells.", "error");
  }
  double target_ss = 0.0;
  for (int cell = 0; cell < n_cells; ++cell) {
    if (!std::isfinite(target_centered[cell])) {
      thisutils::log_message("Target ranks must be finite.", "error");
    }
    target_ss += target_centered[cell] * target_centered[cell];
  }

  Rcpp::NumericVector output(n_features, 0.0);
  double* output_ptr = output.begin();
  const double* matrix_ptr = scaled.begin();
  const double* target_ptr = target_centered.begin();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector<std::pair<double, int> > ordered(n_cells);
    std::vector<double> ranks(n_cells);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int feature = 0; feature < n_features; ++feature) {
      bool finite = true;
      for (int cell = 0; cell < n_cells; ++cell) {
        const double value = matrix_ptr[feature + n_features * cell];
        if (!std::isfinite(value)) {
          finite = false;
          break;
        }
        ordered[cell] = std::make_pair(value, cell);
      }
      if (!finite || target_ss <= 0.0) {
        output_ptr[feature] = 0.0;
        continue;
      }
      std::sort(
        ordered.begin(),
        ordered.end(),
        [](const std::pair<double, int>& left,
           const std::pair<double, int>& right) {
          if (left.first < right.first) {
            return true;
          }
          if (left.first > right.first) {
            return false;
          }
          return left.second < right.second;
        }
      );

      int start = 0;
      while (start < n_cells) {
        int end = start + 1;
        while (end < n_cells && ordered[end].first == ordered[start].first) {
          ++end;
        }
        const double average_rank =
          (static_cast<double>(start + 1) + static_cast<double>(end)) / 2.0;
        for (int position = start; position < end; ++position) {
          ranks[ordered[position].second] = average_rank;
        }
        start = end;
      }

      double rank_sum = 0.0;
      double rank_sumsq = 0.0;
      double numerator = 0.0;
      for (int cell = 0; cell < n_cells; ++cell) {
        const double rank = ranks[cell];
        rank_sum += rank;
        rank_sumsq += rank * rank;
        numerator += rank * target_ptr[cell];
      }
      const double rank_ss = rank_sumsq -
        rank_sum * rank_sum / static_cast<double>(n_cells);
      const double denominator = std::sqrt(rank_ss * target_ss);
      output_ptr[feature] = denominator > 0.0 ?
        numerator / denominator : 0.0;
    }
  }
  return output;
}
