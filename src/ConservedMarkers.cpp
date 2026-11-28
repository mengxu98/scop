#include <Rcpp.h>
#include <thisutils/log_message.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

// CDF of the sum of k independent Uniform(0, 1) variables. The alternating
// polynomial loses precision near the center and especially in the upper
// half; use symmetry and a recurrence with non-negative weights instead.
double uniform_sum_cdf(int k, double x, std::vector<double>& work) {
  if (x <= 0.0) return 0.0;
  if (x >= k) return 1.0;
  if (x == 0.5 * k) return 0.5;

  const bool reflect = x > 0.5 * k;
  if (reflect) x = k - x;
  if (x <= 1.0) {
    const double lower = std::exp(k * std::log(x) - R::lgammafn(k + 1.0));
    return reflect ? 1.0 - lower : lower;
  }

  const int last = static_cast<int>(x);
  // At each order, work[i] stores F_order(x - i). The starting order is 1.
  work.assign(last + 1, 1.0);
  work[last] = x - last;
  for (int order = 2; order <= k; ++order) {
    if (k > 1024 && order % 128 == 0) Rcpp::checkUserInterrupt();
    const double inverse_order = 1.0 / order;
    const int first = std::max(0, last - order + 1);
    // Scores outside the support stay at 1 (left) or are 0 (right).
    // Increasing i keeps work[i + 1] from the preceding order.
    for (int i = first; i <= last; ++i) {
      const double weight = (x - i) * inverse_order;
      const double right = i < last ? work[i + 1] : 0.0;
      work[i] = weight * work[i] + (1.0 - weight) * right;
    }
  }
  return reflect ? 1.0 - work[0] : work[0];
}

}

// [[Rcpp::export]]
NumericVector combine_conserved_pvalues_cpp(
    const NumericMatrix& pvalues,
    const std::string& method) {
  const int rows = pvalues.nrow();
  const int cols = pvalues.ncol();
  NumericVector result(rows, NA_REAL);
  std::vector<double> sump_work;
  const bool sum_method = method == "sump";

  for (int row = 0; row < rows; ++row) {
    if (sum_method && row % 64 == 0) Rcpp::checkUserInterrupt();
    std::vector<double> valid;
    valid.reserve(cols);
    for (int col = 0; col < cols; ++col) {
      const double value = pvalues(row, col);
      if (R_finite(value) && value >= 0.0 && value <= 1.0) {
        valid.push_back(value);
      }
    }
    const int k = valid.size();
    if ((method == "meanp" && k < 4) ||
        (method != "meanp" && k < 2)) {
      continue;
    }

    if (method == "maximump") {
      const double value = *std::max_element(valid.begin(), valid.end());
      result[row] = R::pbeta(value, k, 1.0, true, false);
    } else if (method == "minimump" || method == "wilkinsonp") {
      const double value = *std::min_element(valid.begin(), valid.end());
      result[row] = R::pbeta(value, 1.0, k, true, false);
    } else if (method == "meanp") {
      long double total = 0.0;
      for (double value : valid) {
        total += value;
      }
      const double mean = static_cast<double>(total / k);
      const double z = (0.5 - mean) * std::sqrt(12.0 * k);
      result[row] = R::pnorm(z, 0.0, 1.0, false, false);
    } else if (sum_method) {
      long double total = 0.0;
      for (double value : valid) {
        total += value;
      }
      result[row] = uniform_sum_cdf(k, static_cast<double>(total), sump_work);
    } else if (method == "votep") {
      int positive = 0;
      int negative = 0;
      for (double value : valid) {
        positive += value < 0.5;
        negative += value > 0.5;
      }
      const int votes = positive + negative;
      result[row] = votes == 0 ? 1.0 : R::pbinom(positive - 1, votes, 0.5, false, false);
    } else {
      thisutils::log_message("unsupported conserved-marker p-value method", "error");
    }
  }
  return result;
}
