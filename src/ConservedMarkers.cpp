#include <Rcpp.h>
#include <thisutils/log_message.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector combine_conserved_pvalues_cpp(
    const NumericMatrix& pvalues,
    const std::string& method) {
  const int rows = pvalues.nrow();
  const int cols = pvalues.ncol();
  NumericVector result(rows, NA_REAL);

  for (int row = 0; row < rows; ++row) {
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
    } else if (method == "sump") {
      long double total_extended = 0.0;
      for (double value : valid) {
        total_extended += value;
      }
      const double total = static_cast<double>(total_extended);
      const int terms = static_cast<int>(std::floor(total)) + 1;
      long double combined = 0.0;
      for (int index = 1; index <= terms; ++index) {
        const double base = total - index + 1.0;
        double term = 0.0;
        if (base > 0.0) {
          term = std::exp(
              R::lchoose(k, index - 1) + k * std::log(base) - R::lgammafn(k + 1));
        }
        combined += index % 2 == 1 ? term : -term;
      }
      result[row] = static_cast<double>(combined);
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
