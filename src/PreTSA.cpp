#include <RcppArmadillo.h>
#include <thisutils/log_message.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List pretsa_fit_block_cpp(
    const NumericMatrix& expression,
    const List& bases,
    const List& inverses,
    const IntegerVector& knots) {
  const int features = expression.nrow();
  const int samples = expression.ncol();
  const int models = bases.size();
  if (features < 1 || samples < 2 || models < 1 || inverses.size() != models ||
      knots.size() != models) {
    thisutils::log_message("invalid PreTSA kernel dimensions", "error");
  }

  arma::mat expr(const_cast<double*>(expression.begin()), features, samples, false, true);
  arma::mat fitted(features, samples, arma::fill::zeros);
  arma::vec best_bic(features);
  best_bic.fill(std::numeric_limits<double>::infinity());
  IntegerVector selected_knots(features);
  IntegerVector df1(features);
  IntegerVector df2(features);

  for (int model = 0; model < models; ++model) {
    NumericMatrix basis_r = bases[model];
    NumericMatrix inverse_r = inverses[model];
    if (basis_r.nrow() != samples || inverse_r.nrow() != basis_r.ncol() ||
        inverse_r.ncol() != basis_r.ncol()) {
      thisutils::log_message("incompatible PreTSA basis dimensions", "error");
    }
    arma::mat basis(basis_r.begin(), basis_r.nrow(), basis_r.ncol(), false, true);
    arma::mat inverse(inverse_r.begin(), inverse_r.nrow(), inverse_r.ncol(), false, true);
    arma::mat projection = inverse * basis.t();
    arma::mat beta = projection * expr.t();
    arma::mat fitted_model = (basis * beta).t();
    arma::vec sse = arma::sum(arma::square(expr - fitted_model), 1);
    const double parameters = static_cast<double>(basis.n_cols + 1);
    const double log_samples = std::log(static_cast<double>(samples));

    for (int feature = 0; feature < features; ++feature) {
      const double mse = sse[feature] / static_cast<double>(samples);
      const double bic = static_cast<double>(samples) *
          (1.0 + std::log(6.283185307179586476925286766559) + std::log(mse)) +
          log_samples * parameters;
      if (model == 0 || bic < best_bic[feature]) {
        best_bic[feature] = bic;
        selected_knots[feature] = knots[model];
        df1[feature] = basis.n_cols - 1;
        df2[feature] = samples - basis.n_cols;
        fitted.row(feature) = fitted_model.row(feature);
      }
    }
  }

  return List::create(
      _["fit"] = wrap(fitted),
      _["knotnum"] = selected_knots,
      _["df1"] = df1,
      _["df2"] = df2,
      _["bic"] = wrap(best_bic));
}

// [[Rcpp::export]]
List pretsa_curve_summary_cpp(
    const NumericMatrix& fitted,
    const NumericMatrix& expression,
    const NumericVector& pseudotime) {
  const int features = fitted.nrow();
  const int samples = fitted.ncol();
  if (expression.nrow() != features || expression.ncol() != samples ||
      pseudotime.size() != samples || samples < 1) {
    thisutils::log_message("invalid PreTSA summary dimensions", "error");
  }
  NumericVector peak(features);
  NumericVector valley(features);
  IntegerVector expressed(features);
  std::vector<double> sorted(samples);
  std::vector<double> selected;
  selected.reserve(samples);

  for (int feature = 0; feature < features; ++feature) {
    double minimum = expression(feature, 0);
    int expressed_count = 0;
    for (int sample = 0; sample < samples; ++sample) {
      sorted[sample] = fitted(feature, sample);
      minimum = std::min(minimum, expression(feature, sample));
    }
    for (int sample = 0; sample < samples; ++sample) {
      if (expression(feature, sample) > minimum) {
        ++expressed_count;
      }
    }
    expressed[feature] = expressed_count;
    std::sort(sorted.begin(), sorted.end());

    const double high_position = 0.99 * static_cast<double>(samples - 1);
    const int high_lower = static_cast<int>(std::floor(high_position));
    const int high_upper = static_cast<int>(std::ceil(high_position));
    const double high_fraction = high_position - high_lower;
    const double high_quantile = sorted[high_lower] +
        high_fraction * (sorted[high_upper] - sorted[high_lower]);
    const double low_position = 0.01 * static_cast<double>(samples - 1);
    const int low_lower = static_cast<int>(std::floor(low_position));
    const int low_upper = static_cast<int>(std::ceil(low_position));
    const double low_fraction = low_position - low_lower;
    const double low_quantile = sorted[low_lower] +
        low_fraction * (sorted[low_upper] - sorted[low_lower]);

    selected.clear();
    for (int sample = 0; sample < samples; ++sample) {
      if (fitted(feature, sample) >= high_quantile) {
        selected.push_back(pseudotime[sample]);
      }
    }
    const int high_count = selected.size();
    peak[feature] = high_count % 2 == 1
        ? selected[high_count / 2]
        : (selected[high_count / 2 - 1] + selected[high_count / 2]) / 2.0;

    selected.clear();
    for (int sample = 0; sample < samples; ++sample) {
      if (fitted(feature, sample) <= low_quantile) {
        selected.push_back(pseudotime[sample]);
      }
    }
    const int low_count = selected.size();
    valley[feature] = low_count % 2 == 1
        ? selected[low_count / 2]
        : (selected[low_count / 2 - 1] + selected[low_count / 2]) / 2.0;
  }

  return List::create(
      _["peaktime"] = peak,
      _["valleytime"] = valley,
      _["exp_ncells"] = expressed);
}
