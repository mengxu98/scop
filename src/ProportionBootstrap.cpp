// [[Rcpp::depends(RcppArmadillo, cli)]]
#include <RcppArmadillo.h>
#include <thisutils/cli_progress.h>
#include <cmath>
#include <random>
#include <vector>

using namespace Rcpp;

namespace {

// R's type-7 quantile algorithm on a mutable vector
double quantile_type7(std::vector<double>& x, double prob) {
  int n = static_cast<int>(x.size());
  if (n == 0) return NA_REAL;

  // Count non-NA values
  std::vector<double> valid;
  valid.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (!R_IsNA(x[i]) && !R_IsNaN(x[i]) && std::isfinite(x[i])) {
      valid.push_back(x[i]);
    }
  }

  int n_valid = static_cast<int>(valid.size());
  if (n_valid == 0) return NA_REAL;
  if (n_valid == 1) return valid[0];

  std::sort(valid.begin(), valid.end());

  double index = prob * (n_valid - 1);
  int lo = static_cast<int>(std::floor(index));
  int hi = static_cast<int>(std::ceil(index));
  double gamma = index - lo;

  if (lo < 0) lo = 0;
  if (hi >= n_valid) hi = n_valid - 1;

  return (1.0 - gamma) * valid[lo] + gamma * valid[hi];
}

// Mean ignoring NA/NaN/Inf
double mean_valid(const std::vector<double>& x) {
  double sum = 0.0;
  int n = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (!R_IsNA(x[i]) && !R_IsNaN(x[i]) && std::isfinite(x[i])) {
      sum += x[i];
      ++n;
    }
  }
  return (n > 0) ? sum / n : NA_REAL;
}

// Bootstrap resampling: draw n samples with replacement, compute mean
double bootstrap_sample_mean(
    const std::vector<double>& v,
    int n,
    std::mt19937& rng
) {
  std::uniform_int_distribution<int> dist(0, static_cast<int>(v.size()) - 1);
  double sum = 0.0;
  int count = 0;
  for (int i = 0; i < n; ++i) {
    int idx = dist(rng);
    double val = v[idx];
    if (!R_IsNA(val) && !R_IsNaN(val) && std::isfinite(val)) {
      sum += val;
      ++count;
    }
  }
  return (count > 0) ? sum / count : NA_REAL;
}

}  // namespace

// [[Rcpp::export]]
NumericVector proportion_bootstrap_log2fd(
    NumericVector v1,
    NumericVector v2,
    int n_bootstrap = 1000,
    double pseudocount = 1e-5,
    bool verbose = false
) {
  int n1 = v1.size();
  int n2 = v2.size();

  if (n_bootstrap <= 0 || n1 == 0 || n2 == 0) {
    NumericVector result(1);
    result[0] = NA_REAL;
    return result;
  }

  // Copy to std::vector for faster random access
  std::vector<double> vec1(n1), vec2(n2);
  for (int i = 0; i < n1; ++i) vec1[i] = v1[i];
  for (int i = 0; i < n2; ++i) vec2[i] = v2[i];

  // Use Mersenne Twister with a seed from R's RNG for reproducibility
  // Sample a seed from R's current RNG state
  int seed = static_cast<int>(std::floor(R::runif(0.0, 1.0) * 2147483647.0));
  std::mt19937 rng(seed);

  NumericVector boot(n_bootstrap);
  thisutils::cli_progress progress(
    n_bootstrap,
    verbose,
    "Bootstrap proportion log2 fold differences"
  );

  for (int b = 0; b < n_bootstrap; ++b) {
    if (thisutils::should_check_interrupt(b, n_bootstrap, verbose)) {
      Rcpp::checkUserInterrupt();
    }
    progress.set(b);
    double mean1 = bootstrap_sample_mean(vec1, n1, rng);
    double mean2 = bootstrap_sample_mean(vec2, n2, rng);

    if (R_IsNA(mean1) || R_IsNaN(mean1) || !std::isfinite(mean1) ||
        R_IsNA(mean2) || R_IsNaN(mean2) || !std::isfinite(mean2)) {
      boot[b] = NA_REAL;
    } else {
      boot[b] = std::log2((mean2 + pseudocount) / (mean1 + pseudocount));
    }
  }
  progress.set(n_bootstrap, true);

  return boot;
}

// [[Rcpp::export]]
List proportion_bootstrap_stats(
    NumericVector v1,
    NumericVector v2,
    int n_bootstrap = 1000,
    double pseudocount = 1e-5,
    bool verbose = false
) {
  NumericVector boot = proportion_bootstrap_log2fd(v1, v2, n_bootstrap, pseudocount, verbose);

  // Copy to std::vector for quantile computation
  std::vector<double> boot_vec(boot.size());
  for (int i = 0; i < boot.size(); ++i) boot_vec[i] = boot[i];

  double boot_mean = mean_valid(boot_vec);
  double boot_ci_low = quantile_type7(boot_vec, 0.025);
  double boot_ci_high = quantile_type7(boot_vec, 0.975);

  return List::create(
    _["boot_mean_log2FD"] = boot_mean,
    _["boot_CI_2.5"] = boot_ci_low,
    _["boot_CI_97.5"] = boot_ci_high,
    _["boot_samples"] = boot
  );
}
