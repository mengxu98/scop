#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

namespace {

double gniplr_cor(const arma::rowvec& x, const arma::rowvec& y) {
  const double mx = arma::mean(x);
  const double my = arma::mean(y);
  const arma::rowvec xc = x - mx;
  const arma::rowvec yc = y - my;
  const double sx = std::sqrt(arma::dot(xc, xc));
  const double sy = std::sqrt(arma::dot(yc, yc));
  const double denom = sx * sy;
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  const double out = arma::dot(xc, yc) / denom;
  return std::isfinite(out) ? out : 0.0;
}

arma::mat gniplr_poly_raw(const arma::vec& x, int degree) {
  degree = std::max(1, degree);
  arma::mat out(x.n_elem, degree);
  for (int d = 1; d <= degree; ++d) {
    out.col(d - 1) = arma::pow(x, d);
  }
  return out;
}

double gniplr_soft_threshold(double x, double lambda) {
  if (x > lambda) return x - lambda;
  if (x < -lambda) return x + lambda;
  return 0.0;
}

arma::vec gniplr_project(
  const arma::rowvec& regulator,
  const arma::rowvec& target,
  bool predict_target,
  int degree,
  double alpha
) {
  const arma::uword n = regulator.n_elem;
  std::vector<arma::uword> order(n);
  for (arma::uword i = 0; i < n; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](arma::uword a, arma::uword b) {
    if (regulator[a] == regulator[b]) return target[a] < target[b];
    return regulator[a] < regulator[b];
  });
  arma::vec x(n), y(n), t(n);
  const double min_x = regulator[order.front()];
  const double max_x = regulator[order.back()];
  for (arma::uword i = 0; i < n; ++i) {
    x[i] = regulator[order[i]];
    y[i] = target[order[i]];
    t[i] = n <= 1 ? min_x : min_x + (max_x - min_x) * static_cast<double>(i) /
      static_cast<double>(n - 1);
  }
  degree = std::min(std::max(1, degree), std::max(1, static_cast<int>(n) - 2));
  arma::mat X = gniplr_poly_raw(x, degree);
  arma::rowvec mu = arma::mean(X, 0);
  arma::rowvec sd(X.n_cols, arma::fill::ones);
  for (arma::uword j = 0; j < X.n_cols; ++j) {
    const arma::vec centered = X.col(j) - mu[j];
    const double s = std::sqrt(arma::mean(arma::square(centered)));
    sd[j] = (std::isfinite(s) && s > 0.0) ? s : 1.0;
    X.col(j) = (X.col(j) - mu[j]) / sd[j];
  }
  arma::vec beta(X.n_cols, arma::fill::zeros);
  const double intercept = arma::mean(t);
  arma::vec yc = t - intercept;
  arma::vec resid = yc;
  const double lambda = std::max(0.0, alpha);
  const double inv_n = 1.0 / static_cast<double>(n);
  for (int iter = 0; iter < 1000; ++iter) {
    double max_delta = 0.0;
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      const double old = beta[j];
      resid += X.col(j) * old;
      const double rho = arma::dot(X.col(j), resid) * inv_n;
      const double z = arma::dot(X.col(j), X.col(j)) * inv_n;
      beta[j] = z > 0.0 ? gniplr_soft_threshold(rho, lambda) / z : 0.0;
      resid -= X.col(j) * beta[j];
      max_delta = std::max(max_delta, std::abs(beta[j] - old));
    }
    if (max_delta < 1e-7) break;
  }
  arma::mat P = gniplr_poly_raw(predict_target ? y : x, degree);
  for (arma::uword j = 0; j < P.n_cols; ++j) {
    P.col(j) = (P.col(j) - mu[j]) / sd[j];
  }
  arma::vec pred = intercept + P * beta;
  pred.replace(arma::datum::nan, 0.0);
  pred.replace(arma::datum::inf, 0.0);
  return pred;
}

double gniplr_rss(const arma::vec& y, const arma::mat& X) {
  arma::vec beta = arma::pinv(X) * y;
  arma::vec resid = X * beta - y;
  const double rss = arma::dot(resid, resid);
  return std::isfinite(rss) ? rss : std::numeric_limits<double>::infinity();
}

double gniplr_granger_lag(const arma::vec& a, const arma::vec& b, int lag) {
  const int n = static_cast<int>(a.n_elem);
  if (lag < 1 || n <= (2 * lag + 1)) return 1.0;
  const int rows = n - lag;
  arma::vec y(rows);
  arma::mat full(rows, 1 + 2 * lag, arma::fill::ones);
  arma::mat restricted(rows, 1 + lag, arma::fill::ones);
  for (int r = 0; r < rows; ++r) {
    y[r] = b[r];
    for (int l = 1; l <= lag; ++l) {
      full(r, l) = a[r + l];
      full(r, lag + l) = b[r + l];
      restricted(r, l) = b[r + l];
    }
  }
  const double rssu = gniplr_rss(y, full);
  const double rssr = gniplr_rss(y, restricted);
  const double df2 = static_cast<double>(n - (2 * lag + 1));
  const double denom = rssu / df2;
  if (!std::isfinite(rssu) || !std::isfinite(rssr) || denom <= 0.0) return 1.0;
  const double fval = ((rssr - rssu) / static_cast<double>(lag)) / denom;
  if (!std::isfinite(fval) || fval < 0.0) return 1.0;
  const double p = R::pf(fval, static_cast<double>(lag), df2, false, false);
  return std::isfinite(p) ? p : 1.0;
}

double gniplr_granger(const arma::vec& a, const arma::vec& b, int max_lag) {
  max_lag = std::max(1, std::min(3, max_lag));
  arma::vec aa = (300.0 * a) / 90001.0;
  arma::vec bb = (300.0 * b) / 90001.0;
  double best = 1.0;
  for (int lag = 1; lag <= max_lag; ++lag) {
    best = std::min(best, gniplr_granger_lag(aa, bb, lag));
  }
  return best;
}

} // namespace

// [[Rcpp::export]]
List gniplr_cpp(
  arma::mat expression,
  IntegerVector target_idx,
  double correlation_threshold = 0.3,
  int lasso_degree = 30,
  double lasso_alpha = 0.1,
  int max_lag = 3
) {
  const int gene_num = static_cast<int>(expression.n_rows);
  const int cell_num = static_cast<int>(expression.n_cols);
  max_lag = std::max(1, std::min(3, max_lag));
  const int min_cells = max_lag == 1 ? 4 : (max_lag == 2 ? 6 : 8);
  if (cell_num < min_cells) {
    stop("GNIPLR requires at least %d cells for max_lag=%d", min_cells, max_lag);
  }
  arma::mat corr(gene_num, gene_num, arma::fill::eye);
  for (int i = 0; i < gene_num; ++i) {
    for (int j = i + 1; j < gene_num; ++j) {
      const double c = gniplr_cor(expression.row(i), expression.row(j));
      corr(i, j) = c;
      corr(j, i) = c;
    }
  }
  std::vector<char> is_target(gene_num, 0);
  for (int k = 0; k < target_idx.size(); ++k) {
    const int idx = target_idx[k] - 1;
    if (idx >= 0 && idx < gene_num) is_target[idx] = 1;
  }
  arma::mat grn(gene_num, gene_num, arma::fill::zeros);
  std::vector<int> regulators;
  std::vector<int> targets;
  std::vector<double> pvalues;
  std::vector<double> importance;
  for (int i = 0; i < gene_num; ++i) {
    double mpmax = 0.0;
    for (int j = 0; j < gene_num; ++j) {
      if (i != j) mpmax = std::max(mpmax, std::abs(corr(i, j)));
    }
    if (mpmax <= 0.0) continue;
    for (int j = 0; j < gene_num; ++j) {
      if (i == j || !is_target[j]) continue;
      if (std::abs(corr(i, j)) < mpmax * correlation_threshold) continue;
      arma::vec a = gniplr_project(
        expression.row(i),
        expression.row(j),
        false,
        lasso_degree,
        lasso_alpha
      );
      arma::vec b = gniplr_project(
        expression.row(i),
        expression.row(j),
        true,
        lasso_degree,
        lasso_alpha
      );
      if (arma::approx_equal(a, b, "absdiff", 0.0)) continue;
      const double p = gniplr_granger(a, b, max_lag);
      if (!std::isfinite(p) || p <= 0.0) continue;
      grn(i, j) = p;
      regulators.push_back(i + 1);
      targets.push_back(j + 1);
      pvalues.push_back(p);
      importance.push_back(-std::log10(std::max(
        p,
        std::numeric_limits<double>::min()
      )));
    }
  }
  return List::create(
    _["adjacency"] = DataFrame::create(
      _["regulator"] = regulators,
      _["target"] = targets,
      _["importance"] = importance,
      _["pvalue"] = pvalues,
      _["stringsAsFactors"] = false
    ),
    _["grn_matrix"] = grn
  );
}
