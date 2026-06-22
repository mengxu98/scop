// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "cibersort_libsvm.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct CibersortCoreResult {
  arma::rowvec weights;
  double correlation;
  double rmse;
};

struct LinearSvrTrainingData {
  int n_genes;
  int n_cell_types;
  std::vector<cibersort_svm_node> nodes;
  std::vector<cibersort_svm_node*> row_ptrs;

  explicit LinearSvrTrainingData(const arma::mat& x)
      : n_genes(static_cast<int>(x.n_rows)),
        n_cell_types(static_cast<int>(x.n_cols)),
        nodes(static_cast<std::size_t>(n_genes) *
              static_cast<std::size_t>(n_cell_types + 1)),
        row_ptrs(static_cast<std::size_t>(n_genes)) {
    if (n_genes == 0 || n_cell_types == 0) {
      throw std::runtime_error("empty matrix passed to CIBERSORT SVR");
    }
    for (int i = 0; i < n_genes; ++i) {
      const std::size_t offset =
        static_cast<std::size_t>(i) * static_cast<std::size_t>(n_cell_types + 1);
      row_ptrs[static_cast<std::size_t>(i)] = &nodes[offset];
      for (int j = 0; j < n_cell_types; ++j) {
        nodes[offset + static_cast<std::size_t>(j)].index = j + 1;
        nodes[offset + static_cast<std::size_t>(j)].value = x(i, j);
      }
      nodes[offset + static_cast<std::size_t>(n_cell_types)].index = -1;
      nodes[offset + static_cast<std::size_t>(n_cell_types)].value = 0.0;
    }
  }
};

int cibersort_worker_count(int requested, int n_tasks) {
  if (requested <= 1 || n_tasks <= 1) {
    return 1;
  }
  int cores = std::min(requested, n_tasks);
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    cores = std::min(cores, static_cast<int>(hardware));
  }
  return std::max(1, cores);
}

std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

double sample_sd(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  if (n < 2) {
    return 0.0;
  }
  const double mu = arma::mean(x);
  double ss = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double d = x[i] - mu;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(n - 1));
}

arma::vec standardize_vector(const arma::vec& x) {
  arma::vec out = x;
  const double mu = arma::mean(out);
  const double sd = sample_sd(out);
  out -= mu;
  if (sd > 0.0 && std::isfinite(sd)) {
    out /= sd;
  } else {
    out.zeros();
  }
  return out;
}

double pearson_cor(const arma::vec& x, const arma::vec& y) {
  const arma::uword n = x.n_elem;
  if (n != y.n_elem || n < 2) {
    return NA_REAL;
  }
  const double mx = arma::mean(x);
  const double my = arma::mean(y);
  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double dx = x[i] - mx;
    const double dy = y[i] - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) {
    return NA_REAL;
  }
  return sxy / std::sqrt(sxx * syy);
}

double rmse(const arma::vec& x, const arma::vec& y) {
  if (x.n_elem != y.n_elem || x.n_elem == 0) {
    return NA_REAL;
  }
  double ss = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    const double d = x[i] - y[i];
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(x.n_elem));
}

arma::mat standardize_signature(const arma::mat& x) {
  arma::mat out = x;
  const arma::vec values = arma::vectorise(out);
  const double mu = arma::mean(values);
  const double sd = sample_sd(values);
  if (!(sd > 0.0) || !std::isfinite(sd)) {
    throw std::runtime_error("signature matrix has zero variance");
  }
  out -= mu;
  out /= sd;
  return out;
}

arma::mat quantile_normalize(const arma::mat& y) {
  const int n_rows = static_cast<int>(y.n_rows);
  const int n_cols = static_cast<int>(y.n_cols);
  if (n_rows == 0 || n_cols == 0) {
    return y;
  }

  std::vector<std::vector<std::pair<double, int> > > sorted(n_cols);
  std::vector<double> rank_means(n_rows, 0.0);
  for (int c = 0; c < n_cols; ++c) {
    sorted[c].reserve(n_rows);
    for (int r = 0; r < n_rows; ++r) {
      sorted[c].push_back(std::make_pair(y(r, c), r));
    }
    std::stable_sort(
      sorted[c].begin(),
      sorted[c].end(),
      [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        return a.first < b.first;
      }
    );
    for (int r = 0; r < n_rows; ++r) {
      rank_means[r] += sorted[c][r].first;
    }
  }
  for (int r = 0; r < n_rows; ++r) {
    rank_means[r] /= static_cast<double>(n_cols);
  }

  arma::mat out(n_rows, n_cols);
  for (int c = 0; c < n_cols; ++c) {
    for (int r = 0; r < n_rows; ++r) {
      out(sorted[c][r].second, c) = rank_means[r];
    }
  }
  return out;
}

cibersort_svm_parameter make_svr_parameter(double nu, int n_features) {
  cibersort_svm_parameter param;
  param.svm_type = NU_SVR;
  param.kernel_type = LINEAR;
  param.degree = 3;
  param.gamma = n_features > 0 ? 1.0 / static_cast<double>(n_features) : 0.0;
  param.coef0 = 0.0;
  param.cache_size = 40.0;
  param.eps = 0.001;
  param.C = 1.0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  param.nu = nu;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  return param;
}

arma::rowvec fit_linear_nusvr_weights(
    const LinearSvrTrainingData& training,
    const std::vector<double>& labels,
    double nu
) {
  cibersort_svm_problem prob;
  prob.l = training.n_genes;
  prob.y = const_cast<double*>(labels.data());
  prob.x = const_cast<cibersort_svm_node**>(training.row_ptrs.data());

  cibersort_svm_parameter param = make_svr_parameter(nu, training.n_cell_types);
  const char* check = cibersort_svm_check_parameter(&prob, &param);
  if (check != NULL) {
    throw std::runtime_error(std::string("LIBSVM parameter error: ") + check);
  }

  cibersort_svm_model* model = cibersort_svm_train(&prob, &param);
  if (model == NULL) {
    throw std::runtime_error("LIBSVM failed to train CIBERSORT model");
  }

  arma::rowvec weights(training.n_cell_types, arma::fill::zeros);
  for (int sv = 0; sv < model->l; ++sv) {
    const double coef = model->sv_coef[0][sv];
    const cibersort_svm_node* node = model->SV[sv];
    while (node->index != -1) {
      const int feature = node->index - 1;
      if (feature >= 0 && feature < training.n_cell_types) {
        weights[feature] += coef * node->value;
      }
      ++node;
    }
  }
  cibersort_svm_free_and_destroy_model(&model);
  return weights;
}

arma::rowvec normalize_weights(arma::rowvec weights) {
  for (arma::uword i = 0; i < weights.n_elem; ++i) {
    if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
      weights[i] = 0.0;
    }
  }
  const double total = arma::accu(weights);
  if (total > 0.0 && std::isfinite(total)) {
    weights /= total;
  } else if (weights.n_elem > 0) {
    weights.fill(1.0 / static_cast<double>(weights.n_elem));
  }
  return weights;
}

CibersortCoreResult cibersort_core(
    const arma::mat& x,
    const LinearSvrTrainingData& training,
    const arma::vec& y
) {
  static const double nus[3] = {0.25, 0.5, 0.75};
  CibersortCoreResult best;
  best.weights = arma::rowvec(x.n_cols, arma::fill::zeros);
  best.correlation = NA_REAL;
  best.rmse = std::numeric_limits<double>::infinity();
  std::vector<double> labels(static_cast<std::size_t>(training.n_genes));
  for (int i = 0; i < training.n_genes; ++i) {
    labels[static_cast<std::size_t>(i)] = y[i];
  }

  for (int i = 0; i < 3; ++i) {
    arma::rowvec weights = normalize_weights(fit_linear_nusvr_weights(training, labels, nus[i]));
    const arma::vec pred = x * weights.t();
    const double this_rmse = rmse(pred, y);
    const double this_cor = pearson_cor(pred, y);
    if (std::isfinite(this_rmse) && this_rmse < best.rmse) {
      best.weights = weights;
      best.rmse = this_rmse;
      best.correlation = this_cor;
    }
  }

  if (!std::isfinite(best.rmse)) {
    best.rmse = NA_REAL;
  }
  return best;
}

arma::vec sample_null_vector(
    const arma::vec& pool,
    int n_genes,
    std::uint64_t seed
) {
  const int n_pool = static_cast<int>(pool.n_elem);
  if (n_genes > n_pool) {
    throw std::runtime_error("permutation pool is smaller than gene count");
  }
  std::mt19937_64 rng(seed);
  std::unordered_map<int, int> swapped;
  swapped.reserve(static_cast<std::size_t>(n_genes) * 2U);
  auto value_at = [&swapped](int key) {
    const std::unordered_map<int, int>::const_iterator it = swapped.find(key);
    return it == swapped.end() ? key : it->second;
  };

  arma::vec out(n_genes);
  for (int i = 0; i < n_genes; ++i) {
    std::uniform_int_distribution<int> dist(i, n_pool - 1);
    const int j = dist(rng);
    const int selected = value_at(j);
    swapped[j] = value_at(i);
    out[i] = pool[selected];
  }
  return standardize_vector(out);
}

template <typename Function>
void parallel_for_tasks(int n_tasks, int cores, Function fn) {
  if (n_tasks <= 0) {
    return;
  }
  const int workers = cibersort_worker_count(cores, n_tasks);
  if (workers <= 1) {
    for (int i = 0; i < n_tasks; ++i) {
      fn(i);
    }
    return;
  }

  std::atomic<int> next(0);
  std::exception_ptr error;
  std::mutex error_mutex;
  std::atomic<bool> failed(false);
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (int w = 0; w < workers; ++w) {
    pool.emplace_back([&]() {
      try {
        while (!failed.load()) {
          const int task = next.fetch_add(1);
          if (task >= n_tasks) {
            break;
          }
          fn(task);
        }
      } catch (...) {
        failed.store(true);
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!error) {
          error = std::current_exception();
        }
      }
    });
  }
  for (std::thread& worker : pool) {
    worker.join();
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

std::vector<double> cibersort_null_distribution(
    const arma::mat& x,
    const LinearSvrTrainingData& training,
    const arma::mat& y,
    int perm,
    int cores,
    int seed
) {
  std::vector<double> null_dist(perm);
  if (perm <= 0) {
    return null_dist;
  }
  const arma::vec pool = arma::vectorise(y);
  const int n_genes = static_cast<int>(x.n_rows);
  parallel_for_tasks(perm, cores, [&](int p) {
    const std::uint64_t task_seed = splitmix64(
      static_cast<std::uint64_t>(seed) + static_cast<std::uint64_t>(p + 1)
    );
    const arma::vec yr = sample_null_vector(pool, n_genes, task_seed);
    null_dist[p] = cibersort_core(x, training, yr).correlation;
  });
  std::sort(null_dist.begin(), null_dist.end());
  return null_dist;
}

double cibersort_p_value(double cor_value, const std::vector<double>& null_dist) {
  const int n = static_cast<int>(null_dist.size());
  if (n == 0 || !std::isfinite(cor_value)) {
    return 9999.0;
  }
  std::vector<double>::const_iterator upper =
    std::lower_bound(null_dist.begin(), null_dist.end(), cor_value);
  int best = 0;
  if (upper == null_dist.begin()) {
    best = 0;
  } else if (upper == null_dist.end()) {
    best = n - 1;
  } else {
    const int hi = static_cast<int>(upper - null_dist.begin());
    const int lo = hi - 1;
    best = std::fabs(null_dist[hi] - cor_value) < std::fabs(null_dist[lo] - cor_value)
      ? hi
      : lo;
  }
  return 1.0 - static_cast<double>(best + 1) / static_cast<double>(n);
}

}  // namespace

// [[Rcpp::export]]
List cibersort_cpp(
    NumericMatrix signature,
    NumericMatrix mixture,
    int perm = 0,
    bool QN = true,
    bool absolute = false,
    int cores = 1,
    int seed = 123,
    bool verbose = false
) {
  (void)absolute;
  (void)verbose;
  const int n_genes = signature.nrow();
  const int n_cell_types = signature.ncol();
  const int n_samples = mixture.ncol();
  if (n_genes <= 0 || n_cell_types <= 0 || n_samples <= 0) {
    stop("signature and mixture matrices must be non-empty");
  }
  if (mixture.nrow() != n_genes) {
    stop("signature and mixture matrices must have the same number of rows");
  }
  if (perm < 0) {
    stop("perm must be non-negative");
  }

  arma::mat x(signature.begin(), n_genes, n_cell_types, true);
  arma::mat y(mixture.begin(), n_genes, n_samples, true);
  x.replace(arma::datum::nan, 0.0);
  x.replace(arma::datum::inf, 0.0);
  y.replace(arma::datum::nan, 0.0);
  y.replace(arma::datum::inf, 0.0);

  const double y_max = y.max();
  if (std::isfinite(y_max) && y_max < 50.0) {
    y = arma::exp2(y);
  }
  if (QN) {
    y = quantile_normalize(y);
  }
  const arma::mat x_scaled = standardize_signature(x);
  const LinearSvrTrainingData training(x_scaled);

  if (perm > 0) {
    Rcpp::checkUserInterrupt();
  }
  const std::vector<double> null_dist = cibersort_null_distribution(
    x_scaled,
    training,
    y,
    perm,
    cores,
    seed
  );
  Rcpp::checkUserInterrupt();

  arma::mat fractions(n_samples, n_cell_types, arma::fill::zeros);
  arma::mat stats(n_samples, 3, arma::fill::zeros);
  parallel_for_tasks(n_samples, cores, [&](int s) {
    const arma::vec y_scaled = standardize_vector(y.col(s));
    const CibersortCoreResult result = cibersort_core(x_scaled, training, y_scaled);
    fractions.row(s) = result.weights;
    stats(s, 0) = cibersort_p_value(result.correlation, null_dist);
    stats(s, 1) = result.correlation;
    stats(s, 2) = result.rmse;
  });

  NumericMatrix fraction_out(n_samples, n_cell_types);
  NumericMatrix stats_out(n_samples, 3);
  for (int i = 0; i < n_samples; ++i) {
    for (int j = 0; j < n_cell_types; ++j) {
      fraction_out(i, j) = fractions(i, j);
    }
    stats_out(i, 0) = stats(i, 0);
    stats_out(i, 1) = stats(i, 1);
    stats_out(i, 2) = stats(i, 2);
  }

  List sig_dimnames = signature.attr("dimnames");
  List mix_dimnames = mixture.attr("dimnames");
  if (sig_dimnames.size() >= 2 && !Rf_isNull(sig_dimnames[1]) &&
      mix_dimnames.size() >= 2 && !Rf_isNull(mix_dimnames[1])) {
    fraction_out.attr("dimnames") = List::create(mix_dimnames[1], sig_dimnames[1]);
    stats_out.attr("dimnames") = List::create(
      mix_dimnames[1],
      CharacterVector::create("P-value", "Correlation", "RMSE")
    );
  }

  return List::create(
    _["proportion_matrix"] = fraction_out,
    _["statistics"] = stats_out
  );
}
