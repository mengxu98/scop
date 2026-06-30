// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

static const double SCENIC_FEATURE_THRESHOLD = 1e-7;
static const double SCENIC_PROXY_REL_TOL = 64.0 * std::numeric_limits<double>::epsilon();
static const int SCENIC_RAND_R_MAX = 2147483647;

struct ScenicEdge {
  int regulator;
  int target;
  double importance;
};

struct ScenicSplit {
  int feature;
  double threshold;
  double proxy_gain;
  double importance_gain;
};

struct ScenicTreeNode {
  int feature;
  double threshold;
  double value;
  int left;
  int right;
  int n_rows;
  double sse;
  double importance_gain;
};

struct ScenicGrnProfile {
  double prepare_seconds = 0.0;
  double feature_order_seconds = 0.0;
  double target_setup_seconds = 0.0;
  double sample_seconds = 0.0;
  double build_tree_seconds = 0.0;
  double best_split_seconds = 0.0;
  double oob_predict_seconds = 0.0;
  double update_predict_seconds = 0.0;
  double edge_seconds = 0.0;
  int targets = 0;
  int rounds = 0;
  int trees = 0;
  int nodes = 0;
  int best_split_calls = 0;
  int split_feature_visits = 0;
  int split_threshold_visits = 0;
};

struct ScenicCandidateTrace {
  int trace_node = 1;
  bool collect = false;
  std::vector<int> visit_order;
  std::vector<int> feature;
  std::vector<double> proxy_gain;
  std::vector<double> importance_gain;
  std::vector<double> threshold;
  std::vector<int> n_left;
  std::vector<int> n_right;
};

struct ScenicSparseDgcMatrix {
  IntegerVector p;
  IntegerVector i;
  NumericVector x;
  IntegerVector dim;
  int n_rows;
  int n_cols;

  explicit ScenicSparseDgcMatrix(S4 mat) :
      p(mat.slot("p")),
      i(mat.slot("i")),
      x(mat.slot("x")),
      dim(mat.slot("Dim")) {
    if (!mat.is("dgCMatrix")) {
      stop("expr must be a dgCMatrix.");
    }
    if (dim.size() != 2) {
      stop("expr must have two dimensions.");
    }
    n_rows = dim[0];
    n_cols = dim[1];
    if (p.size() != n_cols + 1) {
      stop("Invalid dgCMatrix column pointer.");
    }
    if (i.size() != x.size()) {
      stop("Invalid dgCMatrix index/value slots.");
    }
  }

  inline int nrow() const { return n_rows; }
  inline int ncol() const { return n_cols; }
};

struct ScenicScopeTimer {
  double* seconds;
  std::chrono::steady_clock::time_point start;

  explicit ScenicScopeTimer(double* seconds_) : seconds(seconds_) {
    if (seconds != nullptr) start = std::chrono::steady_clock::now();
  }

  ~ScenicScopeTimer() {
    if (seconds == nullptr) return;
    const std::chrono::duration<double> elapsed =
      std::chrono::steady_clock::now() - start;
    *seconds += elapsed.count();
  }
};

static bool scenic_proxy_improves(double candidate, double current) {
  if (!R_finite(candidate)) return false;
  if (!R_finite(current) || current <= 0.0) return candidate > current;
  const double scale = std::max(std::fabs(candidate), std::fabs(current));
  const double tol = SCENIC_PROXY_REL_TOL * std::max(1.0, scale);
  return candidate > current + tol;
}

static bool scenic_edge_before(const ScenicEdge& a, const ScenicEdge& b) {
  if (a.target != b.target) return a.target < b.target;
  if (a.importance != b.importance) return a.importance > b.importance;
  return a.regulator < b.regulator;
}

static double scenic_sse_from_sums(double sum, double sum_sq, int n) {
  if (n <= 0) return 0.0;
  const double sse = sum_sq - sum * sum / static_cast<double>(n);
  return sse > 0.0 ? sse : 0.0;
}

static double scenic_tree_feature_value(const NumericMatrix& expr, int row, int feature) {
  double x = expr.begin()[row + feature * expr.nrow()];
  if (!R_finite(x)) x = 0.0;
  return static_cast<double>(static_cast<float>(x));
}

static double scenic_tree_feature_value(const ScenicSparseDgcMatrix& expr, int row, int feature) {
  if (row < 0 || row >= expr.nrow() || feature < 0 || feature >= expr.ncol()) {
    return 0.0;
  }
  const int* p = INTEGER(expr.p);
  const int* i = INTEGER(expr.i);
  const double* x = REAL(expr.x);
  const int begin = p[feature];
  const int end = p[feature + 1];
  const int* found = std::lower_bound(i + begin, i + end, row);
  if (found == i + end || *found != row) return 0.0;
  double value = x[found - i];
  if (!R_finite(value)) value = 0.0;
  return static_cast<double>(static_cast<float>(value));
}

static double scenic_node_mean(
    const std::vector<int>& rows,
    const std::vector<double>& residual) {
  if (rows.empty()) return 0.0;
  double sum = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) sum += residual[rows[i]];
  return sum / static_cast<double>(rows.size());
}

static double scenic_node_mean(
    const std::vector<int>& samples,
    int start,
    int end,
    const std::vector<double>& residual) {
  if (end <= start) return 0.0;
  double sum = 0.0;
  for (int i = start; i < end; ++i) sum += residual[samples[i]];
  return sum / static_cast<double>(end - start);
}

static double scenic_node_sse(
    const std::vector<int>& samples,
    int start,
    int end,
    const std::vector<double>& residual,
    double* sum_out = nullptr,
    double* sum_sq_out = nullptr) {
  double sum = 0.0;
  double sum_sq = 0.0;
  for (int i = start; i < end; ++i) {
    const double y = residual[samples[i]];
    sum += y;
    sum_sq += y * y;
  }
  if (sum_out != nullptr) *sum_out = sum;
  if (sum_sq_out != nullptr) *sum_sq_out = sum_sq;
  return scenic_sse_from_sums(sum, sum_sq, end - start);
}

template <typename Expr>
static void scenic_sort_samples_by_feature(
    const Expr& expr,
    std::vector<int>& samples,
    std::vector<double>& feature_values,
    int start,
    int end,
    int feature) {
  std::vector<std::pair<double, int> > values;
  values.reserve(end - start);
  for (int i = start; i < end; ++i) {
    values.push_back(std::make_pair(scenic_tree_feature_value(expr, samples[i], feature), samples[i]));
  }
  std::sort(values.begin(), values.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  });
  for (int i = start; i < end; ++i) {
    const std::pair<double, int>& value = values[static_cast<std::size_t>(i - start)];
    feature_values[i] = value.first;
    samples[i] = value.second;
  }
}

template <typename Expr>
static void scenic_order_samples_by_feature(
    const Expr& expr,
    const std::vector<std::vector<int> >& feature_orders,
    const std::vector<std::vector<double> >& feature_order_values,
    const std::vector<std::vector<int> >& feature_order_ranks,
    std::vector<int>& samples,
    std::vector<double>& feature_values,
    int start,
    int end,
    int feature,
    std::vector<int>& row_marks,
    int row_mark) {
  if (
      feature < 0 ||
      feature >= static_cast<int>(feature_orders.size()) ||
      feature >= static_cast<int>(feature_order_values.size()) ||
      feature_orders[feature].empty()) {
    scenic_sort_samples_by_feature(expr, samples, feature_values, start, end, feature);
    return;
  }
  int pos = start;
  const std::vector<int>& order = feature_orders[feature];
  const std::vector<double>& values = feature_order_values[feature];
  if (
      static_cast<int>(order.size()) < expr.nrow() &&
      values.size() == order.size()) {
    std::vector<int> original(samples.begin() + start, samples.begin() + end);
    std::vector<std::pair<int, double> > positive_nonzero;
    std::vector<int> emitted_nonzero;
    positive_nonzero.reserve(order.size());
    emitted_nonzero.reserve(order.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
      const int row = order[i];
      if (row_marks[row] != row_mark) continue;
      row_marks[row] = row_mark + 1;
      emitted_nonzero.push_back(row);
      const double value = values[i];
      if (value < -SCENIC_FEATURE_THRESHOLD) {
        feature_values[pos] = value;
        samples[pos] = row;
        ++pos;
      } else {
        positive_nonzero.push_back(std::make_pair(row, value));
      }
    }
    for (std::size_t i = 0; i < original.size(); ++i) {
      const int row = original[i];
      if (row_marks[row] == row_mark) {
        feature_values[pos] = 0.0;
        samples[pos] = row;
        ++pos;
      }
    }
    for (std::size_t i = 0; i < positive_nonzero.size(); ++i) {
      feature_values[pos] = positive_nonzero[i].second;
      samples[pos] = positive_nonzero[i].first;
      ++pos;
    }
    for (std::size_t i = 0; i < emitted_nonzero.size(); ++i) {
      row_marks[emitted_nonzero[i]] = row_mark;
    }
    if (pos != end) {
      scenic_sort_samples_by_feature(expr, samples, feature_values, start, end, feature);
    }
    return;
  }
  if (
      feature < static_cast<int>(feature_order_ranks.size()) &&
      !feature_order_ranks[feature].empty() &&
      (end - start) * 8 < static_cast<int>(order.size())) {
    const std::vector<int>& ranks = feature_order_ranks[feature];
    std::sort(samples.begin() + start, samples.begin() + end, [&](int a, int b) {
      return ranks[a] < ranks[b];
    });
    for (int i = start; i < end; ++i) {
      feature_values[i] = values[ranks[samples[i]]];
    }
    return;
  }
  for (std::size_t i = 0; i < order.size() && pos < end; ++i) {
    const int row = order[i];
    if (row_marks[row] == row_mark) {
      feature_values[pos] = values[i];
      samples[pos] = row;
      ++pos;
    }
  }
  if (pos != end) {
    scenic_sort_samples_by_feature(expr, samples, feature_values, start, end, feature);
  }
}

template <typename Expr>
static int scenic_partition_samples_final(
    const Expr& expr,
    std::vector<int>& samples,
    int start,
    int end,
    int feature,
    double threshold) {
  int p = start;
  int partition_end = end - 1;
  while (p < partition_end) {
    if (scenic_tree_feature_value(expr, samples[p], feature) <= threshold) {
      ++p;
    } else {
      std::swap(samples[p], samples[partition_end]);
      --partition_end;
    }
  }
  int n_left = 0;
  for (int i = start; i < end; ++i) {
    if (scenic_tree_feature_value(expr, samples[i], feature) <= threshold) ++n_left;
  }
  return start + n_left;
}

static int scenic_rand_int(std::uint32_t& split_seed) {
  if (split_seed == 0U) split_seed = 1U;
  split_seed ^= static_cast<std::uint32_t>(split_seed << 13);
  split_seed ^= static_cast<std::uint32_t>(split_seed >> 17);
  split_seed ^= static_cast<std::uint32_t>(split_seed << 5);
  return static_cast<int>(
    split_seed % (static_cast<std::uint32_t>(SCENIC_RAND_R_MAX) + 1U)
  );
}

static double scenic_numpy_random_sample(std::mt19937& rng) {
  const std::uint32_t a = static_cast<std::uint32_t>(rng()) >> 5;
  const std::uint32_t b = static_cast<std::uint32_t>(rng()) >> 6;
  return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) /
    9007199254740992.0;
}

static int scenic_numpy_randint(std::mt19937& rng, int low, int high) {
  const std::uint32_t max_value = static_cast<std::uint32_t>(high - low - 1);
  std::uint32_t mask = max_value;
  mask |= mask >> 1;
  mask |= mask >> 2;
  mask |= mask >> 4;
  mask |= mask >> 8;
  mask |= mask >> 16;
  std::uint32_t value = 0U;
  do {
    value = static_cast<std::uint32_t>(rng()) & mask;
  } while (value > max_value);
  return low + static_cast<int>(value);
}

static std::vector<unsigned char> scenic_sample_mask(
    int n_samples,
    int n_inbag,
    std::mt19937& rng) {
  std::vector<unsigned char> mask(n_samples, 0);
  if (n_inbag >= n_samples) {
    std::fill(mask.begin(), mask.end(), 1);
    return mask;
  }
  int n_bagged = 0;
  for (int i = 0; i < n_samples; ++i) {
    if (scenic_numpy_random_sample(rng) * static_cast<double>(n_samples - i) <
        static_cast<double>(n_inbag - n_bagged)) {
      mask[i] = 1;
      ++n_bagged;
    }
  }
  return mask;
}

static void scenic_sample_mask_fill(
    int n_samples,
    int n_inbag,
    std::mt19937& rng,
    std::vector<unsigned char>& mask,
    std::vector<int>& train_rows) {
  mask.assign(n_samples, 0);
  train_rows.clear();
  train_rows.reserve(n_inbag);
  if (n_inbag >= n_samples) {
    std::fill(mask.begin(), mask.end(), 1);
    train_rows.resize(n_samples);
    std::iota(train_rows.begin(), train_rows.end(), 0);
    return;
  }
  int n_bagged = 0;
  for (int i = 0; i < n_samples; ++i) {
    if (scenic_numpy_random_sample(rng) * static_cast<double>(n_samples - i) <
        static_cast<double>(n_inbag - n_bagged)) {
      mask[i] = 1;
      train_rows.push_back(i);
      ++n_bagged;
    }
  }
}

template <typename Expr>
static std::vector<std::vector<int> > scenic_feature_orders(
    const Expr& expr,
    const std::vector<int>& features,
    std::vector<std::vector<double> >* order_values = nullptr,
    std::vector<std::vector<int> >* order_ranks = nullptr) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  std::vector<std::vector<int> > orders(n_genes);
  if (order_values != nullptr) {
    order_values->clear();
    order_values->resize(n_genes);
  }
  if (order_ranks != nullptr) {
    order_ranks->clear();
    order_ranks->resize(n_genes);
  }
  for (std::size_t fi = 0; fi < features.size(); ++fi) {
    const int g = features[fi];
    if (g < 0 || g >= n_genes || !orders[g].empty()) continue;
    std::vector<int> order(n_samples);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      double xa = scenic_tree_feature_value(expr, a, g);
      double xb = scenic_tree_feature_value(expr, b, g);
      if (xa != xb) return xa < xb;
      return a < b;
    });
    if (order_values != nullptr) {
      std::vector<double> values(n_samples);
      for (int i = 0; i < n_samples; ++i) {
        values[i] = scenic_tree_feature_value(expr, order[i], g);
      }
      (*order_values)[g] = values;
    }
    if (order_ranks != nullptr) {
      std::vector<int> ranks(n_samples);
      for (int i = 0; i < n_samples; ++i) {
        ranks[order[i]] = i;
      }
      (*order_ranks)[g] = ranks;
    }
    orders[g] = order;
  }
  return orders;
}

static std::vector<std::vector<int> > scenic_feature_orders(
    const ScenicSparseDgcMatrix& expr,
    const std::vector<int>& features,
    std::vector<std::vector<double> >* order_values = nullptr,
    std::vector<std::vector<int> >* order_ranks = nullptr) {
  std::vector<std::vector<int> > orders(expr.ncol());
  if (order_values != nullptr) {
    order_values->clear();
    order_values->resize(expr.ncol());
  }
  if (order_ranks != nullptr) {
    order_ranks->clear();
    order_ranks->resize(expr.ncol());
  }
  const int* p = INTEGER(expr.p);
  const int* row_idx = INTEGER(expr.i);
  const double* x = REAL(expr.x);
  for (std::size_t fi = 0; fi < features.size(); ++fi) {
    const int g = features[fi];
    if (g < 0 || g >= expr.ncol() || !orders[g].empty()) continue;
    std::vector<std::pair<double, int> > nonzero;
    nonzero.reserve(p[g + 1] - p[g]);
    for (int ptr = p[g]; ptr < p[g + 1]; ++ptr) {
      double value = x[ptr];
      if (!R_finite(value) || std::fabs(value) <= SCENIC_FEATURE_THRESHOLD) continue;
      nonzero.push_back(std::make_pair(
        static_cast<double>(static_cast<float>(value)),
        row_idx[ptr]
      ));
    }
    std::sort(nonzero.begin(), nonzero.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });
    std::vector<int> order(nonzero.size());
    std::vector<double> values(nonzero.size());
    for (std::size_t i = 0; i < nonzero.size(); ++i) {
      order[i] = nonzero[i].second;
      values[i] = nonzero[i].first;
    }
    orders[g] = order;
    if (order_values != nullptr) {
      (*order_values)[g] = values;
    }
  }
  return orders;
}

template <typename Expr>
static ScenicSplit scenic_best_split(
    const Expr& expr,
    const std::vector<double>& residual,
    std::vector<int>& samples,
    std::vector<double>& feature_values,
    int start,
    int end,
    std::vector<int>& feature_pool,
    std::vector<int>& constant_features,
    int n_known_constants,
    int& n_total_constants,
    int mtry,
    std::uint32_t& split_seed,
    const std::vector<std::vector<int> >& feature_orders,
    const std::vector<std::vector<double> >& feature_order_values,
    const std::vector<std::vector<int> >& feature_order_ranks,
    std::vector<int>& row_marks,
    int& row_mark,
    ScenicGrnProfile* profile,
    ScenicCandidateTrace* candidate_trace = nullptr) {
  if (profile != nullptr) ++profile->best_split_calls;
  ScenicSplit best{-1, 0.0, 0.0, 0.0};
  if (end - start < 2 || feature_pool.empty()) return best;

  double parent_sum = 0.0;
  double parent_sum_sq = 0.0;
  const double parent_sse = scenic_node_sse(samples, start, end, residual, &parent_sum, &parent_sum_sq);
  if (parent_sse <= 0.0) return best;
  ++row_mark;
  if (row_mark == std::numeric_limits<int>::max()) {
    std::fill(row_marks.begin(), row_marks.end(), 0);
    row_mark = 1;
  }
  for (int i = start; i < end; ++i) {
    row_marks[samples[i]] = row_mark;
  }

  const int n_features = static_cast<int>(feature_pool.size());
  int f_i = n_features;
  int n_visited = 0;
  int n_found_constants = 0;
  int n_drawn_constants = 0;
  n_known_constants = std::max(0, std::min(n_known_constants, n_features));
  n_total_constants = n_known_constants;
  while (
      f_i > n_total_constants &&
      (n_visited < mtry ||
       n_visited <= n_found_constants + n_drawn_constants)) {
    int f_j = n_drawn_constants +
      scenic_rand_int(split_seed) % (f_i - n_found_constants - n_drawn_constants);
    if (f_j < n_known_constants) {
      std::swap(feature_pool[n_drawn_constants], feature_pool[f_j]);
      ++n_drawn_constants;
      ++n_visited;
      continue;
    }
    f_j += n_found_constants;
    const int feature = feature_pool[f_j];
    ++n_visited;
    if (profile != nullptr) ++profile->split_feature_visits;

    scenic_order_samples_by_feature(
      expr, feature_orders, feature_order_values, feature_order_ranks,
      samples, feature_values, start, end, feature,
      row_marks, row_mark);
    const double min_x = feature_values[start];
    const double max_x = feature_values[end - 1];
    if (max_x <= min_x + SCENIC_FEATURE_THRESHOLD) {
      std::swap(feature_pool[f_j], feature_pool[n_total_constants]);
      ++n_found_constants;
      ++n_total_constants;
      continue;
    }
    --f_i;
    std::swap(feature_pool[f_i], feature_pool[f_j]);

    double feature_best_proxy = -std::numeric_limits<double>::infinity();
    double feature_best_importance = 0.0;
    double feature_best_threshold = 0.0;
    int feature_best_left = 0;
    int feature_best_right = 0;
    bool have_group = false;
    double current_x = 0.0;
    double group_sum = 0.0;
    double group_sum_sq = 0.0;
    int group_n = 0;
    double left_sum = 0.0;
    double left_sum_sq = 0.0;
    int n_left = 0;

    for (int pos = start; pos < end; ++pos) {
      const int row = samples[pos];
      const double x = feature_values[pos];
      const double y = residual[row];
      if (!have_group) {
        have_group = true;
        current_x = x;
        group_sum = y;
        group_sum_sq = y * y;
        group_n = 1;
        continue;
      }
      if (x <= current_x + SCENIC_FEATURE_THRESHOLD) {
        group_sum += y;
        group_sum_sq += y * y;
        ++group_n;
        continue;
      }

      left_sum += group_sum;
      left_sum_sq += group_sum_sq;
      n_left += group_n;
      const int n_right = end - start - n_left;
      if (n_left < 1 || n_right < 1) continue;
      if (profile != nullptr) ++profile->split_threshold_visits;

      const double right_sum = parent_sum - left_sum;
      const double right_sum_sq = parent_sum_sq - left_sum_sq;
      const double left_sse = scenic_sse_from_sums(left_sum, left_sum_sq, n_left);
      const double right_sse = scenic_sse_from_sums(right_sum, right_sum_sq, n_right);
      const double importance_gain = parent_sse - left_sse - right_sse;
      const double diff =
        static_cast<double>(n_right) * left_sum - static_cast<double>(n_left) * right_sum;
      const double proxy_gain =
        2.0 * diff * diff /
        (static_cast<double>(n_left + n_right) *
         static_cast<double>(n_left) * static_cast<double>(n_right));
      const double threshold = (current_x + x) / 2.0;
      double split_threshold = threshold;
      if (split_threshold == static_cast<double>(x) || !R_finite(split_threshold)) {
        split_threshold = current_x;
      }
      if (scenic_proxy_improves(proxy_gain, feature_best_proxy)) {
        feature_best_proxy = proxy_gain;
        feature_best_importance = importance_gain;
        feature_best_threshold = split_threshold;
        feature_best_left = n_left;
        feature_best_right = n_right;
      }
      if (scenic_proxy_improves(proxy_gain, best.proxy_gain)) {
        best.feature = feature;
        best.threshold = split_threshold;
        best.proxy_gain = proxy_gain;
        best.importance_gain = importance_gain;
      }
      current_x = x;
      group_sum = y;
      group_sum_sq = y * y;
      group_n = 1;
    }
    if (candidate_trace != nullptr && candidate_trace->collect &&
        feature_best_proxy > -std::numeric_limits<double>::infinity()) {
      candidate_trace->visit_order.push_back(n_visited);
      candidate_trace->feature.push_back(feature + 1);
      candidate_trace->proxy_gain.push_back(feature_best_proxy);
      candidate_trace->importance_gain.push_back(feature_best_importance);
      candidate_trace->threshold.push_back(feature_best_threshold);
      candidate_trace->n_left.push_back(feature_best_left);
      candidate_trace->n_right.push_back(feature_best_right);
    }
  }
  for (int i = 0; i < n_known_constants; ++i) {
    feature_pool[i] = constant_features[i];
  }
  for (int i = n_known_constants; i < n_total_constants; ++i) {
    constant_features[i] = feature_pool[i];
  }
  return best;
}

template <typename Expr>
static int scenic_build_tree(
    const Expr& expr,
    const std::vector<double>& residual,
    std::vector<int>& samples,
    std::vector<double>& feature_values,
    int start,
    int end,
    std::vector<int>& feature_pool,
    std::vector<int>& constant_features,
    int n_known_constants,
    int depth,
    int max_depth,
    int mtry,
    std::uint32_t& split_seed,
    const std::vector<std::vector<int> >& feature_orders,
    const std::vector<std::vector<double> >& feature_order_values,
    const std::vector<std::vector<int> >& feature_order_ranks,
    std::vector<int>& row_marks,
    int& row_mark,
    std::vector<ScenicTreeNode>& nodes,
    std::vector<double>& importance,
    ScenicGrnProfile* profile,
    ScenicCandidateTrace* candidate_trace = nullptr) {
  ScenicTreeNode node;
  node.feature = -1;
  node.threshold = 0.0;
  node.value = scenic_node_mean(samples, start, end, residual);
  node.left = -1;
  node.right = -1;
  node.n_rows = end - start;
  double node_sum = 0.0;
  double node_sum_sq = 0.0;
  node.sse = scenic_node_sse(samples, start, end, residual, &node_sum, &node_sum_sq);
  node.importance_gain = 0.0;
  const int node_index = static_cast<int>(nodes.size());
  nodes.push_back(node);
  if (profile != nullptr) ++profile->nodes;

  if (depth >= max_depth || end - start < 2 || feature_pool.empty()) return node_index;

  int n_total_constants = n_known_constants;
  const bool old_collect = candidate_trace != nullptr && candidate_trace->collect;
  if (candidate_trace != nullptr) {
    candidate_trace->collect = (node_index + 1 == candidate_trace->trace_node);
  }
  const ScenicSplit split = scenic_best_split(
    expr, residual, samples, feature_values, start, end, feature_pool, constant_features,
    n_known_constants, n_total_constants, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks,
    row_marks, row_mark, profile, candidate_trace);
  if (candidate_trace != nullptr) {
    candidate_trace->collect = old_collect;
  }
  if (split.feature < 0) return node_index;

  const int split_pos = scenic_partition_samples_final(expr, samples, start, end, split.feature, split.threshold);
  if (split_pos <= start || split_pos >= end) return node_index;

  importance[split.feature] += split.importance_gain;
  const int left = scenic_build_tree(
    expr, residual, samples, feature_values, start, split_pos, feature_pool, constant_features,
    n_total_constants, depth + 1, max_depth, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks, row_marks, row_mark, nodes, importance, profile,
    candidate_trace);
  const int right = scenic_build_tree(
    expr, residual, samples, feature_values, split_pos, end, feature_pool, constant_features,
    n_total_constants, depth + 1, max_depth, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks, row_marks, row_mark, nodes, importance, profile,
    candidate_trace);
  nodes[node_index].feature = split.feature;
  nodes[node_index].threshold = split.threshold;
  nodes[node_index].left = left;
  nodes[node_index].right = right;
  nodes[node_index].importance_gain = split.importance_gain;
  return node_index;
}

template <typename Expr>
static double scenic_predict_tree(
    const Expr& expr,
    const std::vector<ScenicTreeNode>& nodes,
    int row) {
  int node_index = 0;
  while (node_index >= 0) {
    const ScenicTreeNode& node = nodes[node_index];
    if (node.feature < 0 || node.left < 0 || node.right < 0) return node.value;
    double x = scenic_tree_feature_value(expr, row, node.feature);
    node_index = (x <= node.threshold) ? node.left : node.right;
  }
  return 0.0;
}

static void scenic_profile_add_elapsed(
    ScenicGrnProfile* profile,
    double ScenicGrnProfile::*slot,
    const std::chrono::steady_clock::time_point& start) {
  if (profile == nullptr) return;
  const std::chrono::duration<double> elapsed =
    std::chrono::steady_clock::now() - start;
  profile->*slot += elapsed.count();
}

static int scenic_grn_worker_count(int requested, int n_tasks) {
  if (n_tasks <= 1) return 1;
  int cores = requested < 1 ? 1 : requested;
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    cores = std::min(cores, static_cast<int>(hardware));
  }
  return std::max(1, std::min(cores, n_tasks));
}

template <typename Expr>
static void scenic_grnboost_run_targets(
    const Expr& expr,
    const std::vector<int>& regs,
    const std::vector<int>& targets,
    std::size_t target_begin,
    std::size_t target_end,
    const std::vector<double>& means,
    const std::vector<std::vector<int> >& feature_orders,
    const std::vector<std::vector<double> >& feature_order_values,
    const std::vector<std::vector<int> >& feature_order_ranks,
    int n_rounds,
    double learning_rate,
    int max_edges_per_target,
    int max_depth,
    double max_features,
    double subsample,
    int early_stop_window_length,
    int random_seed,
    bool exclude_self,
    ScenicGrnProfile* profile,
    bool check_interrupt,
    std::vector<ScenicEdge>& edges) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  std::vector<double> residual(n_samples);
  std::vector<double> importance(n_genes);
  std::vector<double> tree_importance(n_genes);
  std::vector<double> feature_values(n_samples);

  for (std::size_t ti = target_begin; ti < target_end; ++ti) {
    std::mt19937 rng(static_cast<unsigned int>(random_seed));
    const int target = targets[ti];
    if (profile != nullptr) ++profile->targets;
    std::chrono::steady_clock::time_point target_setup_start;
    if (profile != nullptr) target_setup_start = std::chrono::steady_clock::now();

    std::fill(importance.begin(), importance.end(), 0.0);
    for (int i = 0; i < n_samples; ++i) {
      double y = scenic_tree_feature_value(expr, i, target);
      residual[i] = y - means[target];
    }
    std::vector<int> features;
    features.reserve(regs.size());
    for (std::size_t ri = 0; ri < regs.size(); ++ri) {
      if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
    }
    if (features.empty()) {
      scenic_profile_add_elapsed(
        profile, &ScenicGrnProfile::target_setup_seconds, target_setup_start);
      continue;
    }
    scenic_profile_add_elapsed(
      profile, &ScenicGrnProfile::target_setup_seconds, target_setup_start);
    const int mtry = std::min<int>(
      features.size(),
      std::max(1, static_cast<int>(std::floor(max_features * features.size())))
    );
    std::vector<int> row_marks(n_samples, 0);
    int row_mark = 0;
    std::vector<double> oob_improvements;
    oob_improvements.reserve(n_rounds);
    int actual_rounds = 0;
    double previous_oob_loss = NA_REAL;
    std::vector<int> train_rows;
    train_rows.reserve(sample_size);
    std::vector<unsigned char> is_train(n_samples);
    std::vector<ScenicTreeNode> nodes;
    nodes.reserve((1 << (max_depth + 1)) - 1);
    std::vector<int> feature_pool;
    feature_pool.reserve(features.size());
    std::vector<int> constant_features(features.size());

    for (int round = 0; round < n_rounds; ++round) {
      {
        ScenicScopeTimer sample_timer(
          profile == nullptr ? nullptr : &profile->sample_seconds
        );
        scenic_sample_mask_fill(n_samples, sample_size, rng, is_train, train_rows);
      }

      nodes.clear();
      std::fill(tree_importance.begin(), tree_importance.end(), 0.0);
      std::uint32_t split_seed = static_cast<std::uint32_t>(
        scenic_numpy_randint(rng, 0, SCENIC_RAND_R_MAX)
      );
      feature_pool.assign(features.begin(), features.end());
      std::fill(constant_features.begin(), constant_features.end(), 0);
      {
        ScenicScopeTimer build_timer(
          profile == nullptr ? nullptr : &profile->build_tree_seconds
        );
        scenic_build_tree(
          expr, residual, train_rows, feature_values, 0, static_cast<int>(train_rows.size()),
          feature_pool, constant_features, 0, 0, max_depth, mtry, split_seed,
          feature_orders, feature_order_values, feature_order_ranks,
          row_marks, row_mark, nodes, tree_importance, profile);
      }
      if (nodes.empty()) break;
      for (std::size_t ri = 0; ri < features.size(); ++ri) {
        const int reg = features[ri];
        importance[reg] += tree_importance[reg] / static_cast<double>(train_rows.size());
      }
      if (profile != nullptr) {
        ++profile->rounds;
        ++profile->trees;
      }

      double initial_oob_loss = 0.0;
      double new_oob_loss = 0.0;
      int n_oob = 0;
      if (sample_size < n_samples) {
        ScenicScopeTimer oob_timer(
          profile == nullptr ? nullptr : &profile->oob_predict_seconds
        );
        for (int i = 0; i < n_samples; ++i) {
          if (is_train[i]) continue;
          const double old_resid = residual[i];
          const double pred = scenic_predict_tree(expr, nodes, i);
          const double new_resid = old_resid - learning_rate * pred;
          initial_oob_loss += old_resid * old_resid;
          new_oob_loss += new_resid * new_resid;
          ++n_oob;
        }
      }
      {
        ScenicScopeTimer update_timer(
          profile == nullptr ? nullptr : &profile->update_predict_seconds
        );
        for (int i = 0; i < n_samples; ++i) {
          residual[i] -= learning_rate * scenic_predict_tree(expr, nodes, i);
        }
      }
      ++actual_rounds;

      if (n_oob > 0 && early_stop_window_length > 0) {
        initial_oob_loss /= static_cast<double>(n_oob);
        new_oob_loss /= static_cast<double>(n_oob);
        const double baseline_oob_loss = R_finite(previous_oob_loss) ?
          previous_oob_loss : initial_oob_loss;
        oob_improvements.push_back(baseline_oob_loss - new_oob_loss);
        previous_oob_loss = new_oob_loss;
        if (static_cast<int>(oob_improvements.size()) >= early_stop_window_length) {
          double window_sum = 0.0;
          for (
            int wi = static_cast<int>(oob_improvements.size()) - early_stop_window_length;
            wi < static_cast<int>(oob_improvements.size());
            ++wi
          ) {
            window_sum += oob_improvements[wi];
          }
          if (window_sum / static_cast<double>(early_stop_window_length) < 0.0) break;
        }
      }
    }

    double total_importance = 0.0;
    for (std::size_t ri = 0; ri < features.size(); ++ri) {
      total_importance += importance[features[ri]];
    }
    if (total_importance <= 0.0) continue;

    {
      ScenicScopeTimer edge_timer(
        profile == nullptr ? nullptr : &profile->edge_seconds
      );
      std::vector<ScenicEdge> target_edges;
      target_edges.reserve(features.size());
      const double arboreto_importance_scale = subsample < 1.0 ?
        static_cast<double>(actual_rounds) : 1.0;
      for (std::size_t ri = 0; ri < features.size(); ++ri) {
        const int reg = features[ri];
        if (importance[reg] > 0.0) {
          target_edges.push_back(ScenicEdge{
            reg + 1,
            target + 1,
            importance[reg] * arboreto_importance_scale / total_importance
          });
        }
      }
      std::sort(target_edges.begin(), target_edges.end(), scenic_edge_before);
      const int take = max_edges_per_target > 0 ?
        std::min<int>(max_edges_per_target, target_edges.size()) :
        static_cast<int>(target_edges.size());
      for (int i = 0; i < take; ++i) edges.push_back(target_edges[i]);
    }
    if (check_interrupt) Rcpp::checkUserInterrupt();
  }
}

template <typename Expr>
static DataFrame scenic_grnboost_tree_impl(
    const Expr& expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true,
    ScenicGrnProfile* profile = nullptr,
    int cores = 1) {
  std::chrono::steady_clock::time_point prepare_start;
  if (profile != nullptr) prepare_start = std::chrono::steady_clock::now();
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  if (n_samples < 2 || n_genes < 1) {
    stop("Expression matrix must contain at least two samples and one gene.");
  }
  if (regulator_idx.size() == 0 || target_idx.size() == 0) {
    stop("At least one regulator and one target are required.");
  }
  n_rounds = std::max(1, n_rounds);
  max_depth = std::max(1, max_depth);
  early_stop_window_length = std::max(0, early_stop_window_length);
  if (!R_finite(learning_rate) || learning_rate <= 0.0) learning_rate = 0.01;
  if (!R_finite(max_features) || max_features <= 0.0) max_features = 0.1;
  if (!R_finite(subsample) || subsample <= 0.0 || subsample > 1.0) subsample = 0.9;

  std::vector<int> regs;
  regs.reserve(regulator_idx.size());
  for (int i = 0; i < regulator_idx.size(); ++i) {
    const int idx = regulator_idx[i] - 1;
    if (idx >= 0 && idx < n_genes) regs.push_back(idx);
  }
  std::vector<int> targets;
  targets.reserve(target_idx.size());
  for (int i = 0; i < target_idx.size(); ++i) {
    const int idx = target_idx[i] - 1;
    if (idx >= 0 && idx < n_genes) targets.push_back(idx);
  }
  if (regs.empty() || targets.empty()) {
    stop("Regulator and target indices must refer to expression matrix columns.");
  }

  std::vector<double> means(n_genes, 0.0);
  for (int g = 0; g < n_genes; ++g) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
      double x = scenic_tree_feature_value(expr, i, g);
      sum += x;
    }
    means[g] = sum / static_cast<double>(n_samples);
  }

  std::vector<int> all_rows(n_samples);
  std::iota(all_rows.begin(), all_rows.end(), 0);
  scenic_profile_add_elapsed(profile, &ScenicGrnProfile::prepare_seconds, prepare_start);

  std::vector<ScenicEdge> edges;
  std::vector<std::vector<int> > feature_orders;
  std::vector<std::vector<double> > feature_order_values;
  std::vector<std::vector<int> > feature_order_ranks;
  {
    ScenicScopeTimer order_timer(
      profile == nullptr ? nullptr : &profile->feature_order_seconds
    );
    feature_orders = scenic_feature_orders(expr, regs, &feature_order_values, &feature_order_ranks);
  }

  const int workers = profile == nullptr ?
    scenic_grn_worker_count(cores, static_cast<int>(targets.size())) : 1;
  if (workers <= 1) {
    scenic_grnboost_run_targets(
      expr, regs, targets, 0, targets.size(), means, feature_orders,
      feature_order_values, feature_order_ranks, n_rounds, learning_rate, max_edges_per_target,
      max_depth, max_features, subsample, early_stop_window_length, random_seed,
      exclude_self, profile, true, edges);
  } else {
    std::vector<std::vector<ScenicEdge> > worker_edges(workers);
    std::vector<std::exception_ptr> errors(workers);
    std::vector<std::thread> pool;
    pool.reserve(workers);
    std::atomic<std::size_t> next_target(0);
    const std::size_t block = std::max<std::size_t>(
      1,
      targets.size() / static_cast<std::size_t>(workers * 8)
    );
    for (int worker = 0; worker < workers; ++worker) {
      pool.emplace_back([&, worker]() {
        try {
          while (true) {
            const std::size_t begin = next_target.fetch_add(block);
            if (begin >= targets.size()) break;
            const std::size_t end = std::min(targets.size(), begin + block);
            scenic_grnboost_run_targets(
              expr, regs, targets, begin, end, means, feature_orders,
              feature_order_values, feature_order_ranks, n_rounds, learning_rate, max_edges_per_target,
              max_depth, max_features, subsample, early_stop_window_length,
              random_seed, exclude_self, nullptr, false, worker_edges[worker]);
          }
        } catch (...) {
          errors[worker] = std::current_exception();
        }
      });
    }
    for (std::thread& worker : pool) {
      worker.join();
    }
    for (std::size_t i = 0; i < errors.size(); ++i) {
      if (errors[i]) std::rethrow_exception(errors[i]);
    }
    std::size_t total_edges = 0;
    for (std::size_t i = 0; i < worker_edges.size(); ++i) {
      total_edges += worker_edges[i].size();
    }
    edges.reserve(total_edges);
    for (std::size_t i = 0; i < worker_edges.size(); ++i) {
      edges.insert(edges.end(), worker_edges[i].begin(), worker_edges[i].end());
    }
  }

  {
    ScenicScopeTimer edge_timer(
      profile == nullptr ? nullptr : &profile->edge_seconds
    );
    std::sort(edges.begin(), edges.end(), scenic_edge_before);
  }
  IntegerVector regulator(edges.size());
  IntegerVector target(edges.size());
  NumericVector importance_out(edges.size());
  {
    ScenicScopeTimer edge_timer(
      profile == nullptr ? nullptr : &profile->edge_seconds
    );
    for (std::size_t i = 0; i < edges.size(); ++i) {
      regulator[i] = edges[i].regulator;
      target[i] = edges[i].target;
      importance_out[i] = edges[i].importance;
    }
  }

  return DataFrame::create(
    _["regulator"] = regulator,
    _["target"] = target,
    _["importance"] = importance_out,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
DataFrame grnboost_tree(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  return scenic_grnboost_tree_impl(
    expr, regulator_idx, target_idx, n_rounds, learning_rate, max_edges_per_target,
    max_depth, max_features, subsample, early_stop_window_length, random_seed,
    exclude_self, nullptr);
}

// [[Rcpp::export]]
DataFrame grnboost_tree_parallel(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true,
    int cores = 1) {
  return scenic_grnboost_tree_impl(
    expr, regulator_idx, target_idx, n_rounds, learning_rate, max_edges_per_target,
    max_depth, max_features, subsample, early_stop_window_length, random_seed,
    exclude_self, nullptr, cores);
}

// [[Rcpp::export]]
DataFrame grnboost_tree_sparse(
    S4 expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  ScenicSparseDgcMatrix sparse_expr(expr);
  return scenic_grnboost_tree_impl(
    sparse_expr, regulator_idx, target_idx, n_rounds, learning_rate, max_edges_per_target,
    max_depth, max_features, subsample, early_stop_window_length, random_seed,
    exclude_self, nullptr);
}

// [[Rcpp::export]]
DataFrame grnboost_tree_sparse_parallel(
    S4 expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true,
    int cores = 1) {
  ScenicSparseDgcMatrix sparse_expr(expr);
  return scenic_grnboost_tree_impl(
    sparse_expr, regulator_idx, target_idx, n_rounds, learning_rate, max_edges_per_target,
    max_depth, max_features, subsample, early_stop_window_length, random_seed,
    exclude_self, nullptr, cores);
}

// [[Rcpp::export]]
List grnboost_tree_profile(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    IntegerVector target_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_edges_per_target = 0,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  ScenicGrnProfile profile;
  DataFrame adjacency = scenic_grnboost_tree_impl(
    expr, regulator_idx, target_idx, n_rounds, learning_rate, max_edges_per_target,
    max_depth, max_features, subsample, early_stop_window_length, random_seed,
    exclude_self, &profile);
  NumericVector seconds = NumericVector::create(
    _["prepare"] = profile.prepare_seconds,
    _["feature_order"] = profile.feature_order_seconds,
    _["target_setup"] = profile.target_setup_seconds,
    _["sample"] = profile.sample_seconds,
    _["build_tree"] = profile.build_tree_seconds,
    _["best_split"] = profile.best_split_seconds,
    _["oob_predict"] = profile.oob_predict_seconds,
    _["update_predict"] = profile.update_predict_seconds,
    _["edge_output"] = profile.edge_seconds
  );
  IntegerVector counts = IntegerVector::create(
    _["targets"] = profile.targets,
    _["rounds"] = profile.rounds,
    _["trees"] = profile.trees,
    _["nodes"] = profile.nodes,
    _["best_split_calls"] = profile.best_split_calls,
    _["split_feature_visits"] = profile.split_feature_visits,
    _["split_threshold_visits"] = profile.split_threshold_visits
  );
  return List::create(
    _["adjacency"] = adjacency,
    _["seconds"] = seconds,
    _["counts"] = counts
  );
}

// [[Rcpp::export]]
DataFrame grnboost_tree_round_trace(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    int target_idx,
    IntegerVector trace_regulator_idx,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  n_rounds = std::max(1, n_rounds);
  max_depth = std::max(1, max_depth);
  early_stop_window_length = std::max(0, early_stop_window_length);
  if (!R_finite(learning_rate) || learning_rate <= 0.0) learning_rate = 0.01;
  if (!R_finite(max_features) || max_features <= 0.0) max_features = 0.1;
  if (!R_finite(subsample) || subsample <= 0.0 || subsample > 1.0) subsample = 0.9;

  std::vector<int> regs;
  regs.reserve(regulator_idx.size());
  for (int i = 0; i < regulator_idx.size(); ++i) {
    const int idx = regulator_idx[i] - 1;
    if (idx >= 0 && idx < n_genes) regs.push_back(idx);
  }
  const int target = target_idx - 1;
  if (target < 0 || target >= n_genes || regs.empty()) {
    stop("Trace target and regulator indices must refer to expression matrix columns.");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) {
    stop("No candidate regulators remain for the trace target.");
  }

  std::vector<int> trace_regs;
  trace_regs.reserve(trace_regulator_idx.size());
  for (int i = 0; i < trace_regulator_idx.size(); ++i) {
    const int idx = trace_regulator_idx[i] - 1;
    if (std::find(features.begin(), features.end(), idx) != features.end()) {
      trace_regs.push_back(idx);
    }
  }
  if (trace_regs.empty()) {
    stop("No trace regulators remain after filtering.");
  }

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    residual[i] = y - target_mean;
  }

  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  const int mtry = std::min<int>(
    features.size(),
    std::max(1, static_cast<int>(std::floor(max_features * features.size())))
  );
  std::vector<double> importance(n_genes, 0.0);
  std::vector<double> tree_importance(n_genes, 0.0);
  std::vector<double> feature_values(n_samples);
  std::vector<std::vector<double> > feature_order_values;
  std::vector<std::vector<int> > feature_order_ranks;
  std::vector<std::vector<int> > feature_orders = scenic_feature_orders(
    expr, features, &feature_order_values, &feature_order_ranks);
  std::vector<int> row_marks(n_samples, 0);
  int row_mark = 0;
  std::vector<double> oob_improvements;
  oob_improvements.reserve(n_rounds);
  double previous_oob_loss = NA_REAL;
  std::mt19937 rng(static_cast<unsigned int>(random_seed));

  std::vector<int> round_out;
  std::vector<int> regulator_out;
  std::vector<double> importance_out;
  std::vector<double> raw_importance_out;
  std::vector<double> tree_importance_out;
  std::vector<double> total_importance_out;
  std::vector<int> rank_out;
  std::vector<double> residual_sse_out;
  std::vector<double> oob_improvement_out;

  int actual_rounds = 0;
  for (int round = 0; round < n_rounds; ++round) {
    std::vector<unsigned char> is_train = scenic_sample_mask(n_samples, sample_size, rng);
    std::vector<int> train_rows;
    train_rows.reserve(sample_size);
    for (int i = 0; i < n_samples; ++i) {
      if (is_train[i]) train_rows.push_back(i);
    }

    std::fill(tree_importance.begin(), tree_importance.end(), 0.0);
    std::uint32_t split_seed = static_cast<std::uint32_t>(
      scenic_numpy_randint(rng, 0, SCENIC_RAND_R_MAX)
    );
    std::vector<int> feature_pool(features);
    std::vector<int> constant_features(features.size());
    std::vector<ScenicTreeNode> nodes;
    nodes.reserve((1 << (max_depth + 1)) - 1);
    scenic_build_tree(
      expr, residual, train_rows, feature_values, 0, static_cast<int>(train_rows.size()),
      feature_pool, constant_features, 0, 0, max_depth, mtry, split_seed,
      feature_orders, feature_order_values, feature_order_ranks,
      row_marks, row_mark, nodes, tree_importance, nullptr);
    if (nodes.empty()) break;
    for (std::size_t ri = 0; ri < features.size(); ++ri) {
      const int reg = features[ri];
      importance[reg] += tree_importance[reg] / static_cast<double>(train_rows.size());
    }

    double initial_oob_loss = 0.0;
    double new_oob_loss = 0.0;
    int n_oob = 0;
    if (sample_size < n_samples) {
      for (int i = 0; i < n_samples; ++i) {
        if (is_train[i]) continue;
        const double old_resid = residual[i];
        const double pred = scenic_predict_tree(expr, nodes, i);
        const double new_resid = old_resid - learning_rate * pred;
        initial_oob_loss += old_resid * old_resid;
        new_oob_loss += new_resid * new_resid;
        ++n_oob;
      }
    }

    for (int i = 0; i < n_samples; ++i) {
      residual[i] -= learning_rate * scenic_predict_tree(expr, nodes, i);
    }
    ++actual_rounds;

    double oob_improvement = NA_REAL;
    if (n_oob > 0 && early_stop_window_length > 0) {
      initial_oob_loss /= static_cast<double>(n_oob);
      new_oob_loss /= static_cast<double>(n_oob);
      const double baseline_oob_loss = R_finite(previous_oob_loss) ?
        previous_oob_loss : initial_oob_loss;
      oob_improvement = baseline_oob_loss - new_oob_loss;
      oob_improvements.push_back(oob_improvement);
      previous_oob_loss = new_oob_loss;
    }

    double total_importance = 0.0;
    for (std::size_t ri = 0; ri < features.size(); ++ri) {
      total_importance += importance[features[ri]];
    }
    const double arboreto_importance_scale = subsample < 1.0 ?
      static_cast<double>(actual_rounds) : 1.0;
    double residual_sse = 0.0;
    for (int i = 0; i < n_samples; ++i) residual_sse += residual[i] * residual[i];

    std::vector<std::pair<double, int> > ranked;
    ranked.reserve(features.size());
    for (std::size_t ri = 0; ri < features.size(); ++ri) {
      const int reg = features[ri];
      const double val = total_importance > 0.0 ?
        importance[reg] * arboreto_importance_scale / total_importance : 0.0;
      ranked.push_back(std::make_pair(-val, reg));
    }
    std::sort(ranked.begin(), ranked.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

    for (std::size_t ti = 0; ti < trace_regs.size(); ++ti) {
      const int reg = trace_regs[ti];
      int rank = NA_INTEGER;
      for (std::size_t ri = 0; ri < ranked.size(); ++ri) {
        if (ranked[ri].second == reg) {
          rank = static_cast<int>(ri) + 1;
          break;
        }
      }
      round_out.push_back(actual_rounds);
      regulator_out.push_back(reg + 1);
      importance_out.push_back(total_importance > 0.0 ?
        importance[reg] * arboreto_importance_scale / total_importance : 0.0);
      raw_importance_out.push_back(importance[reg]);
      tree_importance_out.push_back(tree_importance[reg] / static_cast<double>(train_rows.size()));
      total_importance_out.push_back(total_importance);
      rank_out.push_back(rank);
      residual_sse_out.push_back(residual_sse);
      oob_improvement_out.push_back(oob_improvement);
    }

    if (n_oob > 0 && early_stop_window_length > 0 &&
        static_cast<int>(oob_improvements.size()) >= early_stop_window_length) {
      double window_sum = 0.0;
      for (
        int wi = static_cast<int>(oob_improvements.size()) - early_stop_window_length;
        wi < static_cast<int>(oob_improvements.size());
        ++wi
      ) {
        window_sum += oob_improvements[wi];
      }
      if (window_sum / static_cast<double>(early_stop_window_length) < 0.0) break;
    }
  }

  return DataFrame::create(
    _["round"] = round_out,
    _["regulator"] = regulator_out,
    _["importance"] = importance_out,
    _["raw_importance"] = raw_importance_out,
    _["tree_importance"] = tree_importance_out,
    _["total_importance"] = total_importance_out,
    _["rank"] = rank_out,
    _["residual_sse"] = residual_sse_out,
    _["oob_improvement"] = oob_improvement_out,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
DataFrame grnboost_tree_round_nodes(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    int target_idx,
    int trace_round,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  n_rounds = std::max(1, n_rounds);
  trace_round = std::max(1, std::min(trace_round, n_rounds));
  max_depth = std::max(1, max_depth);
  early_stop_window_length = std::max(0, early_stop_window_length);
  if (!R_finite(learning_rate) || learning_rate <= 0.0) learning_rate = 0.01;
  if (!R_finite(max_features) || max_features <= 0.0) max_features = 0.1;
  if (!R_finite(subsample) || subsample <= 0.0 || subsample > 1.0) subsample = 0.9;

  std::vector<int> regs;
  regs.reserve(regulator_idx.size());
  for (int i = 0; i < regulator_idx.size(); ++i) {
    const int idx = regulator_idx[i] - 1;
    if (idx >= 0 && idx < n_genes) regs.push_back(idx);
  }
  const int target = target_idx - 1;
  if (target < 0 || target >= n_genes || regs.empty()) {
    stop("Trace target and regulator indices must refer to expression matrix columns.");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) stop("No candidate regulators remain for the trace target.");

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    residual[i] = y - target_mean;
  }

  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  const int mtry = std::min<int>(
    features.size(),
    std::max(1, static_cast<int>(std::floor(max_features * features.size())))
  );
  std::vector<double> tree_importance(n_genes, 0.0);
  std::vector<double> feature_values(n_samples);
  std::vector<std::vector<double> > feature_order_values;
  std::vector<std::vector<int> > feature_order_ranks;
  std::vector<std::vector<int> > feature_orders = scenic_feature_orders(
    expr, features, &feature_order_values, &feature_order_ranks);
  std::vector<int> row_marks(n_samples, 0);
  int row_mark = 0;
  std::vector<double> oob_improvements;
  oob_improvements.reserve(n_rounds);
  double previous_oob_loss = NA_REAL;
  std::mt19937 rng(static_cast<unsigned int>(random_seed));

  std::vector<ScenicTreeNode> trace_nodes;
  std::vector<int> trace_train_rows;

  int actual_rounds = 0;
  for (int round = 0; round < n_rounds; ++round) {
    std::vector<unsigned char> is_train = scenic_sample_mask(n_samples, sample_size, rng);
    std::vector<int> train_rows;
    train_rows.reserve(sample_size);
    for (int i = 0; i < n_samples; ++i) {
      if (is_train[i]) train_rows.push_back(i);
    }

    std::fill(tree_importance.begin(), tree_importance.end(), 0.0);
    std::uint32_t split_seed = static_cast<std::uint32_t>(
      scenic_numpy_randint(rng, 0, SCENIC_RAND_R_MAX)
    );
    std::vector<int> feature_pool(features);
    std::vector<int> constant_features(features.size());
    std::vector<ScenicTreeNode> nodes;
    nodes.reserve((1 << (max_depth + 1)) - 1);
    scenic_build_tree(
      expr, residual, train_rows, feature_values, 0, static_cast<int>(train_rows.size()),
      feature_pool, constant_features, 0, 0, max_depth, mtry, split_seed,
      feature_orders, feature_order_values, feature_order_ranks,
      row_marks, row_mark, nodes, tree_importance, nullptr);
    if (nodes.empty()) break;
    ++actual_rounds;

    if (actual_rounds == trace_round) {
      trace_nodes = nodes;
      trace_train_rows = train_rows;
      break;
    }

    double initial_oob_loss = 0.0;
    double new_oob_loss = 0.0;
    int n_oob = 0;
    if (sample_size < n_samples) {
      for (int i = 0; i < n_samples; ++i) {
        if (is_train[i]) continue;
        const double old_resid = residual[i];
        const double pred = scenic_predict_tree(expr, nodes, i);
        const double new_resid = old_resid - learning_rate * pred;
        initial_oob_loss += old_resid * old_resid;
        new_oob_loss += new_resid * new_resid;
        ++n_oob;
      }
    }
    for (int i = 0; i < n_samples; ++i) {
      residual[i] -= learning_rate * scenic_predict_tree(expr, nodes, i);
    }
    if (n_oob > 0 && early_stop_window_length > 0) {
      initial_oob_loss /= static_cast<double>(n_oob);
      new_oob_loss /= static_cast<double>(n_oob);
      const double baseline_oob_loss = R_finite(previous_oob_loss) ?
        previous_oob_loss : initial_oob_loss;
      oob_improvements.push_back(baseline_oob_loss - new_oob_loss);
      previous_oob_loss = new_oob_loss;
      if (static_cast<int>(oob_improvements.size()) >= early_stop_window_length) {
        double window_sum = 0.0;
        for (
          int wi = static_cast<int>(oob_improvements.size()) - early_stop_window_length;
          wi < static_cast<int>(oob_improvements.size());
          ++wi
        ) {
          window_sum += oob_improvements[wi];
        }
        if (window_sum / static_cast<double>(early_stop_window_length) < 0.0) break;
      }
    }
  }

  std::vector<int> node_out;
  std::vector<int> feature_out;
  std::vector<int> left_out;
  std::vector<int> right_out;
  std::vector<int> n_rows_out;
  std::vector<double> threshold_out;
  std::vector<double> value_out;
  std::vector<double> sse_out;
  std::vector<double> gain_out;
  node_out.reserve(trace_nodes.size());
  feature_out.reserve(trace_nodes.size());
  left_out.reserve(trace_nodes.size());
  right_out.reserve(trace_nodes.size());
  n_rows_out.reserve(trace_nodes.size());
  threshold_out.reserve(trace_nodes.size());
  value_out.reserve(trace_nodes.size());
  sse_out.reserve(trace_nodes.size());
  gain_out.reserve(trace_nodes.size());

  for (std::size_t i = 0; i < trace_nodes.size(); ++i) {
    const ScenicTreeNode& node = trace_nodes[i];
    node_out.push_back(static_cast<int>(i) + 1);
    feature_out.push_back(node.feature >= 0 ? node.feature + 1 : NA_INTEGER);
    left_out.push_back(node.left >= 0 ? node.left + 1 : NA_INTEGER);
    right_out.push_back(node.right >= 0 ? node.right + 1 : NA_INTEGER);
    n_rows_out.push_back(node.n_rows);
    threshold_out.push_back(node.threshold);
    value_out.push_back(node.value);
    sse_out.push_back(node.sse);
    gain_out.push_back(node.importance_gain);
  }

  return DataFrame::create(
    _["round"] = trace_round,
    _["node"] = node_out,
    _["feature"] = feature_out,
    _["threshold"] = threshold_out,
    _["value"] = value_out,
    _["left"] = left_out,
    _["right"] = right_out,
    _["n_rows"] = n_rows_out,
    _["sse"] = sse_out,
    _["importance_gain"] = gain_out,
    _["train_rows"] = static_cast<int>(trace_train_rows.size()),
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
DataFrame grnboost_tree_node_candidates(
    NumericMatrix expr,
    IntegerVector regulator_idx,
    int target_idx,
    int trace_round,
    int trace_node,
    int n_rounds = 5000,
    double learning_rate = 0.01,
    int max_depth = 3,
    double max_features = 0.1,
    double subsample = 0.9,
    int early_stop_window_length = 25,
    int random_seed = 1234,
    bool exclude_self = true) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  n_rounds = std::max(1, n_rounds);
  trace_round = std::max(1, std::min(trace_round, n_rounds));
  trace_node = std::max(1, trace_node);
  max_depth = std::max(1, max_depth);
  early_stop_window_length = std::max(0, early_stop_window_length);
  if (!R_finite(learning_rate) || learning_rate <= 0.0) learning_rate = 0.01;
  if (!R_finite(max_features) || max_features <= 0.0) max_features = 0.1;
  if (!R_finite(subsample) || subsample <= 0.0 || subsample > 1.0) subsample = 0.9;

  std::vector<int> regs;
  regs.reserve(regulator_idx.size());
  for (int i = 0; i < regulator_idx.size(); ++i) {
    const int idx = regulator_idx[i] - 1;
    if (idx >= 0 && idx < n_genes) regs.push_back(idx);
  }
  const int target = target_idx - 1;
  if (target < 0 || target >= n_genes || regs.empty()) {
    stop("Trace target and regulator indices must refer to expression matrix columns.");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) stop("No candidate regulators remain for the trace target.");

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    residual[i] = y - target_mean;
  }

  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  const int mtry = std::min<int>(
    features.size(),
    std::max(1, static_cast<int>(std::floor(max_features * features.size())))
  );
  std::vector<double> tree_importance(n_genes, 0.0);
  std::vector<double> feature_values(n_samples);
  std::vector<std::vector<double> > feature_order_values;
  std::vector<std::vector<int> > feature_order_ranks;
  std::vector<std::vector<int> > feature_orders = scenic_feature_orders(
    expr, features, &feature_order_values, &feature_order_ranks);
  std::vector<int> row_marks(n_samples, 0);
  int row_mark = 0;
  std::vector<double> oob_improvements;
  oob_improvements.reserve(n_rounds);
  double previous_oob_loss = NA_REAL;
  std::mt19937 rng(static_cast<unsigned int>(random_seed));
  ScenicCandidateTrace candidate_trace;
  candidate_trace.trace_node = trace_node;

  int actual_rounds = 0;
  for (int round = 0; round < n_rounds; ++round) {
    std::vector<unsigned char> is_train = scenic_sample_mask(n_samples, sample_size, rng);
    std::vector<int> train_rows;
    train_rows.reserve(sample_size);
    for (int i = 0; i < n_samples; ++i) {
      if (is_train[i]) train_rows.push_back(i);
    }

    std::fill(tree_importance.begin(), tree_importance.end(), 0.0);
    std::uint32_t split_seed = static_cast<std::uint32_t>(
      scenic_numpy_randint(rng, 0, SCENIC_RAND_R_MAX)
    );
    std::vector<int> feature_pool(features);
    std::vector<int> constant_features(features.size());
    std::vector<ScenicTreeNode> nodes;
    nodes.reserve((1 << (max_depth + 1)) - 1);
    ScenicCandidateTrace* trace_ptr = (round + 1 == trace_round) ? &candidate_trace : nullptr;
    scenic_build_tree(
      expr, residual, train_rows, feature_values, 0, static_cast<int>(train_rows.size()),
      feature_pool, constant_features, 0, 0, max_depth, mtry, split_seed,
      feature_orders, feature_order_values, feature_order_ranks,
      row_marks, row_mark, nodes, tree_importance, nullptr, trace_ptr);
    if (nodes.empty()) break;
    ++actual_rounds;
    if (actual_rounds == trace_round) break;

    double initial_oob_loss = 0.0;
    double new_oob_loss = 0.0;
    int n_oob = 0;
    if (sample_size < n_samples) {
      for (int i = 0; i < n_samples; ++i) {
        if (is_train[i]) continue;
        const double old_resid = residual[i];
        const double pred = scenic_predict_tree(expr, nodes, i);
        const double new_resid = old_resid - learning_rate * pred;
        initial_oob_loss += old_resid * old_resid;
        new_oob_loss += new_resid * new_resid;
        ++n_oob;
      }
    }
    for (int i = 0; i < n_samples; ++i) {
      residual[i] -= learning_rate * scenic_predict_tree(expr, nodes, i);
    }
    if (n_oob > 0 && early_stop_window_length > 0) {
      initial_oob_loss /= static_cast<double>(n_oob);
      new_oob_loss /= static_cast<double>(n_oob);
      const double baseline_oob_loss = R_finite(previous_oob_loss) ?
        previous_oob_loss : initial_oob_loss;
      oob_improvements.push_back(baseline_oob_loss - new_oob_loss);
      previous_oob_loss = new_oob_loss;
      if (static_cast<int>(oob_improvements.size()) >= early_stop_window_length) {
        double window_sum = 0.0;
        for (
          int wi = static_cast<int>(oob_improvements.size()) - early_stop_window_length;
          wi < static_cast<int>(oob_improvements.size());
          ++wi
        ) {
          window_sum += oob_improvements[wi];
        }
        if (window_sum / static_cast<double>(early_stop_window_length) < 0.0) break;
      }
    }
  }

  IntegerVector feature_out(candidate_trace.feature.size());
  for (std::size_t i = 0; i < candidate_trace.feature.size(); ++i) {
    feature_out[i] = candidate_trace.feature[i];
  }
  return DataFrame::create(
    _["round"] = trace_round,
    _["node"] = trace_node,
    _["visit_order"] = candidate_trace.visit_order,
    _["feature"] = feature_out,
    _["proxy_gain"] = candidate_trace.proxy_gain,
    _["importance_gain"] = candidate_trace.importance_gain,
    _["threshold"] = candidate_trace.threshold,
    _["n_left"] = candidate_trace.n_left,
    _["n_right"] = candidate_trace.n_right,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
NumericMatrix scenic_ctx_recovery(IntegerMatrix ranks, NumericVector weights, int rank_threshold) {
  const int n_features = ranks.nrow();
  const int n_targets = ranks.ncol();
  rank_threshold = std::max(1, rank_threshold);

  NumericMatrix out(n_features, rank_threshold);
  std::vector<double> counts(rank_threshold);

  for (int feature = 0; feature < n_features; ++feature) {
    std::fill(counts.begin(), counts.end(), 0.0);
    for (int target = 0; target < n_targets; ++target) {
      const int rank = ranks(feature, target);
      if (rank == NA_INTEGER || rank < 0 || rank >= rank_threshold) continue;
      counts[rank] += weights[target];
    }

    double cumulative = 0.0;
    for (int rank = 0; rank < rank_threshold; ++rank) {
      cumulative += counts[rank];
      out(feature, rank) = cumulative;
    }
  }

  return out;
}

// [[Rcpp::export]]
List scenic_ctx_auc_avg2sd(
    IntegerMatrix ranks,
    int total_genes,
    int rank_threshold,
    int rank_cutoff) {
  const int n_features = ranks.nrow();
  const int n_targets = ranks.ncol();
  rank_threshold = std::min(std::max(1, rank_threshold), std::max(1, total_genes - 1));
  if (rank_cutoff < 1 || rank_cutoff > rank_threshold) {
    stop("auc_threshold and rank_threshold produce an invalid cisTarget rank cutoff");
  }

  NumericVector auc_out(n_features);
  NumericVector avg2sd_out(rank_threshold);
  std::vector<double> recovery_sum_diff(rank_threshold + 1, 0.0);
  std::vector<double> recovery_sumsq_diff(rank_threshold + 1, 0.0);
  std::vector<int> selected;
  const double max_auc = static_cast<double>(rank_cutoff + 1) *
    static_cast<double>(n_targets);

  for (int feature = 0; feature < n_features; ++feature) {
    selected.clear();
    for (int target = 0; target < n_targets; ++target) {
      const int rank = ranks(feature, target);
      if (rank == NA_INTEGER) continue;
      if (rank >= 0 && rank < rank_threshold) selected.push_back(rank);
    }

    if (!selected.empty()) {
      std::sort(selected.begin(), selected.end());
      std::size_t auc_n = 0;
      while (auc_n < selected.size() && selected[auc_n] < rank_cutoff) ++auc_n;
      double auc = 0.0;
      if (auc_n > 0) {
        int previous = selected[0];
        for (std::size_t i = 1; i < auc_n; ++i) {
          auc += static_cast<double>(selected[i] - previous) *
            static_cast<double>(i);
          previous = selected[i];
        }
        auc += static_cast<double>(rank_cutoff - previous) *
          static_cast<double>(auc_n);
        auc_out[feature] = auc / max_auc;
      }

      std::size_t pos = 0;
      double cumulative = 0.0;
      while (pos < selected.size()) {
        const int rank = selected[pos];
        std::size_t next_pos = pos + 1;
        while (next_pos < selected.size() && selected[next_pos] == rank) {
          ++next_pos;
        }
        cumulative += static_cast<double>(next_pos - pos);
        const int next_rank = (next_pos < selected.size()) ?
          selected[next_pos] : rank_threshold;
        if (next_rank > rank) {
          const double cumulative_sq = cumulative * cumulative;
          recovery_sum_diff[rank] += cumulative;
          recovery_sum_diff[next_rank] -= cumulative;
          recovery_sumsq_diff[rank] += cumulative_sq;
          recovery_sumsq_diff[next_rank] -= cumulative_sq;
        }
        pos = next_pos;
      }
    }
  }

  if (n_features > 0) {
    const double denom = static_cast<double>(n_features);
    double recovery_sum = 0.0;
    double recovery_sumsq = 0.0;
    for (int rank = 0; rank < rank_threshold; ++rank) {
      recovery_sum += recovery_sum_diff[rank];
      recovery_sumsq += recovery_sumsq_diff[rank];
      const double mean = recovery_sum / denom;
      const double variance = std::max(
        0.0,
        recovery_sumsq / denom - mean * mean
      );
      avg2sd_out[rank] = mean + 2.0 * std::sqrt(variance);
    }
  }

  return List::create(
    _["auc"] = auc_out,
    _["avg2sd_recovery"] = avg2sd_out
  );
}

// [[Rcpp::export]]
List scenic_ctx_auc_nes(
    IntegerMatrix ranks,
    double nes_threshold,
    int rank_cutoff) {
  const int n_features = ranks.nrow();
  const int n_targets = ranks.ncol();
  rank_cutoff = std::max(1, rank_cutoff);

  std::vector<double> auc_vec(n_features);
  std::vector<double> counts(rank_cutoff);
  std::vector<double> rcc_row(rank_cutoff);

  for (int feature = 0; feature < n_features; ++feature) {
    std::fill(counts.begin(), counts.end(), 0.0);
    for (int target = 0; target < n_targets; ++target) {
      const int rank = ranks(feature, target);
      if (rank == NA_INTEGER || rank < 0 || rank >= rank_cutoff) continue;
      counts[rank] += 1.0;
    }
    double cumulative = 0.0;
    double prefix_sum = 0.0;
    for (int rank = 0; rank < rank_cutoff; ++rank) {
      cumulative += counts[rank];
      rcc_row[rank] = cumulative;
      prefix_sum += cumulative;
    }
    auc_vec[feature] = prefix_sum / (static_cast<double>(rank_cutoff + 1) * n_targets);
  }

  double auc_mean = 0.0;
  for (int i = 0; i < n_features; ++i) auc_mean += auc_vec[i];
  auc_mean /= static_cast<double>(n_features);

  double auc_var = 0.0;
  for (int i = 0; i < n_features; ++i) {
    double d = auc_vec[i] - auc_mean;
    auc_var += d * d;
  }
  auc_var /= static_cast<double>(n_features);
  double auc_sd = std::sqrt(auc_var);

  std::vector<int> enriched_idx;
  if (auc_sd > 0.0) {
    for (int i = 0; i < n_features; ++i) {
      double nes = (auc_vec[i] - auc_mean) / auc_sd;
      if (R_finite(nes) && nes >= nes_threshold) {
        enriched_idx.push_back(i);
      }
    }
  }

  const int n_enriched = static_cast<int>(enriched_idx.size());
  IntegerVector enriched_out(n_enriched);
  NumericVector auc_out(n_enriched);
  NumericVector nes_out(n_enriched);
  NumericMatrix rcc_out(n_enriched, rank_cutoff);
  NumericMatrix avg2stdrcc_out(1, rank_cutoff);

  for (int e = 0; e < n_enriched; ++e) {
    enriched_out[e] = enriched_idx[e] + 1;
    auc_out[e] = auc_vec[enriched_idx[e]];
    nes_out[e] = (auc_vec[enriched_idx[e]] - auc_mean) / auc_sd;
  }

  std::vector<double> col_mean(rank_cutoff, 0.0);
  std::vector<double> col_sq_mean(rank_cutoff, 0.0);
  for (int feature = 0; feature < n_features; ++feature) {
    std::fill(counts.begin(), counts.end(), 0.0);
    for (int target = 0; target < n_targets; ++target) {
      const int rank = ranks(feature, target);
      if (rank == NA_INTEGER || rank < 0 || rank >= rank_cutoff) continue;
      counts[rank] += 1.0;
    }
    double cumulative = 0.0;
    for (int rank = 0; rank < rank_cutoff; ++rank) {
      cumulative += counts[rank];
      col_mean[rank] += cumulative;
      col_sq_mean[rank] += cumulative * cumulative;
    }
  }
  for (int rank = 0; rank < rank_cutoff; ++rank) {
    double cm = col_mean[rank] / static_cast<double>(n_features);
    double cs = col_sq_mean[rank] / static_cast<double>(n_features);
    avg2stdrcc_out(0, rank) = cm + 2.0 * std::sqrt(std::max(0.0, cs - cm * cm));
  }

  for (int e = 0; e < n_enriched; ++e) {
    int fi = enriched_idx[e];
    std::fill(counts.begin(), counts.end(), 0.0);
    for (int target = 0; target < n_targets; ++target) {
      const int rank = ranks(fi, target);
      if (rank == NA_INTEGER || rank < 0 || rank >= rank_cutoff) continue;
      counts[rank] += 1.0;
    }
    double cumulative = 0.0;
    for (int rank = 0; rank < rank_cutoff; ++rank) {
      cumulative += counts[rank];
      rcc_out(e, rank) = cumulative;
    }
  }

  return List::create(
    _["enriched_idx"] = enriched_out,
    _["auc"] = auc_out,
    _["nes"] = nes_out,
    _["rcc"] = rcc_out,
    _["avg2stdrcc"] = avg2stdrcc_out
  );
}

// [[Rcpp::export]]
NumericVector scenicplus_region_gene_cor(
    NumericMatrix atac_log,
    NumericMatrix rna_log,
    IntegerVector region_idx,
    IntegerVector gene_idx) {
  const int n_cells = atac_log.ncol();
  if (rna_log.ncol() != n_cells) {
    stop("ATAC and RNA matrices must have the same number of cells.");
  }
  if (region_idx.size() != gene_idx.size()) {
    stop("Region and gene index vectors must have the same length.");
  }
  NumericVector scores(region_idx.size());
  for (int pair = 0; pair < region_idx.size(); ++pair) {
    const int region = region_idx[pair] - 1;
    const int gene = gene_idx[pair] - 1;
    if (region < 0 || region >= atac_log.nrow() || gene < 0 || gene >= rna_log.nrow()) {
      scores[pair] = NA_REAL;
      continue;
    }
    double region_sum = 0.0;
    double gene_sum = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      double x = atac_log(region, cell);
      double y = rna_log(gene, cell);
      if (!R_finite(x)) x = 0.0;
      if (!R_finite(y)) y = 0.0;
      region_sum += x;
      gene_sum += y;
    }
    const double region_mean = region_sum / static_cast<double>(n_cells);
    const double gene_mean = gene_sum / static_cast<double>(n_cells);
    double numerator = 0.0;
    double region_ss = 0.0;
    double gene_ss = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      double x = atac_log(region, cell);
      double y = rna_log(gene, cell);
      if (!R_finite(x)) x = 0.0;
      if (!R_finite(y)) y = 0.0;
      x -= region_mean;
      y -= gene_mean;
      numerator += x * y;
      region_ss += x * x;
      gene_ss += y * y;
    }
    const double denom = std::sqrt(region_ss * gene_ss);
    scores[pair] = denom > 0.0 ? numerator / denom : NA_REAL;
  }
  return scores;
}

// [[Rcpp::export]]
DataFrame scenicplus_triplets_cpp(
    CharacterVector tf_gene_tf,
    CharacterVector tf_gene_target,
    NumericVector tf_gene_importance,
    CharacterVector region_gene_region,
    CharacterVector region_gene_gene,
    NumericVector region_gene_score,
    CharacterVector tf_region_tf,
    CharacterVector tf_region_region,
    NumericVector tf_region_score) {
  std::unordered_map<std::string, double> tf_gene_weight;
  tf_gene_weight.reserve(tf_gene_tf.size() * 2 + 1);
  for (int i = 0; i < tf_gene_tf.size(); ++i) {
    if (!R_finite(tf_gene_importance[i])) continue;
    const std::string key = as<std::string>(tf_gene_tf[i]) + "\r" + as<std::string>(tf_gene_target[i]);
    if (tf_gene_weight.find(key) == tf_gene_weight.end()) {
      tf_gene_weight[key] = tf_gene_importance[i];
    }
  }

  std::unordered_map<std::string, std::vector<int> > region_gene_rows;
  region_gene_rows.reserve(region_gene_region.size() * 2 + 1);
  for (int i = 0; i < region_gene_region.size(); ++i) {
    region_gene_rows[as<std::string>(region_gene_region[i])].push_back(i);
  }

  std::vector<std::string> out_tf;
  std::vector<std::string> out_region;
  std::vector<std::string> out_gene;
  std::vector<double> out_score;

  for (int tr = 0; tr < tf_region_tf.size(); ++tr) {
    if (!R_finite(tf_region_score[tr])) continue;
    const std::string tf = as<std::string>(tf_region_tf[tr]);
    const std::string region = as<std::string>(tf_region_region[tr]);
    std::unordered_map<std::string, std::vector<int> >::const_iterator rg_it = region_gene_rows.find(region);
    if (rg_it == region_gene_rows.end()) continue;
    const std::vector<int>& rows = rg_it->second;
    for (std::size_t j = 0; j < rows.size(); ++j) {
      const int rg = rows[j];
      if (!R_finite(region_gene_score[rg])) continue;
      const std::string gene = as<std::string>(region_gene_gene[rg]);
      const std::string key = tf + "\r" + gene;
      std::unordered_map<std::string, double>::const_iterator weight_it = tf_gene_weight.find(key);
      if (weight_it == tf_gene_weight.end()) continue;
      out_tf.push_back(tf);
      out_region.push_back(region);
      out_gene.push_back(gene);
      out_score.push_back(std::fabs(region_gene_score[rg]));
    }
  }

  return DataFrame::create(
    _["TF"] = out_tf,
    _["region"] = out_region,
    _["gene"] = out_gene,
    _["score"] = out_score,
    _["stringsAsFactors"] = false
  );
}
