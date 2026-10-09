// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <thisutils/log_message.h>
#include "thread_utils.h"
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
#ifdef _OPENMP
#include <omp.h>
#endif

// Keep split and impurity arithmetic as separate IEEE operations on Clang.
// Fused multiply-add can change equal-gain feature attribution across compilers.
#ifdef __clang__
#pragma clang fp contract(off)
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector scenic_edge_correlation_cpp(
    NumericMatrix expr,
    IntegerVector tf_index,
    IntegerVector target_index,
    int n_threads = 0) {
  const int n_cells = expr.nrow();
  const int n_genes = expr.ncol();
  const int n_edges = tf_index.size();
  if (n_edges != target_index.size()) {
    thisutils::log_message("SCENIC edge index vectors must have the same length", "error");
  }
  NumericVector rho(n_edges, NA_REAL);
  if (n_cells < 2) return rho;

  for (int edge = 0; edge < n_edges; ++edge) {
    const int tf = tf_index[edge] - 1;
    const int target = target_index[edge] - 1;
    if (tf < 0 || tf >= n_genes || target < 0 || target >= n_genes) {
      thisutils::log_message("SCENIC edge index is out of bounds", "error");
    }
  }

  std::vector<double> means(n_genes, 0.0);
  std::vector<double> sumsq(n_genes, 0.0);
  std::vector<bool> valid(n_genes, true);
  const int threads = omp_thread_count(n_threads, std::max(n_genes, n_edges));

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int gene = 0; gene < n_genes; ++gene) {
    double sum = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      const double value = expr(cell, gene);
      if (!std::isfinite(value)) {
        valid[gene] = false;
        break;
      }
      sum += value;
    }
    if (valid[gene]) {
      means[gene] = sum / static_cast<double>(n_cells);
    }
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int gene = 0; gene < n_genes; ++gene) {
    if (!valid[gene]) continue;
    double sumsq_gene = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      const double centered = expr(cell, gene) - means[gene];
      sumsq_gene += centered * centered;
    }
    sumsq[gene] = sumsq_gene;
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int edge = 0; edge < n_edges; ++edge) {
    const int tf = tf_index[edge] - 1;
    const int target = target_index[edge] - 1;
    if (!valid[tf] || !valid[target]) continue;
    double cross = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      const double left = expr(cell, tf) - means[tf];
      const double right = expr(cell, target) - means[target];
      cross += left * right;
    }
    const double scale = std::sqrt(sumsq[tf] * sumsq[target]);
    if (std::isfinite(scale) && scale > 0.0) {
      rho[edge] = cross / scale;
    }
  }
  return rho;
}

static const double SCENIC_FEATURE_THRESHOLD = 1e-7;
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
      thisutils::log_message("expr must be a dgCMatrix.", "error");
    }
    if (dim.size() != 2) {
      thisutils::log_message("expr must have two dimensions.", "error");
    }
    n_rows = dim[0];
    n_cols = dim[1];
    if (p.size() != n_cols + 1) {
      thisutils::log_message("Invalid dgCMatrix column pointer.", "error");
    }
    if (i.size() != x.size()) {
      thisutils::log_message("Invalid dgCMatrix index/value slots.", "error");
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
  return R_finite(candidate) && candidate > current;
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

static double scenic_target_value(const NumericMatrix& expr, int row, int feature) {
  double x = expr.begin()[row + feature * expr.nrow()];
  if (!R_finite(x)) x = 0.0;
  return x;
}

static double scenic_tree_feature_value(const NumericMatrix& expr, int row, int feature) {
  return static_cast<double>(static_cast<float>(scenic_target_value(expr, row, feature)));
}

// Read-only columns keep per-gene candidate matrices inside the existing GBM kernel.
struct ScenicColumnView {
  int rows;
  std::vector<const double*> columns;
  int nrow() const { return rows; }
  int ncol() const { return static_cast<int>(columns.size()); }
};

static double scenic_target_value(const ScenicColumnView& expr, int row, int feature) {
  const double value = expr.columns[feature][row];
  return R_finite(value) ? value : 0.0;
}

static double scenic_tree_feature_value(const ScenicColumnView& expr, int row, int feature) {
  return static_cast<double>(static_cast<float>(scenic_target_value(expr, row, feature)));
}

static double scenic_target_value(const ScenicSparseDgcMatrix& expr, int row, int feature) {
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
  return value;
}

static double scenic_tree_feature_value(const ScenicSparseDgcMatrix& expr, int row, int feature) {
  return static_cast<double>(static_cast<float>(scenic_target_value(expr, row, feature)));
}

struct ScenicSparseTreeView {
  const ScenicSparseDgcMatrix& matrix;
  const std::vector<std::vector<double> >& values;
  const std::vector<std::vector<int> >& ranks;
  int nrow() const { return matrix.nrow(); }
};

template <typename Expr>
static const Expr& scenic_tree_view(
    const Expr& expr,
    const std::vector<std::vector<double> >&,
    const std::vector<std::vector<int> >&) {
  return expr;
}

static ScenicSparseTreeView scenic_tree_view(
    const ScenicSparseDgcMatrix& expr,
    const std::vector<std::vector<double> >& values,
    const std::vector<std::vector<int> >& ranks) {
  return ScenicSparseTreeView{expr, values, ranks};
}

static double scenic_tree_feature_value(const ScenicSparseTreeView& view, int row, int feature) {
  if (feature >= 0 && feature < static_cast<int>(view.ranks.size()) &&
      static_cast<int>(view.ranks[feature].size()) == view.nrow()) {
    const int rank = view.ranks[feature][row];
    return rank < 0 ? 0.0 : view.values[feature][rank];
  }
  return scenic_tree_feature_value(view.matrix, row, feature);
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

// Three-way introsort preserves the reference dense splitter's tie partition order.
static void scenic_sample_swap(std::vector<double>& values, std::vector<int>& samples, int a, int b) {
  std::swap(values[a], values[b]);
  std::swap(samples[a], samples[b]);
}

static void scenic_sample_sift(std::vector<double>& values, std::vector<int>& samples, int offset, int root, int end) {
  while (true) {
    const int child = root * 2 + 1;
    int largest = root;
    if (child < end && values[offset + largest] < values[offset + child]) largest = child;
    if (child + 1 < end && values[offset + largest] < values[offset + child + 1]) largest = child + 1;
    if (largest == root) return;
    scenic_sample_swap(values, samples, offset + root, offset + largest);
    root = largest;
  }
}

static void scenic_sample_sort(std::vector<double>& values, std::vector<int>& samples, int start, int count, int depth) {
  while (count > 1) {
    if (depth <= 0) {
      for (int root = (count - 2) / 2; root >= 0; --root) scenic_sample_sift(values, samples, start, root, count);
      for (int end = count - 1; end > 0; --end) {
        scenic_sample_swap(values, samples, start, start + end);
        scenic_sample_sift(values, samples, start, 0, end);
      }
      return;
    }
    --depth;
    const double a = values[start], b = values[start + count / 2], c = values[start + count - 1];
    const double pivot = a < b ? (b < c ? b : (a < c ? c : a)) : (b < c ? (a < c ? a : c) : b);
    int i = 0, left = 0, right = count;
    while (i < right) {
      if (values[start + i] < pivot) {
        scenic_sample_swap(values, samples, start + i, start + left);
        ++i; ++left;
      } else if (values[start + i] > pivot) {
        --right;
        scenic_sample_swap(values, samples, start + i, start + right);
      } else ++i;
    }
    scenic_sample_sort(values, samples, start, left, depth);
    start += right;
    count -= right;
  }
}

template <typename Expr>
static void scenic_sort_samples_by_feature(
    const Expr& expr, std::vector<int>& samples, std::vector<double>& feature_values,
    int start, int end, int feature) {
  for (int i = start; i < end; ++i) feature_values[i] = scenic_tree_feature_value(expr, samples[i], feature);
  const int count = end - start;
  if (count > 1) scenic_sample_sort(feature_values, samples, start, count, 2 * static_cast<int>(std::log(count)));
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
  // Everything before p is on the left; only the meeting sample remains.
  return p + (p < end && scenic_tree_feature_value(expr, samples[p], feature) <= threshold ? 1 : 0);
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
    std::vector<std::pair<double, int> > nonzero;
    nonzero.reserve(n_samples);
    for (int i = 0; i < n_samples; ++i) {
      const double value = scenic_tree_feature_value(expr, i, g);
      if (!R_finite(value) || static_cast<float>(value) == 0.0f) {
        continue;
      }
      nonzero.push_back(std::make_pair(
        static_cast<double>(static_cast<float>(value)),
        i
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
    if (order_values != nullptr) {
      (*order_values)[g] = values;
    }
    if (order_ranks != nullptr) {
      std::vector<int>& ranks = (*order_ranks)[g];
      ranks.assign(expr.nrow(), -1);
      for (std::size_t i = 0; i < order.size(); ++i) ranks[order[i]] = static_cast<int>(i);
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
      if (!R_finite(value) || static_cast<float>(value) == 0.0f) continue;
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
    if (order_ranks != nullptr) {
      std::vector<int>& ranks = (*order_ranks)[g];
      ranks.assign(expr.nrow(), -1);
      for (std::size_t i = 0; i < order.size(); ++i) ranks[order[i]] = static_cast<int>(i);
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
    double parent_sse,
    double parent_sum,
    ScenicGrnProfile* profile,
    ScenicCandidateTrace* candidate_trace = nullptr) {
  if (profile != nullptr) ++profile->best_split_calls;
  ScenicSplit best{-1, 0.0, -std::numeric_limits<double>::infinity(), 0.0};
  if (end - start < 2 || feature_pool.empty()) return best;

  if (parent_sse <= 0.0) return best;
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

    double feature_best_proxy = -std::numeric_limits<double>::infinity();
    double feature_best_importance = 0.0;
    double feature_best_threshold = 0.0;
    int feature_best_left = 0;
    int feature_best_right = 0;
    double left_sum = 0.0;
    int n_left = 0;
    // Reuse cached float32 values, while sorting the current node sample order.
    const bool cached = feature >= 0 && feature < static_cast<int>(feature_order_ranks.size()) &&
      static_cast<int>(feature_order_ranks[feature].size()) == expr.nrow();
    if (cached) {
      const std::vector<int>& ranks = feature_order_ranks[feature];
      const std::vector<double>& values = feature_order_values[feature];
      for (int i = start; i < end; ++i) {
        const int rank = ranks[samples[i]];
        feature_values[i] = rank < 0 ? 0.0 : values[rank];
      }
      const int count = end - start;
      if (count > 1) scenic_sample_sort(feature_values, samples, start, count, 2 * static_cast<int>(std::log(count)));
    } else {
      scenic_sort_samples_by_feature(expr, samples, feature_values, start, end, feature);
    }
    const int selected_n = end - start;
    const double min_x = feature_values[start], max_x = feature_values[end - 1];
    // Match sequential prefix updates and reverse suffix updates, including
    // floating-point tie decisions. A relative gain tolerance changes TF attribution.
    int previous_pos = start;
    int pos = start;
    while (pos < end) {
      while (pos + 1 < end && static_cast<float>(feature_values[pos + 1]) <=
             static_cast<float>(static_cast<float>(feature_values[pos]) + static_cast<float>(SCENIC_FEATURE_THRESHOLD))) ++pos;
      const int boundary = pos;
      ++pos;
      if (pos == end) continue;
      n_left = pos - start;
      const int n_right = end - pos;
      if (profile != nullptr) ++profile->split_threshold_visits;
      if (pos - previous_pos <= end - pos) {
        for (int row = previous_pos; row < pos; ++row) left_sum += residual[samples[row]];
      } else {
        left_sum = parent_sum;
        for (int row = end - 1; row >= pos; --row) left_sum -= residual[samples[row]];
      }
      previous_pos = pos;
      const double right_sum = parent_sum - left_sum;
      const double diff = static_cast<double>(n_right) * left_sum - static_cast<double>(n_left) * right_sum;
      const double proxy_gain = diff * diff / (static_cast<double>(n_left) * static_cast<double>(n_right));
      const bool feature_improves = scenic_proxy_improves(proxy_gain, feature_best_proxy);
      const bool global_improves = scenic_proxy_improves(proxy_gain, best.proxy_gain);
      if (feature_improves || global_improves) {
        const double importance_gain = proxy_gain / static_cast<double>(n_left + n_right);
        double split_threshold = feature_values[boundary] / 2.0 + feature_values[pos] / 2.0;
        if (split_threshold == feature_values[pos] || !R_finite(split_threshold)) split_threshold = feature_values[boundary];
        if (feature_improves) {
          feature_best_proxy = proxy_gain;
          feature_best_importance = importance_gain;
          feature_best_threshold = split_threshold;
          feature_best_left = n_left;
          feature_best_right = n_right;
        }
        if (global_improves) {
          best.feature = feature;
          best.threshold = split_threshold;
          best.proxy_gain = proxy_gain;
          best.importance_gain = importance_gain;
        }
      }
    }
    if (selected_n < 2 || static_cast<float>(max_x) <= static_cast<float>(static_cast<float>(min_x) + static_cast<float>(SCENIC_FEATURE_THRESHOLD))) {
      std::swap(feature_pool[f_j], feature_pool[n_total_constants]);
      ++n_found_constants;
      ++n_total_constants;
      continue;
    }
    --f_i;
    std::swap(feature_pool[f_i], feature_pool[f_j]);
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
    ScenicCandidateTrace* candidate_trace = nullptr,
    double inherited_impurity = std::numeric_limits<double>::quiet_NaN()) {
  ScenicTreeNode node;
  node.feature = -1;
  node.threshold = 0.0;
  node.left = -1;
  node.right = -1;
  node.n_rows = end - start;
  double node_sum = 0.0;
  double node_sum_sq = 0.0;
  scenic_node_sse(samples, start, end, residual, &node_sum, &node_sum_sq);
  const double count = static_cast<double>(end - start);
  node.value = end > start ? node_sum / count : 0.0;
  const double impurity = std::isnan(inherited_impurity) ?
    node_sum_sq / count - (node_sum / count) * (node_sum / count) : inherited_impurity;
  node.sse = count * impurity;
  node.importance_gain = 0.0;
  const int node_index = static_cast<int>(nodes.size());
  nodes.push_back(node);
  if (profile != nullptr) ++profile->nodes;

  if (depth >= max_depth || end - start < 2 || feature_pool.empty()) return node_index;
  const bool near_zero_leaf =
    node.sse / static_cast<double>(end - start) <= std::numeric_limits<double>::epsilon();
  if (near_zero_leaf) {
    return node_index;
  }
  int n_total_constants = n_known_constants;
  const bool old_collect = candidate_trace != nullptr && candidate_trace->collect;
  if (candidate_trace != nullptr) {
    candidate_trace->collect = (node_index + 1 == candidate_trace->trace_node);
  }
  const ScenicSplit split = scenic_best_split(
    expr, residual, samples, feature_values, start, end, feature_pool, constant_features,
    n_known_constants, n_total_constants, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks,
    row_marks, row_mark, node.sse, node_sum, profile, candidate_trace);
  if (candidate_trace != nullptr) {
    candidate_trace->collect = old_collect;
  }
  if (split.feature < 0) return node_index;

  const auto& partition_view = scenic_tree_view(expr, feature_order_values, feature_order_ranks);
  const int split_pos = scenic_partition_samples_final(
    partition_view, samples, start, end, split.feature, split.threshold);
  if (split_pos <= start || split_pos >= end) return node_index;

  double left_sum = 0.0, left_sq = 0.0;
  if (split_pos - start <= end - split_pos) {
    for (int i = start; i < split_pos; ++i) left_sum += residual[samples[i]];
  } else {
    left_sum = node_sum;
    for (int i = end - 1; i >= split_pos; --i) left_sum -= residual[samples[i]];
  }
  for (int i = start; i < split_pos; ++i) {
    const double value = residual[samples[i]];
    left_sq += value * value;
  }
  const double nl = split_pos - start, nr = end - split_pos;
  const double right_sum = node_sum - left_sum;
  const double left_impurity = left_sq / nl - (left_sum / nl) * (left_sum / nl);
  const double right_impurity = (node_sum_sq - left_sq) / nr - (right_sum / nr) * (right_sum / nr);
  const double importance_gain = count * impurity - nl * left_impurity - nr * right_impurity;
  importance[split.feature] += importance_gain;
  const int left = scenic_build_tree(
    expr, residual, samples, feature_values, start, split_pos, feature_pool, constant_features,
    n_total_constants, depth + 1, max_depth, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks, row_marks, row_mark, nodes, importance, profile,
    candidate_trace, left_impurity);
  const int right = scenic_build_tree(
    expr, residual, samples, feature_values, split_pos, end, feature_pool, constant_features,
    n_total_constants, depth + 1, max_depth, mtry, split_seed, feature_orders,
    feature_order_values, feature_order_ranks, row_marks, row_mark, nodes, importance, profile,
    candidate_trace, right_impurity);
  nodes[node_index].feature = split.feature;
  nodes[node_index].threshold = split.threshold;
  nodes[node_index].left = left;
  nodes[node_index].right = right;
  nodes[node_index].importance_gain = importance_gain;
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

// Accumulate fitted predictions, then derive residuals from the original response.
// The same update is used by production fitting and diagnostic tree traces.
template <typename Expr>
static void scenic_boost_update(
    const Expr& expr, const std::vector<ScenicTreeNode>& nodes,
    const std::vector<double>& response, std::vector<double>& prediction,
    std::vector<double>& residual, const std::vector<unsigned char>& is_train,
    int sample_size, double learning_rate, double& initial_loss, double& new_loss,
    int& n_oob) {
  for (int i = 0; i < expr.nrow(); ++i) {
    const double old_residual = residual[i];
    prediction[i] += learning_rate * scenic_predict_tree(expr, nodes, i);
    residual[i] = response[i] - prediction[i];
    if (sample_size < expr.nrow() && !is_train[i]) {
      initial_loss += old_residual * old_residual;
      new_loss += residual[i] * residual[i];
      ++n_oob;
    }
  }
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
  const auto& prediction_view = scenic_tree_view(expr, feature_order_values, feature_order_ranks);
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  std::vector<double> residual(n_samples);
  std::vector<double> response(n_samples), prediction(n_samples);
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
    const double first_target = scenic_target_value(expr, 0, target);
    bool constant_target = true;
    for (int i = 0; i < n_samples; ++i) {
      double y = scenic_target_value(expr, i, target);
      if (y != first_target) constant_target = false;
      response[i] = y;
      prediction[i] = means[target];
      residual[i] = y - means[target];
    }
    // A constant response cannot yield a split or a nonzero importance. The
    // OOB monitor otherwise repeats all rounds because its improvement is zero.
    if (constant_target) {
      scenic_profile_add_elapsed(
        profile, &ScenicGrnProfile::target_setup_seconds, target_setup_start);
      continue;
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
      {
        ScenicScopeTimer update_timer(
          profile == nullptr ? nullptr : &profile->update_predict_seconds
        );
        scenic_boost_update(prediction_view, nodes, response, prediction, residual, is_train,
                            sample_size, learning_rate, initial_oob_loss, new_oob_loss, n_oob);
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

// Index introsort keeps the deterministic partition order for equal group keys.
template <typename Value>
static std::vector<int> scenic_argsort(const std::vector<Value>& values) {
  std::vector<int> order(values.size());
  std::iota(order.begin(), order.end(), 0);
  if (order.empty()) return order;
  struct Partition { int first, last, depth; };
  int depth = 0;
  for (std::size_t size = values.size(); size > 1; size >>= 1) depth += 2;
  std::vector<Partition> stack{{0, static_cast<int>(order.size()) - 1, depth}};
  auto less = [&](int left, int right) { return values[left] < values[right]; };
  while (!stack.empty()) {
    auto part = stack.back();
    stack.pop_back();
    int first = part.first, last = part.last;
    depth = part.depth;
    if (depth < 0) {
      std::make_heap(order.begin() + first, order.begin() + last + 1, less);
      std::sort_heap(order.begin() + first, order.begin() + last + 1, less);
      continue;
    }
    while (last - first > 15) {
      const int middle = first + (last - first) / 2;
      if (less(order[middle], order[first])) std::swap(order[middle], order[first]);
      if (less(order[last], order[middle])) std::swap(order[last], order[middle]);
      if (less(order[middle], order[first])) std::swap(order[middle], order[first]);
      const Value pivot = values[order[middle]];
      std::swap(order[middle], order[last - 1]);
      int left = first, right = last - 1;
      while (true) {
        do { ++left; } while (values[order[left]] < pivot);
        do { --right; } while (pivot < values[order[right]]);
        if (left >= right) break;
        std::swap(order[left], order[right]);
      }
      std::swap(order[left], order[last - 1]);
      --depth;
      if (left - first < last - left) {
        stack.push_back({left + 1, last, depth});
        last = left - 1;
      } else {
        stack.push_back({first, left - 1, depth});
        first = left + 1;
      }
    }
    for (int i = first + 1; i <= last; ++i) {
      const int value = order[i];
      int j = i;
      while (j > first && less(value, order[j - 1])) {
        order[j] = order[j - 1];
        --j;
      }
      order[j] = value;
    }
  }
  return order;
}

// [[Rcpp::export]]
IntegerVector scenicplus_order_cpp(SEXP values, bool decreasing = false) {
  std::vector<int> order;
  if (TYPEOF(values) == STRSXP) {
    CharacterVector groups(values);
    std::vector<std::string> keys;
    keys.reserve(groups.size());
    for (String group : groups) {
      if (group == NA_STRING) stop("Order keys must not contain missing names");
      keys.push_back(group.get_cstring());
    }
    if (decreasing) std::reverse(keys.begin(), keys.end());
    order = scenic_argsort(keys);
  } else if (TYPEOF(values) == REALSXP || TYPEOF(values) == INTSXP) {
    NumericVector input = as<NumericVector>(values);
    std::vector<double> keys(input.begin(), input.end());
    for (double key : keys) if (!std::isfinite(key)) stop("Order keys must be finite");
    if (decreasing) std::reverse(keys.begin(), keys.end());
    order = scenic_argsort(keys);
  } else {
    stop("Order keys must be numeric or character");
  }
  // Reverse both the input and output for descending order, preserving the
  // reference tie permutation rather than treating equal scores as stable.
  if (decreasing) {
    std::reverse(order.begin(), order.end());
    for (int& index : order) index = order.size() - 1 - index;
  }
  for (int& index : order) ++index;
  return wrap(order);
}

template <typename Function>
static void scenic_parallel_tasks(int n_tasks, int cores, Function run) {
  const int workers = worker_count(cores, n_tasks);
  // Bounded batches allow interrupts without calling R from worker threads.
  for (int first = 0; first < n_tasks; first += workers * 4) {
    const int last = std::min(n_tasks, first + workers * 4);
    std::atomic<int> next(first);
    std::vector<std::exception_ptr> errors(workers);
    std::vector<std::thread> pool;
    if (workers == 1) {
      for (int task = first; task < last; ++task) run(task);
    } else {
      for (int worker = 0; worker < workers; ++worker) pool.emplace_back([&, worker]() {
        try {
          while (true) {
            const int task = next.fetch_add(1);
            if (task >= last) break;
            run(task);
          }
        } catch (...) { errors[worker] = std::current_exception(); }
      });
      for (auto& worker : pool) worker.join();
      for (const auto& error : errors) if (error) std::rethrow_exception(error);
    }
    Rcpp::checkUserInterrupt();
  }
}

// PCG XSL-RR with a four-word seed mixer. Two limbs avoid compiler-specific
// 128-bit types; this stream is used only for SCENIC+ ranking ties.
class ScenicRankRng {
  uint64_t high = 0, low = 0, inc_high, inc_low;
  uint32_t spare = 0;
  bool has_spare = false;
  static uint64_t product_high(uint64_t a, uint64_t b) {
    const uint64_t a0 = static_cast<uint32_t>(a), a1 = a >> 32;
    const uint64_t b0 = static_cast<uint32_t>(b), b1 = b >> 32;
    const uint64_t carry = ((a0 * b0) >> 32) + static_cast<uint32_t>(a1 * b0) +
      static_cast<uint32_t>(a0 * b1);
    return a1 * b1 + (a1 * b0 >> 32) + (a0 * b1 >> 32) + (carry >> 32);
  }
  void step() {
    constexpr uint64_t m0 = 0x4385df649fccf645ULL, m1 = 0x2360ed051fc65da4ULL;
    const uint64_t next_low = low * m0;
    high = high * m0 + low * m1 + product_high(low, m0) + inc_high;
    low = next_low + inc_low;
    high += low < next_low;
  }
  uint32_t next32() {
    if (has_spare) { has_spare = false; return spare; }
    step();
    const uint64_t value = high ^ low;
    const unsigned int rotation = high >> 58;
    const uint64_t rotated = (value >> rotation) | (value << ((64 - rotation) & 63));
    spare = rotated >> 32;
    has_spare = true;
    return static_cast<uint32_t>(rotated);
  }
public:
  explicit ScenicRankRng(uint32_t seed) {
    uint32_t constant = 0x43b0d7e5U;
    auto hash = [&](uint32_t value) {
      value ^= constant;
      constant *= 0x931e8875U;
      value *= constant;
      return value ^ (value >> 16);
    };
    uint32_t pool[4];
    for (int i = 0; i < 4; ++i) pool[i] = hash(i == 0 ? seed : 0U);
    for (int source = 0; source < 4; ++source) {
      for (int target = 0; target < 4; ++target) if (source != target) {
        const uint32_t mixed = 0xca01f9ddU * pool[target] - 0x4973f715U * hash(pool[source]);
        pool[target] = mixed ^ (mixed >> 16);
      }
    }
    uint32_t words[8];
    constant = 0x8b51f9ddU;
    for (int i = 0; i < 8; ++i) {
      uint32_t value = pool[i % 4] ^ constant;
      constant *= 0x58f38dedU;
      value *= constant;
      words[i] = value ^ (value >> 16);
    }
    uint64_t seed_words[4];
    for (int i = 0; i < 4; ++i) seed_words[i] = words[2 * i] | (uint64_t(words[2 * i + 1]) << 32);
    inc_high = (seed_words[2] << 1) | (seed_words[3] >> 63);
    inc_low = (seed_words[3] << 1) | 1U;
    step();
    const uint64_t previous = low;
    low += seed_words[1];
    high += seed_words[0] + (low < previous);
    step();
  }
  uint32_t interval(uint32_t maximum) {
    uint32_t mask = maximum;
    mask |= mask >> 1; mask |= mask >> 2; mask |= mask >> 4;
    mask |= mask >> 8; mask |= mask >> 16;
    uint32_t value;
    do { value = next32() & mask; } while (value > maximum);
    return value;
  }
};

// [[Rcpp::export]]
IntegerMatrix scenicplus_rank_cpp(NumericMatrix expression, int random_seed = 123, int cores = 1) {
  const int features = expression.nrow(), cells = expression.ncol();
  if (features < 1 || cells < 1 || random_seed < 0) stop("Rankings require nonempty expression and a nonnegative seed");
  const double* expr = expression.begin();
  for (R_xlen_t i = 0; i < expression.size(); ++i) if (!std::isfinite(expr[i])) stop("Rankings require finite expression");
  IntegerMatrix result(cells, features);
  int* ranks = result.begin();
  ScenicRankRng rng(static_cast<uint32_t>(random_seed));
  const int batch_size = worker_count(cores, cells) * 4;
  for (int first = 0; first < cells; first += batch_size) {
    const int count = std::min(batch_size, cells - first);
    std::vector<std::vector<int> > permutations(count, std::vector<int>(features));
    // Generate in cell order so scheduling and core count cannot change ties.
    for (auto& permutation : permutations) {
      std::iota(permutation.begin(), permutation.end(), 0);
      for (int i = features - 1; i > 0; --i) std::swap(permutation[i], permutation[rng.interval(i)]);
    }
    scenic_parallel_tasks(count, cores, [&](int task) {
      const int cell = first + task;
      const auto& permutation = permutations[task];
      std::vector<double> values(features);
      for (int feature = 0; feature < features; ++feature) values[feature] = -expr[cell * features + permutation[feature]];
      const auto order = scenic_argsort(values);
      for (int rank = 0; rank < features; ++rank) ranks[cell + cells * permutation[order[rank]]] = rank;
    });
  }
  return result;
}

// [[Rcpp::export]]
List scenicplus_region_gbm_cpp(
    NumericMatrix accessibility,
    NumericMatrix expression,
    List candidates,
    int n_rounds = 500,
    double max_features = 0.1,
    int random_seed = 666,
    int cores = 1) {
  const int n_cells = accessibility.nrow();
  const int n_tasks = expression.ncol();
  if (n_cells < 2 || expression.nrow() != n_cells || candidates.size() != n_tasks) {
    stop("Region GBM requires matching cells and one candidate list per gene");
  }
  std::vector<std::vector<int> > indices(n_tasks);
  for (int task = 0; task < n_tasks; ++task) {
    IntegerVector values = as<IntegerVector>(candidates[task]);
    for (int index : values) {
      if (index == NA_INTEGER || index < 1 || index > accessibility.ncol()) {
        stop("Region GBM candidate index is out of bounds");
      }
      indices[task].push_back(index - 1);
    }
    if (indices[task].empty()) stop("Region GBM candidate lists must not be empty");
  }
  n_rounds = std::max(1, n_rounds);
  if (!R_finite(max_features) || max_features <= 0.0) max_features = 0.1;
  const double* acc = accessibility.begin();
  const double* rna = expression.begin();
  std::vector<std::vector<double> > weights(n_tasks);
  auto run = [&](int task) {
    const int n_regions = indices[task].size();
    ScenicColumnView view{n_cells, {}};
    view.columns.reserve(n_regions + 1);
    for (int region : indices[task]) view.columns.push_back(acc + region * n_cells);
    view.columns.push_back(rna + task * n_cells);
    std::vector<int> regs(n_regions);
    std::iota(regs.begin(), regs.end(), 0);
    std::vector<double> means(n_regions + 1, 0.0);
    const double first_response = scenic_target_value(view, 0, n_regions);
    bool constant_response = true;
    for (int cell = 0; cell < n_cells; ++cell) {
      const double response = scenic_target_value(view, cell, n_regions);
      means[n_regions] += response;
      if (response != first_response) constant_response = false;
    }
    means[n_regions] /= static_cast<double>(n_cells);
    weights[task].assign(n_regions, 0.0);
    // A constant target has no edges; avoid preparing predictor caches for it.
    if (constant_response) return;
    std::vector<std::vector<double> > values;
    std::vector<std::vector<int> > ranks;
    const auto orders = scenic_feature_orders(view, regs, &values, &ranks);
    std::vector<int> targets{n_regions};
    std::vector<ScenicEdge> edges;
    scenic_grnboost_run_targets(view, regs, targets, 0, 1, means, orders,
      values, ranks, n_rounds, 0.01, 0, 3, max_features, 1.0, 25,
      random_seed, true, nullptr, false, edges);
    weights[task].assign(n_regions, 0.0);
    for (const ScenicEdge& edge : edges) weights[task][edge.regulator - 1] = edge.importance;
    const double total = std::accumulate(weights[task].begin(), weights[task].end(), 0.0);
    if (total > 0.0) for (double& value : weights[task]) value /= total;
  };
  scenic_parallel_tasks(n_tasks, cores, run);
  List out(n_tasks);
  for (int task = 0; task < n_tasks; ++task) out[task] = wrap(weights[task]);
  return out;
}

// [[Rcpp::export]]
List scenicplus_gsea_cpp(
    List rankings,
    List hits,
    int n_perm = 1000,
    int random_seed = 666,
    int cores = 1) {
  const int n_tasks = rankings.size();
  if (hits.size() != n_tasks || n_perm < 1) stop("GSEA requires matching modules and positive permutations");
  std::vector<std::vector<double> > weights(n_tasks);
  std::vector<std::vector<unsigned char> > tags(n_tasks);
  for (int task = 0; task < n_tasks; ++task) {
    NumericVector ranking = as<NumericVector>(rankings[task]);
    weights[task].assign(ranking.begin(), ranking.end());
    for (double& value : weights[task]) {
      if (!std::isfinite(value)) stop("GSEA rankings must be finite");
      value = std::fabs(value);
    }
    tags[task].assign(ranking.size(), 0);
    IntegerVector hit = as<IntegerVector>(hits[task]);
    for (int index : hit) {
      if (index == NA_INTEGER || index < 1 || index > ranking.size()) stop("GSEA hit index is out of bounds");
      tags[task][index - 1] = 1;
    }
  }
  std::vector<double> es(n_tasks, NA_REAL), nes(n_tasks, NA_REAL), pval(n_tasks, NA_REAL);
  std::vector<std::vector<int> > leading(n_tasks);
  auto run = [&](int task) {
    const auto& weight = weights[task];
    const auto& tag = tags[task];
    const int n = tag.size();
    const int n_hits = std::accumulate(tag.begin(), tag.end(), 0);
    if (n_hits == 0 || n_hits == n) return;
    auto score = [&](const std::vector<unsigned char>& mask, int* peak) {
      double total = 0.0;
      for (int i = 0; i < n; ++i) if (mask[i]) total += weight[i];
      if (total <= 0.0) return NA_REAL;
      const double miss = 1.0 / static_cast<double>(n - n_hits);
      const double norm = 1.0 / total;
      double running = 0.0, maximum = -std::numeric_limits<double>::infinity();
      double minimum = std::numeric_limits<double>::infinity();
      int max_index = 0, min_index = 0;
      for (int i = 0; i < n; ++i) {
        running += mask[i] ? weight[i] * norm : -miss;
        if (running > maximum) { maximum = running; max_index = i; }
        if (running < minimum) { minimum = running; min_index = i; }
      }
      const bool positive = std::fabs(maximum) > std::fabs(minimum);
      if (peak) *peak = positive ? max_index : min_index;
      return positive ? maximum : minimum;
    };
    int peak = 0;
    es[task] = score(tag, &peak);
    if (!std::isfinite(es[task])) return;
    std::mt19937 seed_rng(static_cast<unsigned int>(random_seed));
    std::mt19937 rng(static_cast<unsigned int>(scenic_numpy_randint(seed_rng, 0, 2147483647)));
    int same_sign = 0, extreme = 0;
    double null_sum = 0.0;
    for (int permutation = 0; permutation < n_perm; ++permutation) {
      auto shuffled = tag;
      for (int i = n - 1; i > 0; --i) {
        const int j = scenic_numpy_randint(rng, 0, i + 1);
        std::swap(shuffled[i], shuffled[j]);
      }
      const double value = score(shuffled, nullptr);
      if ((es[task] >= 0.0 && value >= 0.0) || (es[task] < 0.0 && value < 0.0)) {
        ++same_sign;
        null_sum += value;
        if ((es[task] >= 0.0 && value >= es[task]) || (es[task] < 0.0 && value < es[task])) ++extreme;
      }
    }
    if (same_sign > 0 && null_sum != 0.0) {
      nes[task] = es[task] / (std::fabs(null_sum) / same_sign);
      pval[task] = static_cast<double>(extreme) / same_sign;
      for (int i = 0; i < n; ++i) {
        if (tag[i] && ((es[task] > 0.0 && i <= peak) || (es[task] < 0.0 && i >= peak) || es[task] == 0.0)) {
          leading[task].push_back(i + 1);
        }
      }
    }
  };
  scenic_parallel_tasks(n_tasks, cores, run);
  List edge(n_tasks);
  for (int task = 0; task < n_tasks; ++task) edge[task] = wrap(leading[task]);
  return List::create(_["ES"] = wrap(es), _["NES"] = wrap(nes),
                      _["pval"] = wrap(pval), _["leading_edge"] = edge);
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
    thisutils::log_message("Expression matrix must contain at least two samples and one gene.", "error");
  }
  if (regulator_idx.size() == 0 || target_idx.size() == 0) {
    thisutils::log_message("At least one regulator and one target are required.", "error");
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
    thisutils::log_message("Regulator and target indices must refer to expression matrix columns.", "error");
  }

  std::vector<double> means(n_genes, 0.0);
  for (int g = 0; g < n_genes; ++g) {
    double sum = 0.0;
    for (int i = 0; i < n_samples; ++i) {
      double x = scenic_target_value(expr, i, g);
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
    worker_count(cores, static_cast<int>(targets.size())) : 1;
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
    thisutils::log_message("Trace target and regulator indices must refer to expression matrix columns.", "error");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) {
    thisutils::log_message("No candidate regulators remain for the trace target.", "error");
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
    thisutils::log_message("No trace regulators remain after filtering.", "error");
  }

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  std::vector<double> response(n_samples), prediction(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    response[i] = y;
    prediction[i] = target_mean;
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
  const auto& prediction_view = scenic_tree_view(expr, feature_order_values, feature_order_ranks);
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
    scenic_boost_update(prediction_view, nodes, response, prediction, residual, is_train,
                        sample_size, learning_rate, initial_oob_loss, new_oob_loss, n_oob);
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
    thisutils::log_message("Trace target and regulator indices must refer to expression matrix columns.", "error");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) thisutils::log_message("No candidate regulators remain for the trace target.", "error");

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  std::vector<double> response(n_samples), prediction(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    response[i] = y;
    prediction[i] = target_mean;
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
  const auto& prediction_view = scenic_tree_view(expr, feature_order_values, feature_order_ranks);
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
    scenic_boost_update(prediction_view, nodes, response, prediction, residual, is_train,
                        sample_size, learning_rate, initial_oob_loss, new_oob_loss, n_oob);
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
    thisutils::log_message("Trace target and regulator indices must refer to expression matrix columns.", "error");
  }

  std::vector<int> features;
  features.reserve(regs.size());
  for (std::size_t ri = 0; ri < regs.size(); ++ri) {
    if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
  }
  if (features.empty()) thisutils::log_message("No candidate regulators remain for the trace target.", "error");

  double target_mean = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    target_mean += y;
  }
  target_mean /= static_cast<double>(n_samples);

  std::vector<double> residual(n_samples);
  std::vector<double> response(n_samples), prediction(n_samples);
  for (int i = 0; i < n_samples; ++i) {
    double y = expr(i, target);
    if (!R_finite(y)) y = 0.0;
    response[i] = y;
    prediction[i] = target_mean;
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
  const auto& prediction_view = scenic_tree_view(expr, feature_order_values, feature_order_ranks);
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
    scenic_boost_update(prediction_view, nodes, response, prediction, residual, is_train,
                        sample_size, learning_rate, initial_oob_loss, new_oob_loss, n_oob);
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
  IntegerVector round_out(candidate_trace.feature.size(), trace_round);
  IntegerVector node_out(candidate_trace.feature.size(), trace_node);
  return DataFrame::create(
    _["round"] = round_out,
    _["node"] = node_out,
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
  if (weights.size() != n_targets) {
    thisutils::log_message(
      "SCENIC recovery weights must have one value per rank-matrix column",
      "error"
    );
  }
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
    thisutils::log_message("auc_threshold and rank_threshold produce an invalid cisTarget rank cutoff", "error");
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
      while (auc_n < selected.size() && selected[auc_n] <= rank_cutoff) ++auc_n;
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
NumericVector scenicplus_region_gene_cor(
    NumericMatrix atac_log,
    NumericMatrix rna_log,
    IntegerVector region_idx,
    IntegerVector gene_idx) {
  const int n_cells = atac_log.ncol();
  if (rna_log.ncol() != n_cells) {
    thisutils::log_message("ATAC and RNA matrices must have the same number of cells.", "error");
  }
  if (region_idx.size() != gene_idx.size()) {
    thisutils::log_message("Region and gene index vectors must have the same length.", "error");
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
