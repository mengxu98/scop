// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

static const double SCENIC_FEATURE_THRESHOLD = 1e-12;
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
};

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

static double scenic_node_mean(
    const std::vector<int>& rows,
    const std::vector<double>& residual) {
  if (rows.empty()) return 0.0;
  double sum = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) sum += residual[rows[i]];
  return sum / static_cast<double>(rows.size());
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

static std::vector<std::vector<int> > scenic_feature_orders(const NumericMatrix& expr) {
  const int n_samples = expr.nrow();
  const int n_genes = expr.ncol();
  std::vector<std::vector<int> > orders(n_genes);
  for (int g = 0; g < n_genes; ++g) {
    std::vector<int> order(n_samples);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      double xa = expr(a, g);
      double xb = expr(b, g);
      if (!R_finite(xa)) xa = 0.0;
      if (!R_finite(xb)) xb = 0.0;
      if (xa != xb) return xa < xb;
      return a < b;
    });
    orders[g] = order;
  }
  return orders;
}

static ScenicSplit scenic_best_split(
    const NumericMatrix& expr,
    const std::vector<double>& residual,
    const std::vector<int>& rows,
    std::vector<int>& feature_pool,
    int n_known_constants,
    int& n_total_constants,
    int mtry,
    std::uint32_t& split_seed,
    const std::vector<std::vector<int> >& feature_orders,
    std::vector<unsigned char>& row_in_node,
    std::vector<int>& ordered_rows) {
  ScenicSplit best{-1, 0.0, 0.0, 0.0};
  if (rows.size() < 2 || feature_pool.empty()) return best;

  double parent_sum = 0.0;
  double parent_sum_sq = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    row_in_node[rows[i]] = 1;
    const double y = residual[rows[i]];
    parent_sum += y;
    parent_sum_sq += y * y;
  }
  const double parent_sse = scenic_sse_from_sums(parent_sum, parent_sum_sq, rows.size());
  if (parent_sse <= 0.0) {
    for (std::size_t i = 0; i < rows.size(); ++i) row_in_node[rows[i]] = 0;
    return best;
  }

  ordered_rows.reserve(rows.size());
  const int n_features = static_cast<int>(feature_pool.size());
  int f_i = n_features;
  int n_visited = 0;
  int n_found_constants = 0;
  int n_drawn_constants = 0;
  n_known_constants = std::max(0, std::min(n_known_constants, n_features));
  n_total_constants = n_known_constants;
  while (f_i > n_total_constants && n_visited < mtry) {
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
    ordered_rows.clear();
    const std::vector<int>& order = feature_orders[feature];
    for (std::size_t oi = 0; oi < order.size(); ++oi) {
      if (row_in_node[order[oi]]) ordered_rows.push_back(order[oi]);
    }
    if (ordered_rows.size() < 2) continue;
    double min_x = expr(ordered_rows.front(), feature);
    double max_x = expr(ordered_rows.back(), feature);
    if (!R_finite(min_x)) min_x = 0.0;
    if (!R_finite(max_x)) max_x = 0.0;
    if (max_x <= min_x + SCENIC_FEATURE_THRESHOLD) {
      std::swap(feature_pool[f_j], feature_pool[n_total_constants]);
      ++n_found_constants;
      ++n_total_constants;
      continue;
    }
    --f_i;
    std::swap(feature_pool[f_i], feature_pool[f_j]);

    double left_sum = 0.0;
    double left_sum_sq = 0.0;
    std::size_t pos = 0;
    while (pos + 1 < ordered_rows.size()) {
      double x = expr(ordered_rows[pos], feature);
      if (!R_finite(x)) x = 0.0;
      do {
        const int row = ordered_rows[pos];
        const double y = residual[row];
        left_sum += y;
        left_sum_sq += y * y;
        ++pos;
        if (pos >= ordered_rows.size()) break;
        double next_group_x = expr(ordered_rows[pos], feature);
        if (!R_finite(next_group_x)) next_group_x = 0.0;
        if (next_group_x > x + SCENIC_FEATURE_THRESHOLD) break;
      } while (pos < ordered_rows.size());
      if (pos >= ordered_rows.size()) break;
      double next_x = expr(ordered_rows[pos], feature);
      if (!R_finite(next_x)) next_x = 0.0;

      const int n_left = static_cast<int>(pos);
      const int n_right = static_cast<int>(ordered_rows.size()) - n_left;
      if (n_left < 1 || n_right < 1) continue;

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
      if (proxy_gain > best.proxy_gain) {
        best.feature = feature;
        best.threshold = (x + next_x) / 2.0;
        if (best.threshold == static_cast<double>(next_x) || !R_finite(best.threshold)) {
          best.threshold = x;
        }
        best.proxy_gain = proxy_gain;
        best.importance_gain = importance_gain;
      }
    }
  }
  for (std::size_t i = 0; i < rows.size(); ++i) row_in_node[rows[i]] = 0;
  return best;
}

static int scenic_build_tree(
    const NumericMatrix& expr,
    const std::vector<double>& residual,
    const std::vector<int>& rows,
    std::vector<int>& feature_pool,
    int n_known_constants,
    int depth,
    int max_depth,
    int mtry,
    std::uint32_t& split_seed,
    std::vector<ScenicTreeNode>& nodes,
    std::vector<double>& importance,
    const std::vector<std::vector<int> >& feature_orders,
    std::vector<unsigned char>& row_in_node,
    std::vector<int>& ordered_rows) {
  ScenicTreeNode node;
  node.feature = -1;
  node.threshold = 0.0;
  node.value = scenic_node_mean(rows, residual);
  node.left = -1;
  node.right = -1;
  const int node_index = static_cast<int>(nodes.size());
  nodes.push_back(node);

  if (depth >= max_depth || rows.size() < 2 || feature_pool.empty()) return node_index;

  int n_total_constants = n_known_constants;
  const ScenicSplit split = scenic_best_split(
    expr, residual, rows, feature_pool, n_known_constants, n_total_constants,
    mtry, split_seed, feature_orders, row_in_node, ordered_rows);
  if (split.feature < 0 || split.proxy_gain <= 1e-12) return node_index;

  std::vector<int> left_rows;
  std::vector<int> right_rows;
  left_rows.reserve(rows.size());
  right_rows.reserve(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    double x = expr(rows[i], split.feature);
    if (!R_finite(x)) x = 0.0;
    if (x <= split.threshold) {
      left_rows.push_back(rows[i]);
    } else {
      right_rows.push_back(rows[i]);
    }
  }
  if (left_rows.empty() || right_rows.empty()) return node_index;

  importance[split.feature] += split.importance_gain;
  const int left = scenic_build_tree(
    expr, residual, left_rows, feature_pool, n_total_constants, depth + 1, max_depth, mtry, split_seed, nodes,
    importance, feature_orders, row_in_node, ordered_rows);
  const int right = scenic_build_tree(
    expr, residual, right_rows, feature_pool, n_total_constants, depth + 1, max_depth, mtry, split_seed, nodes,
    importance, feature_orders, row_in_node, ordered_rows);
  nodes[node_index].feature = split.feature;
  nodes[node_index].threshold = split.threshold;
  nodes[node_index].left = left;
  nodes[node_index].right = right;
  return node_index;
}

static double scenic_predict_tree(
    const NumericMatrix& expr,
    const std::vector<ScenicTreeNode>& nodes,
    int row) {
  int node_index = 0;
  while (node_index >= 0) {
    const ScenicTreeNode& node = nodes[node_index];
    if (node.feature < 0 || node.left < 0 || node.right < 0) return node.value;
    double x = expr(row, node.feature);
    if (!R_finite(x)) x = 0.0;
    node_index = (x <= node.threshold) ? node.left : node.right;
  }
  return 0.0;
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
      double x = expr(i, g);
      if (!R_finite(x)) x = 0.0;
      sum += x;
    }
    means[g] = sum / static_cast<double>(n_samples);
  }

  std::vector<int> all_rows(n_samples);
  std::iota(all_rows.begin(), all_rows.end(), 0);
  const std::vector<std::vector<int> > feature_orders = scenic_feature_orders(expr);
  const int sample_size = std::max(1, static_cast<int>(std::floor(subsample * n_samples)));
  std::vector<ScenicEdge> edges;
  std::vector<double> residual(n_samples);
  std::vector<double> importance(n_genes);
  std::vector<unsigned char> row_in_node(n_samples, 0);
  std::vector<int> ordered_rows;

  for (std::size_t ti = 0; ti < targets.size(); ++ti) {
    std::mt19937 rng(static_cast<unsigned int>(random_seed));
    const int target = targets[ti];
    std::vector<int> features;
    features.reserve(regs.size());
    for (std::size_t ri = 0; ri < regs.size(); ++ri) {
      if (!exclude_self || regs[ri] != target) features.push_back(regs[ri]);
    }
    if (features.empty()) continue;

    std::fill(importance.begin(), importance.end(), 0.0);
    for (int i = 0; i < n_samples; ++i) {
      double y = expr(i, target);
      if (!R_finite(y)) y = 0.0;
      residual[i] = y - means[target];
    }
    const int mtry = std::min<int>(
      features.size(),
      std::max(1, static_cast<int>(std::floor(max_features * features.size())))
    );
    std::vector<double> oob_improvements;
    oob_improvements.reserve(n_rounds);
    int actual_rounds = 0;
    double previous_oob_loss = NA_REAL;

    for (int round = 0; round < n_rounds; ++round) {
      std::vector<int> train_rows;
      train_rows.reserve(sample_size);
      std::vector<unsigned char> is_train = scenic_sample_mask(n_samples, sample_size, rng);
      for (int i = 0; i < n_samples; ++i) {
        if (is_train[i]) train_rows.push_back(i);
      }

      std::vector<ScenicTreeNode> nodes;
      nodes.reserve((1 << (max_depth + 1)) - 1);
      std::uint32_t split_seed = static_cast<std::uint32_t>(
        scenic_numpy_randint(rng, 0, SCENIC_RAND_R_MAX)
      );
      std::vector<int> feature_pool(features);
      scenic_build_tree(
        expr, residual, train_rows, feature_pool, 0, 0, max_depth, mtry, split_seed, nodes,
        importance, feature_orders, row_in_node, ordered_rows);
      if (nodes.empty()) break;

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

    std::vector<ScenicEdge> target_edges;
    target_edges.reserve(features.size());
    for (std::size_t ri = 0; ri < features.size(); ++ri) {
      const int reg = features[ri];
      if (importance[reg] > 0.0) {
        target_edges.push_back(ScenicEdge{
          reg + 1,
          target + 1,
          importance[reg] * static_cast<double>(actual_rounds) / total_importance
        });
      }
    }
    std::sort(target_edges.begin(), target_edges.end(), scenic_edge_before);
    const int take = max_edges_per_target > 0 ?
      std::min<int>(max_edges_per_target, target_edges.size()) :
      static_cast<int>(target_edges.size());
    for (int i = 0; i < take; ++i) edges.push_back(target_edges[i]);
    Rcpp::checkUserInterrupt();
  }

  std::sort(edges.begin(), edges.end(), scenic_edge_before);
  IntegerVector regulator(edges.size());
  IntegerVector target(edges.size());
  NumericVector importance_out(edges.size());
  for (std::size_t i = 0; i < edges.size(); ++i) {
    regulator[i] = edges[i].regulator;
    target[i] = edges[i].target;
    importance_out[i] = edges[i].importance;
  }

  return DataFrame::create(
    _["regulator"] = regulator,
    _["target"] = target,
    _["importance"] = importance_out,
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
