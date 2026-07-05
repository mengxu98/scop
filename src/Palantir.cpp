// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <set>
#include <limits>
#include <cstdint>

using namespace Rcpp;

// ═══════════════════════════════════════════════════════════════════════════
// 1. ADAPTIVE ANISOTROPIC KERNEL
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List palantir_compute_kernel_cpp(
    NumericMatrix data,
    IntegerMatrix knn_idx,
    NumericMatrix knn_dist,
    int knn,
    double alpha = 0.0)
{
  int n = data.nrow();
  int k = knn_idx.ncol();

  // Per-cell adaptive bandwidth = distance to (knn/3)-th neighbor
  int adaptive_k = std::max(1, (int)std::floor(knn / 3.0) - 1);
  adaptive_k = std::min(adaptive_k, k - 1);
  NumericVector adaptive_std(n);
  for (int i = 0; i < n; ++i) {
    std::vector<double> row_dists(k);
    for (int j = 0; j < k; ++j) row_dists[j] = knn_dist(i, j);
    std::nth_element(row_dists.begin(), row_dists.begin() + adaptive_k, row_dists.end());
    adaptive_std[i] = row_dists[adaptive_k];
    if (adaptive_std[i] < 1e-10) adaptive_std[i] = 1e-10;
  }

  // Match palantir.utils.compute_kernel():
  // dists /= adaptive_std[x]; W = exp(-dists); kernel = W + W.T.
  std::vector<int> rows, cols;
  std::vector<double> vals;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      int nb = knn_idx(i, j) - 1;
      if (nb < 0 || nb >= n) continue;
      double dij = knn_dist(i, j);
      double w = std::exp(-dij / adaptive_std[i]);
      if (w > 1e-10) {
        rows.push_back(i + 1);
        cols.push_back(nb + 1);
        vals.push_back(w);
      }
    }
  }

  // Add W.T entries. Duplicate triplets are intentionally left for Matrix to
  // sum, matching scipy's sparse addition semantics.
  int ntriplets = vals.size();
  for (int t = 0; t < ntriplets; ++t) {
    int i = rows[t], j = cols[t];
    rows.push_back(j);
    cols.push_back(i);
    vals.push_back(vals[t]);
  }

  if (alpha > 0.0) {
    NumericVector degree(n, 0.0);
    for (int t = 0; t < (int)vals.size(); ++t) degree[rows[t] - 1] += vals[t];
    for (int i = 0; i < n; ++i)
      if (degree[i] > 0.0) degree[i] = std::pow(degree[i], -alpha);
    for (int t = 0; t < (int)vals.size(); ++t) {
      int i = rows[t] - 1, j = cols[t] - 1;
      vals[t] *= degree[i] * degree[j];
    }
  }

  return List::create(
    _["i"] = rows,
    _["j"] = cols,
    _["x"] = vals,
    _["n"] = n,
    _["adaptive_std"] = adaptive_std
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// 2. DIFFUSION MAP NORMALIZATION (Markov normalization)
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List palantir_normalize_kernel_cpp(
    IntegerVector kernel_i,
    IntegerVector kernel_j,
    NumericVector kernel_x,
    int n)
{
  // Row sums of kernel
  int nnz = kernel_i.size();
  NumericVector row_sum(n, 0.0);
  for (int t = 0; t < nnz; ++t) {
    int i = kernel_i[t] - 1;
    row_sum[i] += kernel_x[t];
  }
  for (int i = 0; i < n; ++i)
    if (row_sum[i] < 1e-10) row_sum[i] = 1.0;

  // Markov transition matrix: T(i,j) = K(i,j) / row_sum[i]
  std::vector<int> T_i, T_j;
  std::vector<double> T_x;
  for (int t = 0; t < nnz; ++t) {
    int i = kernel_i[t] - 1;
    double w = kernel_x[t] / row_sum[i];
    if (w > 1e-12) {
      T_i.push_back(i + 1);
      T_j.push_back(kernel_j[t]);
      T_x.push_back(w);
    }
  }

  // Eigenvalues of Laplacian
  // For the diffusion map, we need D^{-1/2} * K * D^{-1/2}
  NumericVector sqrt_inv_row(n);
  for (int i = 0; i < n; ++i)
    sqrt_inv_row[i] = 1.0 / std::sqrt(row_sum[i]);

  std::vector<int> L_i, L_j;
  std::vector<double> L_x;
  for (int t = 0; t < nnz; ++t) {
    int i = kernel_i[t] - 1, j = kernel_j[t] - 1;
    double w = kernel_x[t] * sqrt_inv_row[i] * sqrt_inv_row[j];
    if (w > 1e-12) {
      L_i.push_back(i + 1);
      L_j.push_back(j + 1);
      L_x.push_back(w);
    }
  }

  return List::create(
    _["T_i"] = T_i,
    _["T_j"] = T_j,
    _["T_x"] = T_x,
    _["L_i"] = L_i,
    _["L_j"] = L_j,
    _["L_x"] = L_x,
    _["n"] = n,
    _["row_sum"] = row_sum
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// 3. MULTISCALE SPACE (eigenvalue-scaled diffusion components)
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericMatrix palantir_multiscale_space_cpp(
    NumericMatrix eigenvectors,
    NumericVector eigenvalues)
{
  int n = eigenvectors.nrow();
  int n_eigs = eigenvectors.ncol();

  // Match palantir.utils.determine_multiscale_space(n_eigs = None).
  std::vector<std::pair<double, int> > gaps;
  for (int i = 0; i < n_eigs - 1; ++i) {
    gaps.push_back(std::make_pair(eigenvalues[i] - eigenvalues[i + 1], i));
  }
  std::sort(gaps.begin(), gaps.end(), std::greater<std::pair<double, int> >());
  int n_use = gaps.empty() ? n_eigs : gaps[0].second + 1;
  if (n_use < 3 && gaps.size() > 1) n_use = gaps[1].second + 1;
  if (n_use < 3) n_use = 3;
  n_use = std::min(n_use, n_eigs);

  // Scale: eigval / (1 - eigval) * eigvec
  NumericMatrix ms(n, n_use - 1);
  for (int i = 0; i < n; ++i)
    for (int j = 1; j < n_use; ++j)
      ms(i, j - 1) = eigenvectors(i, j) * eigenvalues[j] / (1.0 - eigenvalues[j]);

  return ms;
}

// ═══════════════════════════════════════════════════════════════════════════
// 4. MAX-MIN WAYPOINT SAMPLING
// ═══════════════════════════════════════════════════════════════════════════

static uint32_t numpy_random_bounded_uint32(std::mt19937& rng, uint32_t max_value) {
  uint32_t mask = max_value;
  mask |= mask >> 1;
  mask |= mask >> 2;
  mask |= mask >> 4;
  mask |= mask >> 8;
  mask |= mask >> 16;
  uint32_t value;
  do {
    value = rng() & mask;
  } while (value > max_value);
  return value;
}

// [[Rcpp::export]]
NumericVector palantir_numpy_random_sample_cpp(int n, int seed = 0) {
  std::mt19937 rng(seed);
  NumericVector out(n);
  const double scale = 1.0 / 9007199254740992.0; // 2^53
  for (int i = 0; i < n; ++i) {
    uint32_t a = rng() >> 5;
    uint32_t b = rng() >> 6;
    out[i] = (static_cast<double>(a) * 67108864.0 + static_cast<double>(b)) * scale;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector palantir_maxmin_waypoints_cpp(
    NumericMatrix ms_data,
    int num_waypoints,
    int seed = 20)
{
  int n = ms_data.nrow();
  int n_cols = ms_data.ncol();
  int min_waypoints = std::max(3, n_cols);
  if (num_waypoints < min_waypoints) num_waypoints = min_waypoints;

  std::mt19937 rng(seed);
  std::set<int> wp_set;

  int per_col = (int)std::floor((double)num_waypoints / n_cols);
  if (per_col < 1) per_col = 1;

  for (int col = 0; col < n_cols; ++col) {
    std::vector<double> vec(n);
    for (int i = 0; i < n; ++i) vec[i] = ms_data(i, col);

    int current = static_cast<int>(numpy_random_bounded_uint32(
      rng,
      static_cast<uint32_t>(n - 1)
    ));
    std::vector<int> iter_set;
    iter_set.push_back(current);
    std::vector<double> min_dists(n);
    for (int i = 0; i < n; ++i) min_dists[i] = std::abs(vec[i] - vec[current]);

    for (int k = 1; k < per_col; ++k) {
      double max_min_dist = -1.0;
      int best = current;
      for (int i = 0; i < n; ++i) {
        double min_dist = min_dists[i];
        if (min_dist > max_min_dist) {
          max_min_dist = min_dist;
          best = i;
        }
      }
      if (best != current) {
        iter_set.push_back(best);
        current = best;
        for (int i = 0; i < n; ++i) {
          double d = std::abs(vec[i] - vec[current]);
          if (d < min_dists[i]) min_dists[i] = d;
        }
      }
    }
    for (int w : iter_set) wp_set.insert(w);
  }

  IntegerVector result(wp_set.begin(), wp_set.end());
  for (int i = 0; i < result.size(); ++i) result[i] += 1;
  return result;
}

// ═══════════════════════════════════════════════════════════════════════════
// 5. PSEUDOTIME VIA SHORTEST PATHS + ITERATIVE REFINEMENT
// ═══════════════════════════════════════════════════════════════════════════

struct DijkstraCmp {
  bool operator()(const std::pair<double,int>& a, const std::pair<double,int>& b) const {
    return a.first > b.first;
  }
};


static void dijkstra_graph(
    int n,
    const std::vector<std::vector<std::pair<int,double>>>& graph,
    int source,
    std::vector<double>& dist)
{
  dist.assign(n, std::numeric_limits<double>::max());
  dist[source] = 0.0;
  std::priority_queue<std::pair<double,int>, std::vector<std::pair<double,int>>, DijkstraCmp> pq;
  pq.push({0.0, source});

  while (!pq.empty()) {
    auto top = pq.top(); pq.pop();
    double d = top.first;
    int u = top.second;
    if (d > dist[u]) continue;
    for (auto& edge : graph[u]) {
      int v = edge.first;
      double nd = d + edge.second;
      if (nd < dist[v]) {
        dist[v] = nd;
        pq.push({nd, v});
      }
    }
  }
}

// [[Rcpp::export]]
List palantir_pseudotime_cpp(
    NumericMatrix ms_data,
    int start_cell,
    IntegerVector waypoints,
    int knn,
    int max_iterations = 25,
    int n_jobs = 1)
{
  int n = ms_data.nrow();
  int d = ms_data.ncol();
  if (start_cell < 0 || start_cell >= n)
    stop("start_cell out of range");

  int n_wp = waypoints.size();

  // Build kNN graph using brute force (small data) or Euclidean distances
  // Since waypoints << n_cells usually, we build graph on ALL cells
  // Compute pairwise distances for kNN
  IntegerMatrix knn_idx_cpp(n, knn);
  NumericMatrix knn_dist_cpp(n, knn);

  for (int i = 0; i < n; ++i) {
    std::vector<std::pair<double,int>> dists;
    dists.reserve(n);
    for (int j = 0; j < n; ++j) {
      double dd = 0.0;
      for (int dim = 0; dim < d; ++dim) {
        double diff = ms_data(i, dim) - ms_data(j, dim);
        dd += diff * diff;
      }
      dists.push_back({dd, j});
    }
    int m = std::min(knn, (int)dists.size());
    std::nth_element(dists.begin(), dists.begin() + m - 1, dists.end());
    std::sort(dists.begin(), dists.begin() + m);
    for (int k = 0; k < m; ++k) {
      knn_idx_cpp(i, k) = dists[k].second + 1;
      knn_dist_cpp(i, k) = std::sqrt(dists[k].first);
    }
    for (int k = m; k < knn; ++k) {
      knn_idx_cpp(i, k) = NA_INTEGER;
      knn_dist_cpp(i, k) = NA_REAL;
    }
  }

  // Build adjacency from kNN
  std::vector<int> adj_i, adj_j;
  std::vector<double> adj_x;
  for (int i = 0; i < n; ++i) {
    for (int k_idx = 0; k_idx < knn; ++k_idx) {
      int j = knn_idx_cpp(i, k_idx) - 1;
      if (j == NA_INTEGER || j < 0) continue;
      adj_i.push_back(i);
      adj_j.push_back(j);
      adj_x.push_back(knn_dist_cpp(i, k_idx));
    }
  }

  // Dijkstra from all waypoints
  int n_sources = n_wp;
  std::vector<int> wp_offsets(n_wp);
  for (int w = 0; w < n_wp; ++w) wp_offsets[w] = waypoints[w] - 1;
  std::vector<std::vector<std::pair<int,double>>> graph(n);
  for (std::size_t t = 0; t < adj_i.size(); ++t) {
    int u = adj_i[t], v = adj_j[t];
    double w = adj_x[t];
    graph[u].push_back({v, w});
    graph[v].push_back({u, w});
  }

  NumericMatrix D(n_wp, n);
  NumericMatrix W(n_wp, n);

  for (int s = 0; s < n_wp; ++s) {
    std::vector<double> dists_v;
    dijkstra_graph(n, graph, wp_offsets[s], dists_v);
    for (int i = 0; i < n; ++i)
      D(s, i) = dists_v[i];
  }

  // Bandwidth for weight matrix
  double d_mean = 0.0;
  const double d_count = static_cast<double>(n_wp) * static_cast<double>(n);
  for (int s = 0; s < n_wp; ++s)
    for (int i = 0; i < n; ++i)
      d_mean += D(s, i);
  d_mean /= d_count;
  double d_var = 0.0;
  for (int s = 0; s < n_wp; ++s) {
    for (int i = 0; i < n; ++i) {
      double diff = D(s, i) - d_mean;
      d_var += diff * diff;
    }
  }
  d_var /= d_count;
  double sdv = std::sqrt(d_var) * 1.06 * std::pow(d_count, -0.2);

  // Weight matrix: W(s,i) = exp(-D(s,i)^2 / (2*sdv^2)).
  // In pandas, `W = W / W.sum()` divides each column by its column sum.
  NumericVector column_weight(n, 0.0);
  for (int s = 0; s < n_wp; ++s) {
    for (int i = 0; i < n; ++i) {
      double w = std::exp(-0.5 * D(s, i) * D(s, i) / (sdv * sdv));
      W(s, i) = w;
      column_weight[i] += w;
    }
  }
  for (int i = 0; i < n; ++i) {
    double denom = column_weight[i];
    if (denom > 0) {
      for (int s = 0; s < n_wp; ++s) W(s, i) /= denom;
    }
  }

  // Start row: waypoint closest to start_cell
  int start_row = -1;
  for (int w = 0; w < n_wp; ++w)
    if (wp_offsets[w] == start_cell) { start_row = w; break; }
  if (start_row < 0) {
    // Find nearest waypoint
    double min_d = 1e100;
    for (int w = 0; w < n_wp; ++w) {
      if (D(w, start_cell) < min_d) { min_d = D(w, start_cell); start_row = w; }
    }
  }

  NumericVector pseudotime(n);
  for (int i = 0; i < n; ++i) pseudotime[i] = D(start_row, i);

  // Iterative refinement
  bool converged = false;
  int iteration = 1;
  while (!converged && iteration < max_iterations) {
    NumericVector new_pt(n, 0.0);
    for (int i = 0; i < n; ++i) {
      double t_i = 0.0;
      for (int s = 0; s < n_wp; ++s) {
        double t_wp = pseudotime[wp_offsets[s]];
        double sign = (pseudotime[i] < t_wp) ? -1.0 : 1.0;
        if (s == start_row) sign = 1.0;
        double P_si = D(s, i) * sign + t_wp;
        t_i += P_si * W(s, i);
      }
      new_pt[i] = t_i;
    }

    // Check convergence
    double corr_num = 0.0, corr_den_x = 0.0, corr_den_y = 0.0;
    double mx = 0.0, my = 0.0;
    for (int i = 0; i < n; ++i) { mx += pseudotime[i]; my += new_pt[i]; }
    mx /= n; my /= n;
    for (int i = 0; i < n; ++i) {
      double dx = pseudotime[i] - mx;
      double dy = new_pt[i] - my;
      corr_num += dx * dy;
      corr_den_x += dx * dx;
      corr_den_y += dy * dy;
    }
    double corr = corr_num / (std::sqrt(corr_den_x) * std::sqrt(corr_den_y) + 1e-12);
    if (corr > 0.9999) converged = true;

    pseudotime = new_pt;
    ++iteration;
  }

  // Normalize to [0, 1]
  double pt_min = pseudotime[0], pt_max = pseudotime[0];
  for (int i = 0; i < n; ++i) {
    if (pseudotime[i] < pt_min) pt_min = pseudotime[i];
    if (pseudotime[i] > pt_max) pt_max = pseudotime[i];
  }
  double pt_range = pt_max - pt_min;
  if (pt_range > 0)
    for (int i = 0; i < n; ++i) pseudotime[i] = (pseudotime[i] - pt_min) / pt_range;

  return List::create(
    _["pseudotime"] = pseudotime,
    _["W"] = W
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// 6. MARKOV CHAIN & ABSORPTION PROBABILITIES
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List palantir_markov_chain_cpp(
    NumericMatrix wp_data,
    int knn,
    NumericVector pseudotime)
{
  int n = wp_data.nrow();
  int d = wp_data.ncol();

  // kNN on waypoints
  IntegerMatrix knn_idx(n, knn);
  NumericMatrix knn_dist(n, knn);

  for (int i = 0; i < n; ++i) {
    std::vector<std::pair<double,int>> dists;
    for (int j = 0; j < n; ++j) {
      double dd = 0.0;
      for (int dim = 0; dim < d; ++dim) {
        double diff = wp_data(i, dim) - wp_data(j, dim);
        dd += diff * diff;
      }
      dists.push_back({dd, j});
    }
    int m = std::min(knn, (int)dists.size());
    std::nth_element(dists.begin(), dists.begin() + m - 1, dists.end());
    std::sort(dists.begin(), dists.begin() + m);
    for (int k = 0; k < m; ++k) {
      knn_idx(i, k) = dists[k].second + 1;
      knn_dist(i, k) = std::sqrt(dists[k].first);
    }
    for (int k = m; k < knn; ++k) {
      knn_idx(i, k) = NA_INTEGER;
      knn_dist(i, k) = NA_REAL;
    }
  }

  // Adaptive bandwidth = distance to (knn/3)-th neighbor
  int adaptive_k = std::max(1, std::min((int)std::floor(knn / 3.0) - 1, knn - 1));
  NumericVector adaptive_std(n);
  for (int i = 0; i < n; ++i) {
    std::vector<double> row_dists(knn);
    for (int j = 0; j < knn; ++j) {
      row_dists[j] = knn_dist(i, j);
      if (!std::isfinite(row_dists[j])) row_dists[j] = 1e10;
    }
    std::nth_element(row_dists.begin(), row_dists.begin() + adaptive_k, row_dists.end());
    adaptive_std[i] = row_dists[adaptive_k];
    if (adaptive_std[i] < 1e-10) adaptive_std[i] = 1e-10;
  }

  // Directed graph: keep edges where neighbor is not too far back in pseudotime
  std::vector<int> T_i, T_j;
  std::vector<double> T_x;

  for (int i = 0; i < n; ++i) {
    double pt_i = pseudotime[i];
    double cutoff = pt_i - adaptive_std[i];
    double row_sum = 0.0;

    std::vector<double> local_w;
    std::vector<int> local_j;

    for (int k_idx = 0; k_idx < knn; ++k_idx) {
      if (knn_idx(i, k_idx) == NA_INTEGER) continue;
      int j = knn_idx(i, k_idx) - 1;
      if (j < 0) continue;
      double dist_ij = knn_dist(i, k_idx);

      // Only keep edges forward in pseudotime (not too far back)
      double pt_j = pseudotime[j];
      if (pt_j < cutoff) continue;

      double w = std::exp(-0.5 * dist_ij * dist_ij /
        (adaptive_std[i] * adaptive_std[j]));
      local_w.push_back(w);
      local_j.push_back(j);
      row_sum += w;
    }

    if (row_sum > 0) {
      for (size_t t = 0; t < local_w.size(); ++t) {
        T_i.push_back(i + 1);
        T_j.push_back(local_j[t] + 1);
        T_x.push_back(local_w[t] / row_sum);
      }
    } else {
      // Self-loop if no valid neighbors
      T_i.push_back(i + 1);
      T_j.push_back(i + 1);
      T_x.push_back(1.0);
    }
  }

  return List::create(
    _["T_i"] = T_i,
    _["T_j"] = T_j,
    _["T_x"] = T_x,
    _["n"] = n
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// 7. TERMINAL STATE DETECTION via Markov chain eigen analysis
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List palantir_terminal_states_cpp(
    IntegerVector T_i,
    IntegerVector T_j,
    NumericVector T_x,
    int n,
    NumericMatrix wp_data)
{
  // This is done in R using RSpectra or eigen for eigenvalue decomposition
  // This C++ function computes the ranks and cutoff threshold
  // Called from R after eigen decomposition

  // R will:
  // 1. Compute leading eigenvector of T (using RSpectra)
  // 2. Pass eigenvalues/vectors here
  // 3. We compute cutoff and terminal state assignment

  // This is a stub — the heavy lifting (eigendecomposition) is done in R
  // Placeholder for the cutoff computation

  return List::create(
    _["n"] = n
  );
}

// ═══════════════════════════════════════════════════════════════════════════
// 8. ABSORPTION PROBABILITIES via (I-Q) solve
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericMatrix palantir_absorption_cpp(
    IntegerVector T_i,
    IntegerVector T_j,
    NumericVector T_x,
    int n,
    IntegerVector terminal_state_indices)
{
  // Terminal state indices (0-based from caller, but we receive 1-based)
  int n_abs = terminal_state_indices.size();
  std::set<int> abs_set;
  for (int t = 0; t < n_abs; ++t) abs_set.insert(terminal_state_indices[t] - 1);

  // Identify transient states
  std::vector<int> trans;
  for (int i = 0; i < n; ++i)
    if (abs_set.find(i) == abs_set.end()) trans.push_back(i);

  int n_trans = trans.size();
  if (n_trans == 0 || n_abs == 0) {
    NumericMatrix bp(n, n_abs);
    for (int i = 0; i < n_abs; ++i)
      bp(terminal_state_indices[i] - 1, i) = 1.0;
    return bp;
  }

  // Map old indices to new indices in I-Q
  std::vector<int> new_idx(n, -1);
  for (int i = 0; i < n_trans; ++i) new_idx[trans[i]] = i;

  // Build I-Q matrix (dense — acceptable since n = num waypoints, typically < 2000)
  arma::mat I_Q(n_trans, n_trans, arma::fill::zeros);
  arma::mat R(n_trans, n_abs, arma::fill::zeros);

  for (int t = 0; t < (int)T_i.size(); ++t) {
    int i = T_i[t] - 1, j = T_j[t] - 1;
    double val = T_x[t];
    if (new_idx[i] >= 0 && new_idx[j] >= 0) {
      // Both transient: add to I-Q
      I_Q(new_idx[i], new_idx[j]) -= val;
    } else if (new_idx[i] >= 0 && abs_set.find(j) != abs_set.end()) {
      // i = transient, j = absorbing: add to R
      int abs_col = 0;
      for (int a = 0; a < n_abs; ++a)
        if (terminal_state_indices[a] - 1 == j) { abs_col = a; break; }
      R(new_idx[i], abs_col) += val;
    }
  }

  // Set diagonal of I-Q to 1
  for (int i = 0; i < n_trans; ++i) I_Q(i, i) += 1.0;

  // Solve (I-Q) * B = R using Armadillo
  arma::mat B;
  bool solve_ok = arma::solve(B, I_Q, R);
  if (!solve_ok) {
    // Fallback: pseudoinverse
    B = arma::pinv(I_Q) * R;
  }

  // Assemble full absorption matrix
  NumericMatrix bp(n, n_abs);
  for (int i = 0; i < n_abs; ++i)
    bp(terminal_state_indices[i] - 1, i) = 1.0;

  for (int i = 0; i < n_trans; ++i) {
    for (int a = 0; a < n_abs; ++a) {
      double p = B(i, a);
      if (p < 0) p = 0.0;
      bp(trans[i], a) = p;
    }
  }

  // Row-normalize
  for (int i = 0; i < n; ++i) {
    double rs = 0.0;
    for (int a = 0; a < n_abs; ++a) rs += bp(i, a);
    if (rs > 0)
      for (int a = 0; a < n_abs; ++a) bp(i, a) /= rs;
  }

  return bp;
}
