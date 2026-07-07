// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <limits>
#include <random>
#include <unordered_map>

using namespace Rcpp;
using namespace arma;

static std::uint64_t phate_edge_key(int row, int col) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(row)) << 32) |
    static_cast<std::uint32_t>(col);
}

// [[Rcpp::export]]
List phate_graphtools_affinity_data_cpp(
    NumericMatrix data,
    int knn,
    double decay,
    double thresh = 1e-4,
    int knn_max = -1)
{
  const int n_cells = data.nrow();
  const int n_features = data.ncol();
  if (n_cells < 1) stop("data must contain at least one row");
  if (n_features < 1) stop("data must contain at least one column");
  knn = std::max(1, std::min(knn, std::max(1, n_cells - 1)));
  const int max_neighbors = (knn_max > 0) ?
    std::min(knn_max, std::max(1, n_cells - 1)) :
    std::max(1, n_cells - 1);

  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> vals;
  rows.reserve(static_cast<std::size_t>(n_cells) * static_cast<std::size_t>(max_neighbors + 1));
  cols.reserve(rows.capacity());
  vals.reserve(rows.capacity());

  std::vector<std::pair<double, int> > distances;
  distances.reserve(std::max(0, n_cells - 1));
  for (int i = 0; i < n_cells; ++i) {
    distances.clear();
    for (int j = 0; j < n_cells; ++j) {
      if (i == j) continue;
      double dist_sq = 0.0;
      for (int f = 0; f < n_features; ++f) {
        const double diff = data(i, f) - data(j, f);
        dist_sq += diff * diff;
      }
      distances.push_back(std::make_pair(std::sqrt(std::max(0.0, dist_sq)), j));
    }
    const int take = std::min(max_neighbors, static_cast<int>(distances.size()));
    if (take <= 0) continue;
    std::partial_sort(
      distances.begin(),
      distances.begin() + take,
      distances.end(),
      [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        if (a.first == b.first) return a.second < b.second;
        return a.first < b.first;
      }
    );
    const int bw_pos = std::min(knn, take) - 1;
    double bandwidth = distances[static_cast<std::size_t>(bw_pos)].first;
    if (!std::isfinite(bandwidth) || bandwidth <= 0.0) {
      bandwidth = std::numeric_limits<double>::epsilon();
    }
    for (int p = 0; p < take; ++p) {
      const double d = distances[static_cast<std::size_t>(p)].first;
      const int j = distances[static_cast<std::size_t>(p)].second;
      const double weight = std::exp(-std::pow(d / bandwidth, decay));
      if (std::isfinite(weight) && weight >= thresh) {
        rows.push_back(i);
        cols.push_back(j);
        vals.push_back(weight);
      }
    }
  }
  for (int i = 0; i < n_cells; ++i) {
    rows.push_back(i);
    cols.push_back(i);
    vals.push_back(1.0);
  }

  std::unordered_map<std::uint64_t, double> directed;
  directed.reserve(vals.size() * 2);
  for (std::size_t e = 0; e < vals.size(); ++e) {
    const std::uint64_t key = phate_edge_key(rows[e], cols[e]);
    directed[key] += vals[e];
  }
  std::unordered_map<std::uint64_t, double> sym;
  sym.reserve(directed.size() * 2);
  for (std::unordered_map<std::uint64_t, double>::const_iterator it = directed.begin();
       it != directed.end();
       ++it) {
    const int row = static_cast<int>(it->first >> 32);
    const int col = static_cast<int>(it->first & 0xffffffffU);
    const std::uint64_t reverse_key = phate_edge_key(col, row);
    double reverse_value = 0.0;
    std::unordered_map<std::uint64_t, double>::const_iterator rit = directed.find(reverse_key);
    if (rit != directed.end()) reverse_value = rit->second;
    sym[it->first] = 0.5 * (it->second + reverse_value);
    if (sym.find(reverse_key) == sym.end()) {
      sym[reverse_key] = 0.5 * (reverse_value + it->second);
    }
  }

  std::vector<int> sym_rows;
  std::vector<int> sym_cols;
  std::vector<double> sym_vals;
  sym_rows.reserve(sym.size());
  sym_cols.reserve(sym.size());
  sym_vals.reserve(sym.size());
  for (std::unordered_map<std::uint64_t, double>::const_iterator it = sym.begin();
       it != sym.end();
       ++it) {
    const double v = it->second;
    if (!std::isfinite(v) || v <= 0.0) continue;
    sym_rows.push_back(static_cast<int>(it->first >> 32));
    sym_cols.push_back(static_cast<int>(it->first & 0xffffffffU));
    sym_vals.push_back(v);
  }

  return List::create(
    _["affinity_rows"] = IntegerVector(sym_rows.begin(), sym_rows.end()),
    _["affinity_cols"] = IntegerVector(sym_cols.begin(), sym_cols.end()),
    _["affinity_vals"] = NumericVector(sym_vals.begin(), sym_vals.end()),
    _["n_cells"] = n_cells
  );
}

// ── 1. α-Decaying kernel affinity ─────────────────────────────────────────
// Builds a sparse affinity matrix from k-NN distances using the α-decaying
// kernel from PHATE (Moon, van Dijk et al., Nature Biotechnology 2019).
//
// For each cell i and neighbor j, affinity K(i,j) = exp(-(d_{ij} / σ_i)^α)
// where σ_i is the local bandwidth (distance to k-th neighbor, or median).
// The result is symmetrized: A = (K + K^T) / 2, then row-normalized.
//
// [[Rcpp::export]]
List phate_affinity_cpp(
    NumericMatrix knn_dist,
    IntegerMatrix knn_idx,
    double alpha_decay = 1.0,
    int bandwidth_k = -1)
{
  const int n_cells = knn_dist.nrow();
  const int k = knn_dist.ncol();
  const int bw = (bandwidth_k > 0 && bandwidth_k <= k) ? bandwidth_k : k;

  std::vector<int> rows, cols;
  std::vector<double> vals;
  rows.reserve(n_cells * k * 2);
  cols.reserve(n_cells * k * 2);
  vals.reserve(n_cells * k * 2);

  for (int i = 0; i < n_cells; ++i) {
    // Local bandwidth σ_i = distance to bw-th neighbor
    double sigma = 1.0;
    int valid_neighbors = 0;
    for (int j = 0; j < bw; ++j) {
      if (knn_idx(i, j) != NA_INTEGER && std::isfinite(knn_dist(i, j))) {
        sigma = knn_dist(i, j);
        valid_neighbors++;
      }
    }
    if (sigma <= 0.0 || !std::isfinite(sigma)) sigma = 1.0;
    // If alpha=1, this is just a Gaussian kernel scaled by 1/sigma_i
    // If alpha>1, steeper decay with distance

    for (int j = 0; j < k; ++j) {
      int nb = knn_idx(i, j);
      if (nb == NA_INTEGER) continue;
      nb -= 1;  // 0-based
      if (nb < 0 || nb >= n_cells || nb == i) continue;
      double d = knn_dist(i, j);
      if (!std::isfinite(d) || d <= 0.0) continue;

      double d_scaled = d / sigma;
      double weight = std::exp(-std::pow(d_scaled, alpha_decay));
      if (weight < 1e-4) continue;

      // Upper triangle for symmetric matrix
      rows.push_back(i);
      cols.push_back(nb);
      vals.push_back(weight);
    }
  }

  // Build sparse triplet list (symmetric, so we add both directions)
  int nnz = static_cast<int>(vals.size());
  std::vector<int> sym_rows, sym_cols;
  std::vector<double> sym_vals;
  sym_rows.reserve(nnz * 2);
  sym_cols.reserve(nnz * 2);
  sym_vals.reserve(nnz * 2);
  for (int e = 0; e < nnz; ++e) {
    sym_rows.push_back(rows[e]);
    sym_cols.push_back(cols[e]);
    sym_vals.push_back(vals[e]);
    // Symmetrize
    if (rows[e] != cols[e]) {
      sym_rows.push_back(cols[e]);
      sym_cols.push_back(rows[e]);
      sym_vals.push_back(vals[e]);
    }
  }

  return List::create(
    _["affinity_rows"] = IntegerVector(sym_rows.begin(), sym_rows.end()),
    _["affinity_cols"] = IntegerVector(sym_cols.begin(), sym_cols.end()),
    _["affinity_vals"] = NumericVector(sym_vals.begin(), sym_vals.end()),
    _["n_cells"] = n_cells
  );
}

// ── 2. Diffusion operator ──────────────────────────────────────────────────
// From sparse affinity triplets, computes the Markov transition matrix P
// via row-normalization, then powers it to diffusion time t.
//
// Returns P^t × (random walk starting distribution) for MDS input,
// or the powered transition matrix in log-space for potential distances.
//
// [[Rcpp::export]]
NumericMatrix phate_diffusion_operator_cpp(
    IntegerVector rows,
    IntegerVector cols,
    NumericVector vals,
    int n_cells,
    int t_max = 10)
{
  // Build row-normalized transition matrix P (dense for now, can be optimized)
  mat P = zeros<mat>(n_cells, n_cells);
  vec row_sums = zeros<vec>(n_cells);

  int nnz = rows.size();
  for (int e = 0; e < nnz; ++e) {
    int i = rows[e];
    int j = cols[e];
    double v = vals[e];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      P(i, j) += v;
      row_sums(i) += v;
    }
  }

  // Row normalize
  for (int i = 0; i < n_cells; ++i) {
    if (row_sums(i) > 1e-15) {
      for (int j = 0; j < n_cells; ++j) {
        P(i, j) /= row_sums(i);
      }
    } else {
      P(i, i) = 1.0;
    }
  }

  // Power iteration: P^t for t = 2..t_max
  mat Pt = P;
  for (int t = 1; t < t_max; ++t) {
    Pt = Pt * P;
  }

  // Convert to log-space: log(P^t + epsilon) for potential distance
  const double eps = 1e-10;
  mat logPt = zeros<mat>(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      logPt(i, j) = std::log(std::max(Pt(i, j), eps));
    }
  }

  // Return log-transformed diffusion probabilities
  NumericMatrix result(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      result(i, j) = logPt(i, j);
    }
  }
  return result;
}

// ── 3. Potential distance ─────────────────────────────────────────────────
// Computes D(i,j) = ||log P^t_{i,·} − log P^t_{j,·}||
// (Euclidean distance between rows of log-transformed diffusion matrix)
//
// Uses landmark-based approximation if n_landmarks < n_cells.
//
// [[Rcpp::export]]
NumericMatrix phate_potential_distance_cpp(
    NumericMatrix log_transition,
    int n_landmarks = -1)
{
  const int n = log_transition.nrow();
  if (log_transition.ncol() != n)
    stop("log_transition must be a square matrix");
  const int use_landmarks = (n_landmarks > 0 && n_landmarks < n) ? n_landmarks : n;

  mat logP(log_transition.begin(), n, n, false);

  // If landmark MDS: select landmarks randomly and compute distances to landmarks
  if (use_landmarks < n) {
    std::vector<int> landmark_idx(n);
    for (int i = 0; i < n; ++i) landmark_idx[i] = i;
    std::mt19937 rng(42);
    std::shuffle(landmark_idx.begin(), landmark_idx.end(), rng);
    landmark_idx.resize(use_landmarks);

    mat landmark_dist = zeros<mat>(n, use_landmarks);
    for (int i = 0; i < n; ++i) {
      for (int m = 0; m < use_landmarks; ++m) {
        int lm = landmark_idx[m];
        double dist = 0.0;
        for (int d = 0; d < n; ++d) {
          double diff = logP(i, d) - logP(lm, d);
          dist += diff * diff;
        }
        landmark_dist(i, m) = std::sqrt(std::max(dist, 0.0));
      }
    }

    NumericMatrix D(n, n);
    for (int i = 0; i < n; ++i) {
      D(i, i) = 0.0;
      for (int j = i + 1; j < n; ++j) {
        double dist = 0.0;
        for (int m = 0; m < use_landmarks; ++m) {
          double diff = landmark_dist(i, m) - landmark_dist(j, m);
          dist += diff * diff;
        }
        D(i, j) = std::sqrt(std::max(dist, 0.0));
        D(j, i) = D(i, j);
      }
    }
    return D;
  }

  // Full pairwise distance matrix
  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      double dist = 0.0;
      for (int d = 0; d < n; ++d) {
        double diff = logP(i, d) - logP(j, d);
        dist += diff * diff;
      }
      D(i, j) = std::sqrt(std::max(dist, 0.0));
      D(j, i) = D(i, j);
    }
  }
  return D;
}

// ── 4. Classical MDS (Metric MDS) ─────────────────────────────────────────
// Given a pairwise distance matrix D, computes the top n_components
// eigenvectors of the double-centered squared distance matrix.
//
// Returns: embedding coordinates (n × n_components)
//
// [[Rcpp::export]]
NumericMatrix phate_metric_mds_cpp(
    NumericMatrix D,
    int n_components = 2)
{
  const int n = D.nrow();
  if (D.ncol() != n)
    stop("D must be a square distance matrix");
  if (n < 1)
    stop("D must contain at least one row");
  if (n_components < 1)
    stop("n_components must be at least 1");
  if (n_components > n) n_components = n;

  mat dist(D.begin(), n, n, false);

  // Double-center the squared distance matrix: B = -0.5 * J * D² * J
  mat D2 = dist % dist;  // element-wise square

  // Centering matrix J = I - (1/n) * 1*1^T
  vec one_vec = arma::ones<vec>(n);
  mat J = eye<mat>(n, n) - (one_vec * one_vec.t()) / static_cast<double>(n);

  mat B = -0.5 * J * D2 * J;

  // Eigendecomposition of B (symmetric)
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, B);

  // Sort eigenvalues descending (Armadillo returns ascending)
  // Take top n_components
  NumericMatrix embedding(n, n_components);
  for (int c = 0; c < n_components; ++c) {
    int src = n - 1 - c;  // descending order
    double ev = eigval(src);
    if (ev <= 0.0) {
      // Fill remaining columns with zeros
      for (int i = 0; i < n; ++i) embedding(i, c) = 0.0;
    } else {
      double scale = std::sqrt(ev);
      for (int i = 0; i < n; ++i) {
        embedding(i, c) = eigvec(i, src) * scale;
      }
    }
  }

  return embedding;
}

// ── 5. Optimal diffusion time via Von Neumann entropy ─────────────────────
// Finds the optimal t by minimizing the rate of change of
// Von Neumann entropy of the powered diffusion operator.
// Returns: optimal t (1 ≤ t ≤ t_max)
//
// [[Rcpp::export]]
int phate_find_optimal_t_cpp(
    IntegerVector rows,
    IntegerVector cols,
    NumericVector vals,
    int n_cells,
    int t_max = 20)
{
  if (t_max < 1) t_max = 1;

  // Build P
  mat P = zeros<mat>(n_cells, n_cells);
  vec row_sums = zeros<vec>(n_cells);
  int nnz = rows.size();
  for (int e = 0; e < nnz; ++e) {
    int i = rows[e], j = cols[e];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      P(i, j) += vals[e];
      row_sums(i) += vals[e];
    }
  }
  for (int i = 0; i < n_cells; ++i) {
    if (row_sums(i) > 1e-15) {
      for (int j = 0; j < n_cells; ++j) P(i, j) /= row_sums(i);
    } else {
      P(i, i) = 1.0;
    }
  }

  // Compute VNE at each t
  mat Pt = P;
  std::vector<double> vne_vals(t_max);
  for (int t = 0; t < t_max; ++t) {
    if (t > 0) Pt = Pt * P;

    // Von Neumann entropy: -sum(λ * log(λ))
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, (Pt + Pt.t()) / 2.0);  // symmetrize for stability
    double vne = 0.0;
    for (int i = 0; i < n_cells; ++i) {
      double lambda = std::abs(eigval(i));
      if (lambda > 1e-15) vne -= lambda * std::log(lambda);
    }
    vne_vals[t] = vne;
  }

  // Find knee point: max second derivative (minimum rate of VNE change)
  int optimal_t = 1;
  double min_curvature = std::numeric_limits<double>::max();
  for (int t = 1; t < t_max - 1; ++t) {
    double curvature = std::abs(vne_vals[t+1] - 2.0 * vne_vals[t] + vne_vals[t-1]);
    if (curvature < min_curvature) {
      min_curvature = curvature;
      optimal_t = t + 1;
    }
  }

  return std::max(1, optimal_t);
}
