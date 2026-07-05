// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <random>

using namespace Rcpp;
using namespace arma;

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
