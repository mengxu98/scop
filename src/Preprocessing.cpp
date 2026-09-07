#include <Rcpp.h>
#include <thisutils/log_message.h>
#include "thread_utils.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <limits>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ── 1. Filter genes (matching scv.pp.filter_genes exactly) ──────────────────

// [[Rcpp::export]]
IntegerVector scanpy_filter_genes_cpp(
    NumericMatrix spliced,     // genes × cells
    NumericMatrix unspliced,
    int min_counts = 3,
    int min_counts_u = 3,
    int n_threads = 0)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    thisutils::log_message("spliced and unspliced must have identical dimensions", "error");

  IntegerVector keep(n_genes, 1);
  const int threads = omp_thread_count(n_threads, n_genes);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  #endif
  for (int g = 0; g < n_genes; ++g) {
    double sum_s = 0.0, sum_u = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      sum_s += spliced(g, c);
      sum_u += unspliced(g, c);
    }
    // scvelo does TWO PASSES: first filter on spliced, then on unspliced
    // Each pass keeps if sum >= min_counts
    if (sum_s < min_counts || sum_u < min_counts_u)
      keep[g] = 0;
  }
  return keep;
}

// ── 2. Normalize per cell (matching scv.pp.normalize_per_cell, without log1p) ──

// [[Rcpp::export]]
List scanpy_normalize_cpp(
    NumericMatrix spliced,     // genes × cells, ALREADY FILTERED
    NumericMatrix unspliced,   // genes × cells, ALREADY FILTERED
    NumericVector initial_spliced_totals,  // per-cell totals BEFORE filtering (length = n_cells)
    NumericVector initial_unspliced_totals)
{
  int n_genes = spliced.nrow();
  int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    thisutils::log_message("spliced and unspliced must have identical dimensions", "error");
  if (initial_spliced_totals.size() != n_cells || initial_unspliced_totals.size() != n_cells)
    thisutils::log_message("initial totals must match n_cells", "error");

  // Median of pre-filtering totals (matching scvelo's get_initial_size)
  std::vector<double> sorted_s(n_cells), sorted_u(n_cells);
  for (int c = 0; c < n_cells; ++c) {
    sorted_s[c] = initial_spliced_totals[c];
    sorted_u[c] = initial_unspliced_totals[c];
  }
  std::sort(sorted_s.begin(), sorted_s.end());
  std::sort(sorted_u.begin(), sorted_u.end());
  // numpy.median averages the two middle values for even-length input;
  // scvelo's normalize_per_cell uses np.median(get_initial_size(...)).
  const bool even = n_cells % 2 == 0;
  const int mid = n_cells / 2;
  double median_s = even
    ? 0.5 * (sorted_s[mid - 1] + sorted_s[mid])
    : sorted_s[mid];
  double median_u = even
    ? 0.5 * (sorted_u[mid - 1] + sorted_u[mid])
    : sorted_u[mid];
  if (median_s <= 0) median_s = 1.0;
  if (median_u <= 0) median_u = 1.0;

  // Normalize in place and return the input matrices, avoiding two dense
  // 84 MB-class copies for large single-cell datasets.
  for (int c = 0; c < n_cells; ++c) {
    double scale_s = initial_spliced_totals[c] > 0
      ? median_s / initial_spliced_totals[c] : 1.0;
    double scale_u = initial_unspliced_totals[c] > 0
      ? median_u / initial_unspliced_totals[c] : 1.0;
    for (int g = 0; g < n_genes; ++g) {
      spliced(g, c) *= scale_s;
      unspliced(g, c) *= scale_u;
    }
  }

  return List::create(
    _["spliced_norm"] = spliced,
    _["unspliced_norm"] = unspliced
  );
}


// ── 3. KNN via exact brute-force Euclidean (deterministic, matching scanpy) ──

// [[Rcpp::export]]
List scanpy_knn_cpp(
    NumericMatrix coords,
    int n_neighbors = 10,
    bool exclude_self = true,
    int n_threads = 0)
{
  int n_cells = coords.nrow();
  int n_dims = coords.ncol();
  if (n_neighbors > n_cells) n_neighbors = n_cells;
  const int threads = omp_thread_count(n_threads, n_cells);

  IntegerMatrix idx(n_cells, n_neighbors);
  NumericMatrix dist(n_cells, n_neighbors);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int i = 0; i < n_cells; ++i) {
    std::vector<std::pair<double, int>> dists;
    dists.reserve(n_cells);

    for (int j = 0; j < n_cells; ++j) {
      if (exclude_self && j == i) continue;
      double d = 0.0;
      for (int dim = 0; dim < n_dims; ++dim) {
        double delta = coords(i, dim) - coords(j, dim);
        d += delta * delta;
      }
      dists.push_back({d, j});
    }

    auto dist_less = [](const std::pair<double, int>& a,
                        const std::pair<double, int>& b) {
      if (a.first < b.first) return true;
      if (a.first > b.first) return false;
      return a.second < b.second;
    };
    const int keep = std::min(n_neighbors, static_cast<int>(dists.size()));
    if (keep > 0) {
      std::partial_sort(
        dists.begin(),
        dists.begin() + keep,
        dists.end(),
        dist_less
      );
    }

    for (int k = 0; k < keep; ++k) {
      idx(i, k) = dists[k].second + 1;
      dist(i, k) = std::sqrt(dists[k].first);
    }
    for (int k = keep; k < n_neighbors; ++k) {
      idx(i, k) = NA_INTEGER;
      dist(i, k) = NA_REAL;
    }
  }

  return List::create(_["idx"] = idx, _["dist"] = dist);
}

