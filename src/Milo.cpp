#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

using namespace Rcpp;

static double median_in_place(std::vector<double>& values) {
  if (values.empty()) {
    return NA_REAL;
  }
  const std::size_t n = values.size();
  const std::size_t mid = n / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  const double upper = values[mid];
  if (n % 2 == 1) {
    return upper;
  }
  std::nth_element(values.begin(), values.begin() + mid - 1, values.begin() + mid);
  return 0.5 * (values[mid - 1] + upper);
}

// [[Rcpp::export]]
NumericMatrix milo_neighborhood_medians_cpp(NumericMatrix coords, IntegerMatrix knn_idx) {
  const int n_seeds = knn_idx.nrow();
  const int k = knn_idx.ncol();
  const int n_cells = coords.nrow();
  const int n_dims = coords.ncol();
  NumericMatrix medians(n_seeds, n_dims);
  std::vector<double> values;
  values.reserve(k);

  for (int seed = 0; seed < n_seeds; ++seed) {
    for (int dim = 0; dim < n_dims; ++dim) {
      values.clear();
      for (int col = 0; col < k; ++col) {
        const int idx = knn_idx(seed, col);
        if (idx == NA_INTEGER) {
          continue;
        }
        const int cell = idx - 1;
        if (cell < 0 || cell >= n_cells) {
          continue;
        }
        const double value = coords(cell, dim);
        if (R_finite(value)) {
          values.push_back(value);
        }
      }
      medians(seed, dim) = median_in_place(values);
    }
  }

  return medians;
}

// [[Rcpp::export]]
List milo_nhood_counts_cpp(
  IntegerMatrix knn_idx,
  IntegerVector sampled_vertices,
  IntegerVector sample_id,
  int n_samples,
  NumericVector k_dist
) {
  const int n_cells = knn_idx.nrow();
  const int k = knn_idx.ncol();
  const int n_nhoods = sampled_vertices.size();
  if (sample_id.size() != n_cells) {
    stop("sample_id length must match knn_idx rows");
  }
  if (k_dist.size() != n_cells) {
    stop("k_dist length must match knn_idx rows");
  }
  if (n_samples < 1) {
    stop("n_samples must be positive");
  }

  std::vector<std::vector<int>> adjacency(n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    adjacency[cell].reserve(k * 2 + 1);
    adjacency[cell].push_back(cell);
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int col = 0; col < k; ++col) {
      const int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) {
        continue;
      }
      const int nb0 = nb - 1;
      if (nb0 < 0 || nb0 >= n_cells || nb0 == cell) {
        continue;
      }
      adjacency[cell].push_back(nb0);
      adjacency[nb0].push_back(cell);
    }
  }

  IntegerMatrix counts(n_nhoods, n_samples);
  IntegerVector nhood_sizes(n_nhoods);
  NumericVector nhood_k_dist(n_nhoods);
  List members_out(n_nhoods);
  std::unordered_set<int> members;

  for (int nh = 0; nh < n_nhoods; ++nh) {
    const int seed = sampled_vertices[nh] - 1;
    if (seed < 0 || seed >= n_cells) {
      stop("sampled_vertices must be 1-based cell indices");
    }
    members.clear();
    for (int member : adjacency[seed]) {
      members.insert(member);
    }
    nhood_sizes[nh] = members.size();
    nhood_k_dist[nh] = k_dist[seed];
    for (int member : members) {
      const int sid = sample_id[member];
      if (sid == NA_INTEGER) {
        continue;
      }
      const int sid0 = sid - 1;
      if (sid0 < 0 || sid0 >= n_samples) {
        stop("sample_id must be 1-based sample indices within n_samples");
      }
      counts(nh, sid0) += 1;
    }
    std::vector<int> member_vec(members.begin(), members.end());
    std::sort(member_vec.begin(), member_vec.end());
    IntegerVector member_r(member_vec.size());
    for (std::size_t i = 0; i < member_vec.size(); ++i) {
      member_r[i] = member_vec[i] + 1;
    }
    members_out[nh] = member_r;
  }

  return List::create(
    _["counts"] = counts,
    _["nhood_size"] = nhood_sizes,
    _["k_distance"] = nhood_k_dist,
    _["sampled_vertices"] = sampled_vertices,
    _["members"] = members_out
  );
}

// [[Rcpp::export]]
NumericVector milo_weighted_fdr_cpp(NumericVector pvalues, NumericVector weights) {
  const int n = pvalues.size();
  if (weights.size() != n) {
    stop("weights length must match pvalues length");
  }

  std::vector<int> valid;
  valid.reserve(n);
  double weight_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(pvalues[i]) && R_finite(pvalues[i])) {
      double w = weights[i];
      if (!R_finite(w) || w <= 0.0) {
        w = 1.0;
      }
      weight_sum += w;
      valid.push_back(i);
    }
  }

  NumericVector out(n, NA_REAL);
  if (valid.empty()) {
    return out;
  }

  std::sort(valid.begin(), valid.end(), [&pvalues](int a, int b) {
    if (pvalues[a] < pvalues[b]) return true;
    if (pvalues[a] > pvalues[b]) return false;
    return a < b;
  });

  std::vector<double> cumulative(valid.size());
  double csum = 0.0;
  for (std::size_t rank = 0; rank < valid.size(); ++rank) {
    double w = weights[valid[rank]];
    if (!R_finite(w) || w <= 0.0) {
      w = 1.0;
    }
    csum += w;
    cumulative[rank] = weight_sum * pvalues[valid[rank]] / csum;
  }

  double running_min = R_PosInf;
  for (int rank = static_cast<int>(valid.size()) - 1; rank >= 0; --rank) {
    running_min = std::min(running_min, cumulative[rank]);
    out[valid[rank]] = std::min(1.0, running_min);
  }

  return out;
}
