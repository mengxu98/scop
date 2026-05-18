#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;

struct PagaEdge {
  int i;
  int j;
  double weight;
};

struct PagaDsu {
  std::vector<int> parent;
  std::vector<int> rank;

  explicit PagaDsu(int n) : parent(n), rank(n, 0) {
    std::iota(parent.begin(), parent.end(), 0);
  }

  int find(int x) {
    if (parent[x] != x) {
      parent[x] = find(parent[x]);
    }
    return parent[x];
  }

  bool unite(int a, int b) {
    int root_a = find(a);
    int root_b = find(b);
    if (root_a == root_b) {
      return false;
    }
    if (rank[root_a] < rank[root_b]) {
      std::swap(root_a, root_b);
    }
    parent[root_b] = root_a;
    if (rank[root_a] == rank[root_b]) {
      ++rank[root_a];
    }
    return true;
  }
};

// [[Rcpp::export]]
List paga_connectivities_cpp(IntegerMatrix knn_idx, IntegerVector groups, int n_groups) {
  const int n_cells = groups.size();
  if (knn_idx.nrow() != n_cells) {
    stop("knn_idx rows must match groups length");
  }
  if (n_groups < 1) {
    stop("n_groups must be positive");
  }

  NumericMatrix directed_edges(n_groups, n_groups);
  NumericVector group_sizes(n_groups);

  for (int cell = 0; cell < n_cells; ++cell) {
    const int group = groups[cell] - 1;
    if (group < 0 || group >= n_groups) {
      stop("groups must be 1-based integers within n_groups");
    }
    group_sizes[group] += 1.0;
  }

  for (int cell = 0; cell < n_cells; ++cell) {
    const int from_group = groups[cell] - 1;
    for (int col = 0; col < knn_idx.ncol(); ++col) {
      const int neighbor = knn_idx(cell, col);
      if (neighbor == NA_INTEGER) {
        continue;
      }
      const int neighbor0 = neighbor - 1;
      if (neighbor0 < 0 || neighbor0 >= n_cells || neighbor0 == cell) {
        continue;
      }
      const int to_group = groups[neighbor0] - 1;
      if (to_group < 0 || to_group >= n_groups) {
        continue;
      }
      directed_edges(from_group, to_group) += 1.0;
    }
  }

  NumericVector edge_totals(n_groups);
  for (int i = 0; i < n_groups; ++i) {
    for (int j = 0; j < n_groups; ++j) {
      edge_totals[i] += directed_edges(i, j);
    }
  }

  NumericMatrix connectivities(n_groups, n_groups);
  NumericMatrix expected_edges(n_groups, n_groups);

  const double denom = static_cast<double>(std::max(n_cells - 1, 1));
  for (int i = 0; i < n_groups; ++i) {
    for (int j = i + 1; j < n_groups; ++j) {
      const double observed = directed_edges(i, j) + directed_edges(j, i);
      if (observed <= 0.0) {
        continue;
      }
      const double expected =
        (edge_totals[i] * group_sizes[j] + edge_totals[j] * group_sizes[i]) / denom;
      double scaled = expected > 0.0 ? observed / expected : 1.0;
      if (scaled > 1.0) {
        scaled = 1.0;
      }
      connectivities(i, j) = scaled;
      connectivities(j, i) = scaled;
      expected_edges(i, j) = expected;
      expected_edges(j, i) = expected;
    }
  }

  std::vector<PagaEdge> edges;
  edges.reserve(static_cast<std::size_t>(n_groups * (n_groups - 1) / 2));
  for (int i = 0; i < n_groups; ++i) {
    for (int j = i + 1; j < n_groups; ++j) {
      const double weight = connectivities(i, j);
      if (weight > 0.0 && std::isfinite(weight)) {
        PagaEdge edge;
        edge.i = i;
        edge.j = j;
        edge.weight = weight;
        edges.push_back(edge);
      }
    }
  }

  std::sort(edges.begin(), edges.end(), [](const PagaEdge& a, const PagaEdge& b) {
    if (a.weight > b.weight) {
      return true;
    }
    if (a.weight < b.weight) {
      return false;
    }
    if (a.i < b.i) {
      return true;
    }
    if (a.i > b.i) {
      return false;
    }
    return a.j < b.j;
  });

  NumericMatrix connectivities_tree(n_groups, n_groups);
  PagaDsu dsu(n_groups);
  int used_edges = 0;
  for (const PagaEdge& edge : edges) {
    if (dsu.unite(edge.i, edge.j)) {
      connectivities_tree(edge.i, edge.j) = edge.weight;
      connectivities_tree(edge.j, edge.i) = edge.weight;
      ++used_edges;
      if (used_edges == n_groups - 1) {
        break;
      }
    }
  }

  return List::create(
    _["connectivities"] = connectivities,
    _["connectivities_tree"] = connectivities_tree,
    _["expected_n_edges_random"] = expected_edges,
    _["group_sizes"] = group_sizes,
    _["directed_edges"] = directed_edges
  );
}
