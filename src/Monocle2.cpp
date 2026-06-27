#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

using namespace Rcpp;

namespace {

double squared_distance_columns(const NumericMatrix& x, int a, int b) {
  double out = 0.0;
  for (int r = 0; r < x.nrow(); ++r) {
    const double diff = x(r, a) - x(r, b);
    out += diff * diff;
  }
  return out;
}

double euclidean_distance_columns(const NumericMatrix& x, int a, int b) {
  return std::sqrt(squared_distance_columns(x, a, b));
}

double euclidean_distance_vector_column(
    const std::vector<double>& point,
    const NumericMatrix& x,
    int col) {
  double out = 0.0;
  for (int r = 0; r < x.nrow(); ++r) {
    const double diff = point[r] - x(r, col);
    out += diff * diff;
  }
  return std::sqrt(out);
}

} // namespace

// [[Rcpp::export]]
List monocle2_order_from_mst_cpp(
    NumericMatrix distances,
    IntegerMatrix edges,
    int root_cell = 1) {
  const int n_cells = distances.nrow();
  if (n_cells < 1) stop("distances must contain at least one cell");
  if (distances.ncol() != n_cells) stop("distances must be a square matrix");
  if (edges.ncol() != 2) stop("edges must have two columns");

  std::vector<std::vector<int>> adj(n_cells);
  for (int i = 0; i < edges.nrow(); ++i) {
    const int from = edges(i, 0) - 1;
    const int to = edges(i, 1) - 1;
    if (from < 0 || from >= n_cells || to < 0 || to >= n_cells || from == to) {
      continue;
    }
    adj[from].push_back(to);
    adj[to].push_back(from);
  }
  for (int i = 0; i < n_cells; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
  }

  int root = root_cell - 1;
  if (root < 0 || root >= n_cells) root = 0;

  std::vector<int> parent(n_cells, -1);
  std::vector<int> order;
  order.reserve(n_cells);
  std::vector<std::pair<int, int>> stack;
  stack.push_back(std::make_pair(root, 0));
  parent[root] = -2;
  while (!stack.empty()) {
    int node = stack.back().first;
    int& next_idx = stack.back().second;
    if (next_idx == 0) order.push_back(node);
    if (next_idx >= static_cast<int>(adj[node].size())) {
      stack.pop_back();
      continue;
    }
    const int next = adj[node][next_idx++];
    if (parent[next] != -1) continue;
    parent[next] = node;
    stack.push_back(std::make_pair(next, 0));
  }

  NumericVector pseudotime(n_cells);
  IntegerVector state(n_cells);
  CharacterVector parent_out(n_cells);
  int curr_state = 1;
  for (int i = 0; i < n_cells; ++i) {
    state[i] = 1;
    parent_out[i] = NA_STRING;
  }

  for (int idx = 0; idx < static_cast<int>(order.size()); ++idx) {
    const int node = order[idx];
    const int par = parent[node];
    if (par >= 0) {
      pseudotime[node] = pseudotime[par] + distances(node, par);
      if (static_cast<int>(adj[par].size()) > 2) {
        ++curr_state;
      }
      parent_out[node] = std::to_string(par + 1);
    } else {
      pseudotime[node] = 0.0;
    }
    state[node] = curr_state;
  }
  parent[root] = -1;

  IntegerVector order_out(order.size());
  for (int i = 0; i < static_cast<int>(order.size()); ++i) {
    order_out[i] = order[i] + 1;
  }

  return List::create(
    _["pseudotime"] = pseudotime,
    _["state"] = state,
    _["parent"] = parent_out,
    _["order"] = order_out,
    _["root_cell"] = root + 1,
    _["state_count"] = curr_state
  );
}

// [[Rcpp::export]]
List monocle2_order_from_weighted_edges_cpp(
    int n_cells,
    IntegerMatrix edges,
    NumericVector weights,
    int root_cell = 1) {
  if (n_cells < 1) stop("n_cells must be positive");
  if (edges.ncol() != 2) stop("edges must have two columns");
  if (weights.size() != edges.nrow()) stop("weights length must match edge rows");

  std::vector<std::vector<std::pair<int, double>>> adj(n_cells);
  for (int i = 0; i < edges.nrow(); ++i) {
    const int from = edges(i, 0) - 1;
    const int to = edges(i, 1) - 1;
    if (from < 0 || from >= n_cells || to < 0 || to >= n_cells || from == to) {
      continue;
    }
    const double weight = weights[i];
    adj[from].push_back(std::make_pair(to, weight));
    adj[to].push_back(std::make_pair(from, weight));
  }
  for (int i = 0; i < n_cells; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
  }

  int root = root_cell - 1;
  if (root < 0 || root >= n_cells) root = 0;

  std::vector<int> parent(n_cells, -1);
  std::vector<double> parent_weight(n_cells, 0.0);
  std::vector<int> order;
  order.reserve(n_cells);
  std::vector<std::pair<int, int>> stack;
  stack.push_back(std::make_pair(root, 0));
  parent[root] = -2;
  while (!stack.empty()) {
    int node = stack.back().first;
    int& next_idx = stack.back().second;
    if (next_idx == 0) order.push_back(node);
    if (next_idx >= static_cast<int>(adj[node].size())) {
      stack.pop_back();
      continue;
    }
    const int next = adj[node][next_idx].first;
    const double weight = adj[node][next_idx].second;
    ++next_idx;
    if (parent[next] != -1) continue;
    parent[next] = node;
    parent_weight[next] = weight;
    stack.push_back(std::make_pair(next, 0));
  }

  NumericVector pseudotime(n_cells);
  IntegerVector state(n_cells);
  CharacterVector parent_out(n_cells);
  int curr_state = 1;
  for (int i = 0; i < n_cells; ++i) {
    state[i] = 1;
    parent_out[i] = NA_STRING;
  }

  for (int idx = 0; idx < static_cast<int>(order.size()); ++idx) {
    const int node = order[idx];
    const int par = parent[node];
    if (par >= 0) {
      pseudotime[node] = pseudotime[par] + parent_weight[node];
      if (static_cast<int>(adj[par].size()) > 2) {
        ++curr_state;
      }
      parent_out[node] = std::to_string(par + 1);
    } else {
      pseudotime[node] = 0.0;
    }
    state[node] = curr_state;
  }

  IntegerVector order_out(order.size());
  for (int i = 0; i < static_cast<int>(order.size()); ++i) {
    order_out[i] = order[i] + 1;
  }

  return List::create(
    _["pseudotime"] = pseudotime,
    _["state"] = state,
    _["parent"] = parent_out,
    _["order"] = order_out,
    _["root_cell"] = root + 1,
    _["state_count"] = curr_state
  );
}

// [[Rcpp::export]]
List monocle2_project_cells_to_mst_cpp(
    NumericMatrix z,
    NumericMatrix y,
    IntegerMatrix graph_edges,
    IntegerVector closest_vertex) {
  const int dims = z.nrow();
  const int n_cells = z.ncol();
  const int n_centers = y.ncol();
  if (y.nrow() != dims) stop("z and y must have the same number of rows");
  if (closest_vertex.size() != n_cells) {
    stop("closest_vertex length must equal the number of cells");
  }
  if (graph_edges.ncol() != 2) stop("graph_edges must have two columns");

  std::vector<std::vector<int>> center_adj(n_centers);
  for (int i = 0; i < graph_edges.nrow(); ++i) {
    const int from = graph_edges(i, 0) - 1;
    const int to = graph_edges(i, 1) - 1;
    if (from < 0 || from >= n_centers || to < 0 || to >= n_centers || from == to) {
      continue;
    }
    center_adj[from].push_back(to);
    center_adj[to].push_back(from);
  }
  std::vector<int> is_tip(n_centers, 0);
  for (int i = 0; i < n_centers; ++i) {
    std::sort(center_adj[i].begin(), center_adj[i].end());
    is_tip[i] = center_adj[i].size() == 1;
  }

  NumericMatrix projected(dims, n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    int center = closest_vertex[cell] - 1;
    if (center < 0 || center >= n_centers) center = 0;

    std::vector<double> best(dims, 0.0);
    double best_distance = std::numeric_limits<double>::infinity();
    for (int neighbor : center_adj[center]) {
      double ab_squared = 0.0;
      double ap_dot_ab = 0.0;
      for (int r = 0; r < dims; ++r) {
        const double ab = y(r, neighbor) - y(r, center);
        const double ap = z(r, cell) - y(r, center);
        ab_squared += ab * ab;
        ap_dot_ab += ap * ab;
      }
      double t = 0.0;
      if (ab_squared > 0.0) {
        t = ap_dot_ab / ab_squared;
      }
      if (!is_tip[center]) {
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
      }

      std::vector<double> candidate(dims, 0.0);
      double distance = 0.0;
      for (int r = 0; r < dims; ++r) {
        candidate[r] = y(r, center) + t * (y(r, neighbor) - y(r, center));
        const double diff = z(r, cell) - candidate[r];
        distance += diff * diff;
      }
      distance = std::sqrt(distance);
      if (distance < best_distance) {
        best_distance = distance;
        best = candidate;
      }
    }
    if (!std::isfinite(best_distance)) {
      for (int r = 0; r < dims; ++r) best[r] = y(r, center);
    }
    for (int r = 0; r < dims; ++r) projected(r, cell) = best[r];
  }

  double min_dist = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n_cells; ++i) {
    for (int j = i + 1; j < n_cells; ++j) {
      const double dist = euclidean_distance_columns(projected, i, j);
      if (dist > 0.0 && dist < min_dist) min_dist = dist;
    }
  }
  if (!std::isfinite(min_dist)) min_dist = 0.0;

  std::vector<double> best_weight(n_cells, std::numeric_limits<double>::infinity());
  std::vector<int> parent(n_cells, -1);
  std::vector<int> in_tree(n_cells, 0);
  best_weight[0] = 0.0;

  for (int iter = 0; iter < n_cells; ++iter) {
    int node = -1;
    double node_weight = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_cells; ++i) {
      if (!in_tree[i] && best_weight[i] < node_weight) {
        node = i;
        node_weight = best_weight[i];
      }
    }
    if (node < 0) break;
    in_tree[node] = 1;

    for (int other = 0; other < n_cells; ++other) {
      if (in_tree[other] || other == node) continue;
      double weight = euclidean_distance_columns(projected, node, other);
      weight += min_dist;
      if (weight < best_weight[other]) {
        best_weight[other] = weight;
        parent[other] = node;
      }
    }
  }

  IntegerMatrix mst_edges(std::max(0, n_cells - 1), 2);
  NumericVector mst_weights(std::max(0, n_cells - 1));
  int edge_idx = 0;
  for (int node = 1; node < n_cells; ++node) {
    if (parent[node] < 0) continue;
    mst_edges(edge_idx, 0) = parent[node] + 1;
    mst_edges(edge_idx, 1) = node + 1;
    mst_weights[edge_idx] = best_weight[node];
    ++edge_idx;
  }
  if (edge_idx == 0) {
    mst_edges = IntegerMatrix(0, 2);
    mst_weights = NumericVector(0);
  } else if (edge_idx != n_cells - 1) {
    mst_edges = mst_edges(Range(0, edge_idx - 1), Range(0, 1));
    mst_weights = mst_weights[Range(0, edge_idx - 1)];
  }

  return List::create(
    _["projected"] = projected,
    _["edges"] = mst_edges,
    _["weights"] = mst_weights,
    _["closest_vertex"] = closest_vertex,
    _["min_dist"] = min_dist
  );
}

// [[Rcpp::export]]
int monocle2_select_root_by_state_cpp(
    NumericMatrix coords,
    IntegerVector candidate_cells,
    NumericVector pseudotime,
    IntegerVector closest_vertex,
    bool use_min_pseudotime = false) {
  const int n_total = coords.ncol();
  const int n_candidates = candidate_cells.size();
  if (n_candidates < 1) stop("candidate_cells must not be empty");
  if (pseudotime.size() != n_total) stop("pseudotime length must match coords columns");
  if (closest_vertex.size() != n_total) stop("closest_vertex length must match coords columns");

  if (n_candidates == 1) {
    const int cell = candidate_cells[0] - 1;
    if (cell < 0 || cell >= n_total) stop("candidate_cells contains an invalid index");
    return closest_vertex[cell];
  }

  std::vector<int> cells(n_candidates);
  for (int i = 0; i < n_candidates; ++i) {
    cells[i] = candidate_cells[i] - 1;
    if (cells[i] < 0 || cells[i] >= n_total) {
      stop("candidate_cells contains an invalid index");
    }
  }

  std::vector<double> best_weight(n_candidates, std::numeric_limits<double>::infinity());
  std::vector<int> parent(n_candidates, -1);
  std::vector<int> in_tree(n_candidates, 0);
  best_weight[0] = 0.0;

  for (int iter = 0; iter < n_candidates; ++iter) {
    int node = -1;
    double node_weight = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_candidates; ++i) {
      if (!in_tree[i] && best_weight[i] < node_weight) {
        node = i;
        node_weight = best_weight[i];
      }
    }
    if (node < 0) break;
    in_tree[node] = 1;

    for (int other = 0; other < n_candidates; ++other) {
      if (in_tree[other] || other == node) continue;
      const double weight = euclidean_distance_columns(coords, cells[node], cells[other]);
      if (weight < best_weight[other]) {
        best_weight[other] = weight;
        parent[other] = node;
      }
    }
  }

  std::vector<std::vector<std::pair<int, double>>> adj(n_candidates);
  for (int node = 1; node < n_candidates; ++node) {
    if (parent[node] < 0) continue;
    const double weight = best_weight[node];
    adj[node].push_back(std::make_pair(parent[node], weight));
    adj[parent[node]].push_back(std::make_pair(node, weight));
  }

  auto farthest_from = [&](int start, std::vector<int>* parent_out) {
    std::vector<double> dist(n_candidates, -1.0);
    std::vector<int> par(n_candidates, -1);
    std::vector<int> stack;
    stack.push_back(start);
    dist[start] = 0.0;
    while (!stack.empty()) {
      const int node = stack.back();
      stack.pop_back();
      for (const auto& edge : adj[node]) {
        const int next = edge.first;
        if (dist[next] >= 0.0) continue;
        dist[next] = dist[node] + edge.second;
        par[next] = node;
        stack.push_back(next);
      }
    }
    int farthest = start;
    double farthest_dist = dist[start];
    for (int i = 0; i < n_candidates; ++i) {
      if (dist[i] > farthest_dist) {
        farthest = i;
        farthest_dist = dist[i];
      }
    }
    if (parent_out != nullptr) {
      *parent_out = par;
    }
    return farthest;
  };

  const int endpoint_a = farthest_from(0, nullptr);
  std::vector<int> diameter_parent;
  const int endpoint_b = farthest_from(endpoint_a, &diameter_parent);

  std::vector<int> diameter_path;
  int node = endpoint_b;
  while (node >= 0) {
    diameter_path.push_back(node);
    if (node == endpoint_a) break;
    node = diameter_parent[node];
  }
  if (diameter_path.empty()) {
    diameter_path.push_back(endpoint_b);
  }

  int selected = diameter_path[0];
  double selected_pt = pseudotime[cells[selected]];
  for (int idx : diameter_path) {
    const double pt = pseudotime[cells[idx]];
    if ((use_min_pseudotime && pt < selected_pt) ||
        (!use_min_pseudotime && pt > selected_pt)) {
      selected = idx;
      selected_pt = pt;
    }
  }

  return closest_vertex[cells[selected]];
}
