#include <Rcpp.h>
#include "velocity_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
List paga_connectivities_cpp(IntegerMatrix knn_idx, IntegerVector groups, int n_groups) {
  const int n_cells = groups.size();
  if (knn_idx.nrow() != n_cells) stop("knn_idx rows must match groups length");
  if (n_groups < 1) stop("n_groups must be positive");

  NumericMatrix directed_edges(n_groups, n_groups);
  NumericVector group_sizes(n_groups);
  NumericVector edge_totals(n_groups);
  for (int cell = 0; cell < n_cells; ++cell) {
    const int group = groups[cell] - 1;
    if (group < 0 || group >= n_groups) stop("groups must be 1-based integers within n_groups");
    group_sizes[group] += 1.0;
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    const int from_group = groups[cell] - 1;
    for (int col = 0; col < knn_idx.ncol(); ++col) {
      const int neighbor = knn_idx(cell, col);
      if (neighbor == NA_INTEGER) continue;
      const int neighbor0 = neighbor - 1;
      if (neighbor0 < 0 || neighbor0 >= n_cells || neighbor0 == cell) continue;
      const int to_group = groups[neighbor0] - 1;
      if (to_group < 0 || to_group >= n_groups) continue;
      directed_edges(from_group, to_group) += 1.0;
      edge_totals[from_group] += 1.0;
    }
  }

  NumericMatrix connectivities(n_groups, n_groups);

  const double total_edges = std::accumulate(edge_totals.begin(), edge_totals.end(), 0.0);
  const double denom = std::max(total_edges, 1.0);
  for (int i = 0; i < n_groups; ++i) {
    for (int j = i + 1; j < n_groups; ++j) {
      const double observed = directed_edges(i, j) + directed_edges(j, i);
      if (observed <= 0.0) continue;
      const double expected =
        (edge_totals[i] * group_sizes[j] + edge_totals[j] * group_sizes[i]) / denom;
      double scaled = expected > 0.0 ? observed / expected : 1.0;
      if (scaled > 1.0) scaled = 1.0;
      connectivities(i, j) = scaled;
      connectivities(j, i) = scaled;
    }
  }

  std::vector<scop_util::Edge> edges;
  edges.reserve(n_groups * (n_groups - 1) / 2);
  for (int i = 0; i < n_groups; ++i) {
    for (int j = i + 1; j < n_groups; ++j) {
      double w = connectivities(i, j);
      if (w > 0.0 && std::isfinite(w)) edges.push_back({ i, j, w });
    }
  }
  std::sort(edges.begin(), edges.end(), [](const scop_util::Edge& a, const scop_util::Edge& b) {
    if (a.weight > b.weight) return true;
    if (a.weight < b.weight) return false;
    if (a.i < b.i) return true;
    if (a.i > b.i) return false;
    return a.j < b.j;
  });

  NumericMatrix connectivities_tree = scop_util::build_mst_matrix(edges, n_groups);

  return List::create(
    _["connectivities"] = connectivities,
    _["connectivities_tree"] = connectivities_tree,
    _["group_sizes"] = group_sizes
  );
}

// ── 2. PAGA diffusion pseudotime (improved: uses R eigen() for accuracy) ──

// [[Rcpp::export]]
List paga_diffusion_pseudotime_cpp(
    NumericMatrix connectivities,
    IntegerVector root_group,
    int n_dcs = 10,
    int n_branchings = 0,
    NumericVector group_sizes = NumericVector(),
    double min_group_size = 0.01)
{
  const int n_groups = connectivities.nrow();
  if (connectivities.ncol() != n_groups)
    stop("connectivities must be square");

  std::vector<double> row_sum(n_groups, 0.0);
  for (int i = 0; i < n_groups; ++i)
    for (int j = 0; j < n_groups; ++j)
      row_sum[i] += connectivities(i, j);

  // Transition matrix (column-major)
  std::vector<double> T(n_groups * n_groups, 0.0);
  for (int i = 0; i < n_groups; ++i) {
    if (row_sum[i] < 1e-10) {
      row_sum[i] = 1.0;
      T[i * n_groups + i] = 1.0;
    }
  }
  for (int i = 0; i < n_groups; ++i)
    for (int j = 0; j < n_groups; ++j)
      if (connectivities(i, j) > 0)
        T[i + j * n_groups] = connectivities(i, j) / row_sum[i];

  // sqrt(D) normalized symmetric matrix
  NumericMatrix Tsym(n_groups, n_groups);
  for (int i = 0; i < n_groups; ++i) {
    double di = row_sum[i] > 0 ? std::sqrt(row_sum[i]) : 1.0;
    for (int j = 0; j < n_groups; ++j) {
      double dj = row_sum[j] > 0 ? std::sqrt(row_sum[j]) : 1.0;
      Tsym(i, j) = di * T[i + j * n_groups] / dj;
    }
  }

  // Use R's eigen() for accurate eigendecomposition (replaces manual power iteration)
  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(Tsym, Named("symmetric", true));
  NumericVector evals_c = eig["values"];
  NumericMatrix evecs_c = eig["vectors"];

  // Sort by eigenvalue magnitude (descending)
  std::vector<std::pair<double, int>> pairs;
  for (int i = 0; i < n_groups; ++i)
    pairs.push_back({std::abs(evals_c[i]), i});
  std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

  int k = std::min(n_dcs + 1, n_groups);
  std::vector<int> use_idx;
  for (int i = 0; i < k; ++i) {
    double ev = evals_c[pairs[i].second];
    if (std::abs(ev - 1.0) > 0.01 && ev > 0.01)
      use_idx.push_back(pairs[i].second);
  }
  int n_use = (int)use_idx.size();
  if (n_use < 1) n_use = 1;
  if (n_use > n_dcs) n_use = n_dcs;
  if ((int)use_idx.size() > n_use) use_idx.resize(n_use);
  if (use_idx.empty() && pairs.size() > 1) use_idx.push_back(pairs[1].second);

  NumericMatrix dc(n_groups, n_use);
  NumericVector eigvals_out(k);
  for (int i = 0; i < k; ++i) {
    int orig = pairs[i].second;
    eigvals_out[i] = evals_c[orig];
  }
  for (int g = 0; g < n_groups; ++g) {
    for (int d = 0; d < n_use; ++d) {
      int comp = use_idx[d];
      double di = row_sum[g] > 0 ? 1.0 / std::sqrt(row_sum[g]) : 1.0;
      dc(g, d) = di * evecs_c(g, comp);
    }
  }

  // Branching detection (if n_branchings > 0)
  int n_branches_found = 0;
  if (n_branchings > 0 && n_use >= 2) {
    // Detect branches by finding groups far from the main trajectory
    // in the second+ diffusion component
    for (int d = 1; d < n_use && n_branches_found < n_branchings; ++d) {
      std::vector<double> comp_vals(n_groups);
      for (int g = 0; g < n_groups; ++g) comp_vals[g] = dc(g, d);
      // Check if this component separates groups into distinct branches
      double mean_val = 0;
      for (int g = 0; g < n_groups; ++g) mean_val += comp_vals[g];
      mean_val /= n_groups;
      // A branch is detected if the component has large spread
      double spread = 0;
      for (int g = 0; g < n_groups; ++g) spread += (comp_vals[g] - mean_val) * (comp_vals[g] - mean_val);
      spread /= n_groups;
      if (spread > 0.01 * (1.0 / n_groups)) n_branches_found++;
    }
  }

  NumericVector pseudotime(n_groups);
  if (min_group_size > 0 && group_sizes.size() == n_groups) {
    double total_cells = 0.0;
    for (int g = 0; g < n_groups; ++g) total_cells += group_sizes[g];
    double min_cells = total_cells * min_group_size;
    for (int g = 0; g < n_groups; ++g)
      if (group_sizes[g] < min_cells) pseudotime[g] = NA_REAL;
  }
  int root = 0;
  for (int i = 0; i < root_group.size(); ++i) {
    int rg = root_group[i] - 1;
    if (rg >= 0 && rg < n_groups) { root = rg; break; }
  }
  for (int g = 0; g < n_groups; ++g) {
    double dist = 0.0;
    for (int d = 0; d < n_use; ++d)
      dist += (dc(g, d) - dc(root, d)) * (dc(g, d) - dc(root, d));
    pseudotime[g] = std::sqrt(dist);
  }
  double pmax = 0.0;
  for (int g = 0; g < n_groups; ++g)
    if (pseudotime[g] > pmax) pmax = pseudotime[g];
  if (pmax > 0)
    for (int g = 0; g < n_groups; ++g) pseudotime[g] /= pmax;

  return List::create(
    _["pseudotime"] = pseudotime,
    _["diffusion_components"] = dc,
    _["diffusion_eigenvalues"] = eigvals_out,
    _["root_group"] = root + 1,
    _["n_branchings_found"] = n_branches_found
  );
}

// ── 3. PAGA velocity transitions (group-level) ───────────────────────────────

// [[Rcpp::export]]
List paga_velocity_transitions_cpp(
    NumericMatrix velocity_embedding,  // cells x dims
    IntegerMatrix knn_idx,             // cells x k (1-based)
    IntegerVector groups,              // 1-based
    int n_groups,
    double softmax_scale = 4.0)
{
  const int n_cells = velocity_embedding.nrow();
  const int n_dims = velocity_embedding.ncol();
  const int n_neighbors = knn_idx.ncol();

  if (groups.size() != n_cells)
    stop("groups length must match number of cells");

  // Build group-level transition matrix
  NumericMatrix transitions(n_groups, n_groups);
  NumericVector group_sizes(n_groups);
  for (int i = 0; i < n_cells; ++i) {
    int g = groups[i] - 1;
    if (g < 0 || g >= n_groups) continue;
    group_sizes[g] += 1.0;
  }

  for (int cell = 0; cell < n_cells; ++cell) {
    int from_group = groups[cell] - 1;
    if (from_group < 0 || from_group >= n_groups) continue;

    double vn = 0.0;
    for (int d = 0; d < n_dims; ++d)
      vn += velocity_embedding(cell, d) * velocity_embedding(cell, d);
    vn = std::sqrt(vn);
    if (vn < 1e-10) continue;

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      int to_group = groups[nb] - 1;
      if (to_group < 0 || to_group >= n_groups) continue;

      double dot = 0.0, dn = 0.0;
      for (int d = 0; d < n_dims; ++d) {
        double delta = velocity_embedding(nb, d) - velocity_embedding(cell, d);
        dot += velocity_embedding(cell, d) * delta;
        dn += delta * delta;
      }
      dn = std::sqrt(dn);
      if (dn < 1e-10) continue;

      double cosine = dot / (vn * dn);
      if (cosine > 0 && std::isfinite(cosine)) {
        double weight = std::exp(cosine * softmax_scale);
        transitions(from_group, to_group) += weight;
        row_sum += weight;
      }
    }
  }

  // Normalize rows
  for (int i = 0; i < n_groups; ++i) {
    double rs = 0.0;
    for (int j = 0; j < n_groups; ++j) rs += transitions(i, j);
    if (rs > 0)
      for (int j = 0; j < n_groups; ++j) transitions(i, j) /= rs;
  }

  // Build MST of transitions for tree
  std::vector<scop_util::Edge> edges;
  edges.reserve(n_groups * (n_groups - 1) / 2);
  for (int i = 0; i < n_groups; ++i)
    for (int j = i + 1; j < n_groups; ++j) {
      double w = (transitions(i, j) + transitions(j, i)) / 2.0;
      if (w > 0) edges.push_back({i, j, w});
    }
  std::sort(edges.begin(), edges.end(), [](const scop_util::Edge& a, const scop_util::Edge& b) {
    return a.weight > b.weight;
  });

  NumericMatrix transitions_tree = scop_util::build_mst_matrix(edges, n_groups);

  return List::create(
    _["transitions_confidence"] = transitions,
    _["transitions_confidence_tree"] = transitions_tree,
    _["group_sizes"] = group_sizes
  );
}

// ── 4. PAGA root cell selection ────────────────────────────────────────────────

// [[Rcpp::export]]
IntegerVector paga_root_cell_cpp(
    NumericMatrix embedding,    // cells x dims (e.g., UMAP)
    IntegerVector groups,       // 1-based
    int root_group)             // which group is root (1-based)
{
  const int n_cells = embedding.nrow();
  const int n_dims = embedding.ncol();

  // Find centroid of root group
  std::vector<double> centroid(n_dims, 0.0);
  int count = 0;
  for (int i = 0; i < n_cells; ++i) {
    if (groups[i] == root_group) {
      for (int d = 0; d < n_dims; ++d)
        centroid[d] += embedding(i, d);
      ++count;
    }
  }
  if (count == 0) count = 1;
  for (int d = 0; d < n_dims; ++d) centroid[d] /= count;

  // Find cell in root_group closest to centroid
  int best_cell = 0;
  double best_dist = std::numeric_limits<double>::max();
  for (int i = 0; i < n_cells; ++i) {
    if (groups[i] != root_group) continue;
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d) {
      double delta = embedding(i, d) - centroid[d];
      dist += delta * delta;
    }
    if (dist < best_dist) {
      best_dist = dist;
      best_cell = i;
    }
  }

  // Return all cells in root_group sorted by distance to centroid (1-based)
  std::vector<std::pair<double, int>> dist_idx;
  for (int i = 0; i < n_cells; ++i) {
    if (groups[i] != root_group) continue;
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d) {
      double delta = embedding(i, d) - centroid[d];
      dist += delta * delta;
    }
    dist_idx.push_back({dist, i + 1});  // 1-based
  }
  std::sort(dist_idx.begin(), dist_idx.end());

  IntegerVector result(dist_idx.size());
  for (size_t i = 0; i < dist_idx.size(); ++i)
    result[i] = dist_idx[i].second;

  return result;  // 1-based cell indices, closest to centroid first
}
