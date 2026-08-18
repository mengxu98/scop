#include <Rcpp.h>
#include <thisutils/log_message.h>
#include "velocity_utils.h"

using namespace Rcpp;

namespace {

constexpr int kMaxDenseDptCells = 4000;

void check_dense_dpt_size(int n_cells) {
  if (n_cells > kMaxDenseDptCells) {
    thisutils::log_message(
      "cell-level DPT currently requires dense O(n^2) memory; reduce the dataset size or skip cell-level DPT",
      "error"
    );
  }
}

NumericMatrix build_gauss_connectivities(
    IntegerMatrix knn_idx,
    NumericMatrix knn_dist)
{
  const int n_cells = knn_idx.nrow();
  const int n_neighbors = knn_idx.ncol();
  check_dense_dpt_size(n_cells);
  if (knn_dist.nrow() != n_cells || knn_dist.ncol() != n_neighbors) {
    thisutils::log_message("knn_idx and knn_dist must have the same dimensions", "error");
  }

  NumericMatrix W(n_cells, n_cells);
  std::fill(W.begin(), W.end(), 0.0);
  std::vector<double> sigma_sq(n_cells, 1e-20);

  for (int i = 0; i < n_cells; ++i) {
    double max_dsq = 0.0;
    for (int col = 0; col < n_neighbors; ++col) {
      const double d = knn_dist(i, col);
      if (R_FINITE(d)) {
        max_dsq = std::max(max_dsq, d * d);
      }
    }
    sigma_sq[i] = std::max(max_dsq / 4.0, 1e-20);
  }

  for (int i = 0; i < n_cells; ++i) {
    const double sigma_i = std::sqrt(sigma_sq[i]);
    for (int col = 0; col < n_neighbors; ++col) {
      int j = knn_idx(i, col);
      if (j == NA_INTEGER) continue;
      j -= 1;
      if (j < 0 || j >= n_cells || j == i) continue;
      const double d = knn_dist(i, col);
      if (!R_FINITE(d)) continue;
      const double dsq = d * d;
      const double sigma_j = std::sqrt(sigma_sq[j]);
      const double den = sigma_sq[i] + sigma_sq[j];
      const double weight =
        std::sqrt(2.0 * sigma_i * sigma_j / den) * std::exp(-dsq / den);
      W(i, j) = weight;
    }
  }

  for (int i = 0; i < n_cells; ++i) {
    for (int col = 0; col < n_neighbors; ++col) {
      int j = knn_idx(i, col);
      if (j == NA_INTEGER) continue;
      j -= 1;
      if (j < 0 || j >= n_cells || j == i) continue;
      if (W(j, i) == 0.0 && W(i, j) > 0.0) {
        W(j, i) = W(i, j);
      }
    }
  }

  return W;
}

NumericMatrix transition_symmetric(const NumericMatrix& W) {
  const int n_cells = W.nrow();
  check_dense_dpt_size(n_cells);
  if (W.ncol() != n_cells) {
    thisutils::log_message("connectivities must be a square matrix", "error");
  }

  std::vector<double> q(n_cells, 0.0);
  for (int j = 0; j < n_cells; ++j) {
    for (int i = 0; i < n_cells; ++i) {
      q[j] += W(i, j);
    }
    if (q[j] < 1e-20) q[j] = 1e-20;
  }

  NumericMatrix K(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      K(i, j) = W(i, j) / (q[i] * q[j]);
    }
  }

  std::vector<double> z(n_cells, 0.0);
  for (int j = 0; j < n_cells; ++j) {
    for (int i = 0; i < n_cells; ++i) {
      z[j] += K(i, j);
    }
    z[j] = std::sqrt(std::max(z[j], 1e-20));
  }

  NumericMatrix Tsym(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      Tsym(i, j) = K(i, j) / (z[i] * z[j]);
    }
  }
  return Tsym;
}

List dpt_from_transition(
    const NumericMatrix& transitions_sym,
    int root_cell,
    int n_dcs)
{
  const int n_cells = transitions_sym.nrow();
  check_dense_dpt_size(n_cells);
  if (root_cell < 1 || root_cell > n_cells) {
    thisutils::log_message("root_cell must be a valid 1-based cell index", "error");
  }

  const int n_comps = std::min(n_dcs, n_cells - 1);
  if (n_comps < 1) {
    thisutils::log_message("n_dcs must be at least 1", "error");
  }

  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(transitions_sym, Named("symmetric", true));
  NumericVector evals_all = eig["values"];
  NumericMatrix evecs_all = eig["vectors"];

  std::vector<std::pair<double, int>> pairs;
  pairs.reserve(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    pairs.push_back({evals_all[i], i});
  }
  std::sort(
    pairs.begin(),
    pairs.end(),
    std::greater<std::pair<double, int>>()
  );

  NumericVector evals(n_comps);
  NumericMatrix evecs(n_cells, n_comps);
  for (int j = 0; j < n_comps; ++j) {
    const int idx = pairs[j].second;
    evals[j] = evals_all[idx];
    for (int i = 0; i < n_cells; ++i) {
      evecs(i, j) = evecs_all(i, idx);
    }
  }

  NumericVector pseudotime(n_cells);
  const int root = root_cell - 1;
  for (int cell = 0; cell < n_cells; ++cell) {
    double row = 0.0;
    for (int j = 0; j < n_comps; ++j) {
      const double ev = evals[j];
      const double diff = evecs(cell, j) - evecs(root, j);
      if (ev < 0.9994) {
        row += std::pow(ev / (1.0 - ev) * diff, 2.0);
      } else {
        row += diff * diff;
      }
    }
    pseudotime[cell] = std::sqrt(row);
  }

  double pmax = 0.0;
  for (int i = 0; i < n_cells; ++i) {
    if (R_FINITE(pseudotime[i]) && pseudotime[i] > pmax) {
      pmax = pseudotime[i];
    }
  }
  if (pmax > 0.0) {
    for (int i = 0; i < n_cells; ++i) {
      if (R_FINITE(pseudotime[i])) {
        pseudotime[i] /= pmax;
      }
    }
  }

  return List::create(
    _["pseudotime"] = pseudotime,
    _["diffusion_components"] = evecs,
    _["diffusion_eigenvalues"] = evals,
    _["root_cell"] = root_cell
  );
}

}  // namespace

// [[Rcpp::export]]
List paga_connectivities_cpp(IntegerMatrix knn_idx, IntegerVector groups, int n_groups) {
  const int n_cells = groups.size();
  if (knn_idx.nrow() != n_cells) thisutils::log_message("knn_idx rows must match groups length", "error");
  if (n_groups < 1) thisutils::log_message("n_groups must be positive", "error");

  NumericMatrix directed_edges(n_groups, n_groups);
  NumericMatrix expected_n_edges_random(n_groups, n_groups);
  NumericVector group_sizes(n_groups);
  NumericVector edge_totals(n_groups);
  for (int cell = 0; cell < n_cells; ++cell) {
    const int group = groups[cell] - 1;
    if (group < 0 || group >= n_groups) thisutils::log_message("groups must be 1-based integers within n_groups", "error");
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

  // Scanpy PAGA v1.2 uses (n_cells - 1) as the random-null denominator:
  // expected = (es[i] * ns[j] + es[j] * ns[i]) / (n - 1)
  const double denom = std::max(static_cast<double>(n_cells - 1), 1.0);
  for (int i = 0; i < n_groups; ++i) {
    for (int j = i + 1; j < n_groups; ++j) {
      const double observed = directed_edges(i, j) + directed_edges(j, i);
      if (observed <= 0.0) continue;
      const double expected =
        (edge_totals[i] * group_sizes[j] + edge_totals[j] * group_sizes[i]) / denom;
      expected_n_edges_random(i, j) = expected;
      expected_n_edges_random(j, i) = expected;
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
    _["expected_n_edges_random"] = expected_n_edges_random,
    _["group_sizes"] = group_sizes,
    _["directed_edges"] = directed_edges
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
    thisutils::log_message("connectivities must be square", "error");

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

// ── 2b. Scanpy-compatible connectivities and cell-level DPT ──

// [[Rcpp::export]]
NumericMatrix gauss_connectivities_cpp(
    IntegerMatrix knn_idx,
    NumericMatrix knn_dist)
{
  return build_gauss_connectivities(knn_idx, knn_dist);
}

// [[Rcpp::export]]
List dpt_from_connectivities_cpp(
    NumericMatrix connectivities,
    int root_cell,
    int n_dcs = 10)
{
  if (connectivities.nrow() < 2) {
    thisutils::log_message("at least two cells are required for cell-level DPT", "error");
  }
  NumericMatrix transitions_sym = transition_symmetric(connectivities);
  return dpt_from_transition(transitions_sym, root_cell, n_dcs);
}

// [[Rcpp::export]]
List cell_dpt_pseudotime_cpp(
    IntegerMatrix knn_idx,
    NumericMatrix knn_dist,
    int root_cell,
    int n_dcs = 10)
{
  NumericMatrix connectivities = build_gauss_connectivities(knn_idx, knn_dist);
  return dpt_from_transition(
    transition_symmetric(connectivities),
    root_cell,
    n_dcs
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
    thisutils::log_message("groups length must match number of cells", "error");

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
