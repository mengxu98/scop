#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <limits>
#include <string>

using namespace Rcpp;

// ── 1. Filter genes (matching scv.pp.filter_genes exactly) ──────────────────

// [[Rcpp::export]]
IntegerVector scvelo_filter_genes_scanpy_cpp(
    NumericMatrix spliced,     // genes × cells
    NumericMatrix unspliced,
    int min_counts = 3,
    int min_counts_u = 3)
{
  int n_genes = spliced.nrow();
  int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");

  IntegerVector keep(n_genes, 1);
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
List scvelo_normalize_scanpy_cpp(
    NumericMatrix spliced,     // genes × cells, ALREADY FILTERED
    NumericMatrix unspliced,   // genes × cells, ALREADY FILTERED
    NumericVector initial_spliced_totals,  // per-cell totals BEFORE filtering (length = n_cells)
    NumericVector initial_unspliced_totals)
{
  int n_genes = spliced.nrow();
  int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");
  if (initial_spliced_totals.size() != n_cells || initial_unspliced_totals.size() != n_cells)
    stop("initial totals must match n_cells");

  // Median of pre-filtering totals (matching scvelo's get_initial_size)
  std::vector<double> sorted_s(n_cells), sorted_u(n_cells);
  for (int c = 0; c < n_cells; ++c) {
    sorted_s[c] = initial_spliced_totals[c];
    sorted_u[c] = initial_unspliced_totals[c];
  }
  std::sort(sorted_s.begin(), sorted_s.end());
  std::sort(sorted_u.begin(), sorted_u.end());
  double median_s = sorted_s[n_cells / 2];
  double median_u = sorted_u[n_cells / 2];
  if (median_s <= 0) median_s = 1.0;
  if (median_u <= 0) median_u = 1.0;

  NumericMatrix ns(n_genes, n_cells);
  NumericMatrix nu(n_genes, n_cells);

  for (int c = 0; c < n_cells; ++c) {
    double scale_s = initial_spliced_totals[c] > 0
      ? median_s / initial_spliced_totals[c] : 1.0;
    double scale_u = initial_unspliced_totals[c] > 0
      ? median_u / initial_unspliced_totals[c] : 1.0;
    for (int g = 0; g < n_genes; ++g) {
      ns(g, c) = spliced(g, c) * scale_s;
      nu(g, c) = unspliced(g, c) * scale_u;
    }
  }

  return List::create(
    _["spliced_norm"] = ns,
    _["unspliced_norm"] = nu
  );
}

// ── 3. PCA via covariance_eigh (matching scanpy sc.pp.pca with svd_solver='covariance_eigh') ──

// This is DETERMINISTIC — no random seed needed.
// Matches scanpy's pca(n_pcs, zero_center=True) up to eigenvector sign flips.

// [[Rcpp::export]]
List scvelo_pca_scanpy_cpp(
    NumericMatrix X,           // genes × cells — will be transposed to cells×genes
    int n_pcs = 30)
{
  int n_genes = X.nrow();
  int n_cells = X.ncol();

  int k = std::min(std::min(n_pcs, n_cells), n_genes);
  if (k < 1) k = 1;

  // Step 1: Center each gene (column-center in scanpy space: cells×genes)
  // In our matrix X (genes×cells), we center per ROW (per gene)
  NumericMatrix Xc(n_genes, n_cells);
  for (int g = 0; g < n_genes; ++g) {
    double mean_g = 0.0;
    for (int c = 0; c < n_cells; ++c) mean_g += X(g, c);
    mean_g /= n_cells;
    for (int c = 0; c < n_cells; ++c) Xc(g, c) = X(g, c) - mean_g;
  }

// Step 2: Choose method based on n_genes vs n_cells
  // When n_genes >> n_cells, use dual form (Gram matrix = Xc @ Xc^T, n_cells × n_cells)
  NumericMatrix scores(n_cells, k);
  NumericMatrix loadings(n_genes, k);
  NumericVector variance(k);
  NumericVector variance_ratio(k);

  bool use_gram = (n_genes > n_cells);

  if (use_gram) {
    // Dual form: G = Xc @ Xc^T  (n_cells × n_cells)
    NumericMatrix G(n_cells, n_cells);
    for (int i = 0; i < n_cells; ++i) {
      for (int j = i; j < n_cells; ++j) {
        double s = 0.0;
        for (int g = 0; g < n_genes; ++g)
          s += Xc(g, i) * Xc(g, j);
        G(i, j) = s;
        G(j, i) = s;
      }
    }

    Environment base("package:base");
    Function eigen_fun = base["eigen"];
    List eig = eigen_fun(G, Named("symmetric", true));
    NumericVector evals_g = eig["values"];
    NumericMatrix evecs_g = eig["vectors"];

    // Sort by descending eigenvalue
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < n_cells; ++i)
      if (evals_g[i] > 1e-12)
        pairs.push_back({evals_g[i], i});
    std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

    int k_use = std::min(k, (int)pairs.size());
    if (k_use < 1) k_use = 1;

    for (int comp = 0; comp < k_use; ++comp) {
      int orig = pairs[comp].second;
      double eval = evals_g[orig];
      variance[comp] = eval;
      variance_ratio[comp] = 0.0;  // will compute below

      // Score: PCA coordinates = left singular vector * singular value.
      for (int c = 0; c < n_cells; ++c)
        scores(c, comp) = evecs_g(c, orig) * std::sqrt(eval);

      // Loading: Xc^T @ v / sqrt(eval)  — projection back to gene space
      double norm_sq = 0.0;
      double inv_norm = eval > 1e-12 ? 1.0 / std::sqrt(eval) : 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double s = 0.0;
        for (int c = 0; c < n_cells; ++c)
          s += Xc(g, c) * evecs_g(c, orig);
        loadings(g, comp) = s * inv_norm;
        norm_sq += loadings(g, comp) * loadings(g, comp);
      }
      // Normalize loadings to unit length
      if (norm_sq > 0) {
        double inv_n = 1.0 / std::sqrt(norm_sq);
        for (int g = 0; g < n_genes; ++g) loadings(g, comp) *= inv_n;
      }
    }

    // Fill remaining components with 0
    for (int comp = k_use; comp < k; ++comp) {
      variance[comp] = 0.0;
      variance_ratio[comp] = 0.0;
    }

    // Compute variance ratio from eigenvalues
    double total_eval = 0.0;
    for (int i = 0; i < n_cells; ++i)
      if (evals_g[i] > 0) total_eval += evals_g[i];
    if (total_eval > 0) {
      for (int comp = 0; comp < k_use; ++comp)
        variance_ratio[comp] = variance[comp] / total_eval * 100.0;
    }

    k = k_use;
  } else {
    // Primal form: Cov = Xc @ Xc^T / (n_cells - 1)  (n_genes × n_genes)
    NumericMatrix cov(n_genes, n_genes);
    for (int i = 0; i < n_genes; ++i) {
      for (int j = i; j < n_genes; ++j) {
        double s = 0.0;
        for (int c = 0; c < n_cells; ++c)
          s += Xc(i, c) * Xc(j, c);
        s /= (n_cells - 1);
        cov(i, j) = s;
        cov(j, i) = s;
      }
    }

    Environment base("package:base");
    Function eigen_fun = base["eigen"];
    List eig = eigen_fun(cov, Named("symmetric", true));
    NumericVector evals = eig["values"];
    NumericMatrix evecs = eig["vectors"];

    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < n_genes; ++i)
      pairs.push_back({evals[i], i});
    std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

    double total_var = 0.0;
    for (int i = 0; i < n_genes; ++i) {
      double v = evals[i];
      if (v > 0) total_var += v;
    }
    if (total_var <= 0) total_var = 1.0;

    for (int comp = 0; comp < k; ++comp) {
      int orig = pairs[comp].second;
      double eval = evals[orig];
      variance[comp] = eval;
      variance_ratio[comp] = eval / total_var * 100.0;

      for (int g = 0; g < n_genes; ++g)
        loadings(g, comp) = evecs(g, orig);

      for (int c = 0; c < n_cells; ++c) {
        double s = 0.0;
        for (int g = 0; g < n_genes; ++g)
          s += Xc(g, c) * evecs(g, orig);
        scores(c, comp) = s;
      }
    }
  }

  return List::create(
    _["scores"] = scores,
    _["loadings"] = loadings,
    _["variance"] = variance,
    _["variance_ratio"] = variance_ratio,
    _["n_pcs"] = k
  );
}

// ── 4. KNN via exact brute-force Euclidean (deterministic, matching scanpy) ──

// [[Rcpp::export]]
List scvelo_knn_scanpy_cpp(
    NumericMatrix coords,      // cells × dims
    int n_neighbors = 10,
    bool exclude_self = true)
{
  int n_cells = coords.nrow();
  int n_dims = coords.ncol();
  int k_actual = n_neighbors + (exclude_self ? 1 : 0);
  if (k_actual > n_cells) k_actual = n_cells;

  IntegerMatrix idx(n_cells, n_neighbors);
  NumericMatrix dist(n_cells, n_neighbors);

  for (int i = 0; i < n_cells; ++i) {
    // Compute distances to all other cells
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

    // Sort by distance (ascending), tie-break by index for determinism
    std::sort(dists.begin(), dists.end(),
      [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
        if (a.first < b.first) return true;
        if (a.first > b.first) return false;
        return a.second < b.second;  // tie-break by index
      });

    // Store top k
    for (int k = 0; k < n_neighbors && k < (int)dists.size(); ++k) {
      idx(i, k) = dists[k].second + 1;  // 1-based
      dist(i, k) = std::sqrt(dists[k].first);
    }
    // Fill remaining with NA
    for (int k = dists.size(); k < n_neighbors; ++k) {
      idx(i, k) = NA_INTEGER;
      dist(i, k) = NA_REAL;
    }
  }

  return List::create(_["idx"] = idx, _["dist"] = dist);
}

// ── 5. Filter genes with shared counts (min_shared_counts) ──────────────────

// [[Rcpp::export]]
IntegerVector scvelo_filter_genes_shared_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced,
    int min_shared_counts = 30)
{
  int n_genes = spliced.nrow();
  int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");

  IntegerVector keep(n_genes, 1);
  for (int g = 0; g < n_genes; ++g) {
    int shared_count = 0;
    for (int c = 0; c < n_cells; ++c) {
      if (spliced(g, c) >= min_shared_counts && unspliced(g, c) >= min_shared_counts)
        ++shared_count;
    }
    // Keep gene if at least some cells have shared counts
    // scvelo's filter_genes_shared keeps genes where both S and U
    // have >= min_shared_counts in at least one cell
    if (shared_count == 0) keep[g] = 0;
  }
  return keep;
}

// ── 6. Full scanpy-compatible preprocessing pipeline ───────────────────────

// [[Rcpp::export]]
List scvelo_preprocess_scanpy_cpp(
    NumericMatrix spliced,         // genes × cells
    NumericMatrix unspliced,
    int n_pcs = 30,
    int n_neighbors = 10,
    int min_counts = 3,
    int min_counts_u = 3)
{
  int n_genes = spliced.nrow();
  int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");

  // Step 1: Compute initial per-cell totals (before filtering)
  NumericVector initial_spliced_totals(n_cells);
  NumericVector initial_unspliced_totals(n_cells);
  for (int c = 0; c < n_cells; ++c) {
    double ss = 0.0, su = 0.0;
    for (int g = 0; g < n_genes; ++g) {
      ss += spliced(g, c);
      su += unspliced(g, c);
    }
    initial_spliced_totals[c] = ss;
    initial_unspliced_totals[c] = su;
  }

  // Step 2: Filter genes
  IntegerVector keep = scvelo_filter_genes_scanpy_cpp(spliced, unspliced, min_counts, min_counts_u);

  int n_keep = 0;
  for (int g = 0; g < n_genes; ++g) n_keep += keep[g] > 0;
  if (n_keep < 2) stop("Too few genes pass filtering");

  NumericMatrix spliced_f(n_keep, n_cells);
  NumericMatrix unspliced_f(n_keep, n_cells);
  {
    int idx = 0;
    for (int g = 0; g < n_genes; ++g) {
      if (keep[g] > 0) {
        for (int c = 0; c < n_cells; ++c) {
          spliced_f(idx, c) = spliced(g, c);
          unspliced_f(idx, c) = unspliced(g, c);
        }
        ++idx;
      }
    }
  }

  // Step 3: Normalize (NO log1p — matching scv.pp.normalize_per_cell)
  List normed = scvelo_normalize_scanpy_cpp(spliced_f, unspliced_f,
    initial_spliced_totals, initial_unspliced_totals);
  NumericMatrix spliced_n = normed["spliced_norm"];
  NumericMatrix unspliced_n = normed["unspliced_norm"];

  // Step 4: PCA on log1p(normalized spliced) — matching sc.pp.log1p(adata.X) then sc.pp.pca
  NumericMatrix spliced_for_pca(n_keep, n_cells);
  for (int g = 0; g < n_keep; ++g)
    for (int c = 0; c < n_cells; ++c)
      spliced_for_pca(g, c) = std::log1p(spliced_n(g, c));

  List pca = scvelo_pca_scanpy_cpp(spliced_for_pca, n_pcs);
  NumericMatrix pca_scores = pca["scores"];
  int n_pcs_out = pca["n_pcs"];

  // Step 5: KNN in PCA space. scanpy stores n_neighbors - 1 non-self
  // distances when querying the training data itself.
  int knn_nonself = std::max(1, std::min(n_neighbors - 1, n_cells - 1));
  List knn = scvelo_knn_scanpy_cpp(pca_scores, knn_nonself, true);

  return List::create(
    _["spliced_norm"] = spliced_n,
    _["unspliced_norm"] = unspliced_n,
    _["pca_scores"] = pca_scores,
    _["pca_loadings"] = pca["loadings"],
    _["pca_variance"] = pca["variance"],
    _["pca_variance_ratio"] = pca["variance_ratio"],
    _["knn_idx"] = knn["idx"],
    _["knn_dist"] = knn["dist"],
    _["n_genes_filtered"] = n_keep,
    _["n_pcs"] = n_pcs_out
  );
}
