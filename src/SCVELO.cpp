#include <Rcpp.h>
#include "velocity_utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// ── 1. Filter genes ───────────────────────────────────────────────────────────

// [[Rcpp::export]]
IntegerVector scvelo_filter_genes_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced,
    int min_counts = 3,
    int min_counts_u = 3)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");

  IntegerVector keep(n_genes, 1);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 16)
  #endif
  for (int g = 0; g < n_genes; ++g) {
    double sum_s = 0.0, sum_u = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      sum_s += spliced(g, c);
      sum_u += unspliced(g, c);
    }
    if (sum_s < min_counts || sum_u < min_counts_u)
      keep[g] = 0;
  }
  return keep;
}


// ── 2. Normalize per cell + log1p ─────────────────────────────────────────────

// [[Rcpp::export]]
List scvelo_normalize_log_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");

  NumericMatrix ns(n_genes, n_cells);
  NumericMatrix nu(n_genes, n_cells);

  for (int c = 0; c < n_cells; ++c) {
    double cell_sum = 0.0;
    for (int g = 0; g < n_genes; ++g)
      cell_sum += spliced(g, c) + unspliced(g, c);
    double scale = cell_sum > 0 ? 1e4 / cell_sum : 1.0;
    for (int g = 0; g < n_genes; ++g) {
      ns(g, c) = std::log1p(spliced(g, c) * scale);
      nu(g, c) = std::log1p(unspliced(g, c) * scale);
    }
  }
  return List::create(_["spliced_norm"] = ns, _["unspliced_norm"] = nu);
}


// ── 3. Compute moments (KNN smoothing) ────────────────────────────────────────

// [[Rcpp::export]]
List scvelo_moments_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced,
    IntegerMatrix knn_idx)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  const int n_neighbors = knn_idx.ncol();

  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");

  NumericMatrix Ms(n_genes, n_cells);
  NumericMatrix Mu(n_genes, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    int count = 1;
    for (int gene = 0; gene < n_genes; ++gene) {
      Ms(gene, cell) = spliced(gene, cell);
      Mu(gene, cell) = unspliced(gene, cell);
    }
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;
      ++count;
      for (int gene = 0; gene < n_genes; ++gene) {
        Ms(gene, cell) += spliced(gene, nb);
        Mu(gene, cell) += unspliced(gene, nb);
      }
    }
    double inv = 1.0 / count;
    for (int gene = 0; gene < n_genes; ++gene) {
      Ms(gene, cell) *= inv;
      Mu(gene, cell) *= inv;
    }
  }
  return List::create(_["Ms"] = Ms, _["Mu"] = Mu);
}

// [[Rcpp::export]]
List scvelo_moments_connectivities_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced,
    IntegerMatrix knn_idx)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  const int n_neighbors = knn_idx.ncol();

  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");

  std::vector<std::vector<int> > adj(n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    adj[cell].push_back(cell);
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;
      adj[cell].push_back(nb);
      adj[nb].push_back(cell);
    }
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    std::sort(adj[cell].begin(), adj[cell].end());
    adj[cell].erase(std::unique(adj[cell].begin(), adj[cell].end()), adj[cell].end());
  }
  std::vector<int> offsets(n_cells + 1, 0);
  for (int cell = 0; cell < n_cells; ++cell)
    offsets[cell + 1] = offsets[cell] + static_cast<int>(adj[cell].size());
  std::vector<int> neighbors(offsets[n_cells]);
  for (int cell = 0; cell < n_cells; ++cell)
    std::copy(adj[cell].begin(), adj[cell].end(), neighbors.begin() + offsets[cell]);

  NumericMatrix Ms(n_genes, n_cells);
  NumericMatrix Mu(n_genes, n_cells);
  NumericMatrix Mss(n_genes, n_cells);
  NumericMatrix Mus(n_genes, n_cells);
  const double* spliced_ptr = REAL(spliced);
  const double* unspliced_ptr = REAL(unspliced);
  double* Ms_ptr = REAL(Ms);
  double* Mu_ptr = REAL(Mu);
  double* Mss_ptr = REAL(Mss);
  double* Mus_ptr = REAL(Mus);

  for (int cell = 0; cell < n_cells; ++cell) {
    const int start = offsets[cell];
    const int end = offsets[cell + 1];
    const double inv = start == end ? 1.0 : 1.0 / static_cast<double>(end - start);
    const int out_base = cell * n_genes;
    for (int pos = start; pos < end; ++pos) {
      const int nb = neighbors[pos];
      const int in_base = nb * n_genes;
      for (int gene = 0; gene < n_genes; ++gene) {
        const int in_idx = in_base + gene;
        const int out_idx = out_base + gene;
        const double s = spliced_ptr[in_idx];
        const double u = unspliced_ptr[in_idx];
        Ms_ptr[out_idx] += s * inv;
        Mu_ptr[out_idx] += u * inv;
        Mss_ptr[out_idx] += s * s * inv;
        Mus_ptr[out_idx] += s * u * inv;
      }
    }
  }

  return List::create(_["Ms"] = Ms, _["Mu"] = Mu, _["Mss"] = Mss, _["Mus"] = Mus);
}

// [[Rcpp::export]]
List scvelo_second_order_moments_cpp(
    NumericMatrix spliced,
    NumericMatrix unspliced,
    IntegerMatrix knn_idx)
{
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  const int n_neighbors = knn_idx.ncol();

  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells)
    stop("spliced and unspliced must have identical dimensions");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");

  std::vector<std::vector<int> > adj(n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    adj[cell].push_back(cell);
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;
      adj[cell].push_back(nb);
      adj[nb].push_back(cell);
    }
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    std::sort(adj[cell].begin(), adj[cell].end());
    adj[cell].erase(std::unique(adj[cell].begin(), adj[cell].end()), adj[cell].end());
  }

  NumericMatrix Mss(n_genes, n_cells);
  NumericMatrix Mus(n_genes, n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    const double inv = adj[cell].empty() ? 1.0 : 1.0 / static_cast<double>(adj[cell].size());
    for (int nb : adj[cell]) {
      for (int gene = 0; gene < n_genes; ++gene) {
        const double s = spliced(gene, nb);
        Mss(gene, cell) += s * s * inv;
        Mus(gene, cell) += s * unspliced(gene, nb) * inv;
      }
    }
  }

  return List::create(_["Mss"] = Mss, _["Mus"] = Mus);
}


// ── 4. Deterministic velocity + embedding ─────────────────────────────────────

// [[Rcpp::export]]
List scvelo_deterministic_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    IntegerMatrix knn_idx,
    NumericMatrix embedding,
    bool fit_offset = false,
    double perc = 0.0)
{
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  const int n_dims = embedding.ncol();

  if (Mu.nrow() != n_genes || Mu.ncol() != n_cells)
    stop("Ms and Mu must have identical dimensions");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");
  if (embedding.nrow() != n_cells)
    stop("embedding rows must match number of cells");

  // --- Estimate gamma per gene ---
  NumericVector gamma(n_genes);
  NumericVector gamma_r2(n_genes);
  IntegerVector velocity_genes(n_genes, 1);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 16)
  #endif
  for (int g = 0; g < n_genes; ++g) {
    // Collect valid (s, u) pairs
    std::vector<double> s_vals, u_vals;
    s_vals.reserve(n_cells);
    u_vals.reserve(n_cells);
    for (int c = 0; c < n_cells; ++c) {
      double s = Ms(g, c);
      double u = Mu(g, c);
      if (s > 0 && u >= 0) {
        s_vals.push_back(s);
        u_vals.push_back(u);
      }
    }
    int nv = (int)s_vals.size();
    if (nv < 2) { gamma[g] = 0.0; gamma_r2[g] = 0.0; velocity_genes[g] = 0; continue; }

    // Quantile filtering if perc > 0: match Python parity_bridge.py
    // Python: norm = s/max(s) + u/max(u), keep top perc%
    if (perc > 0.0 && perc < 100.0) {
      // Compute normalized (s/u) metric
      double s_max = 0.0, u_max = 0.0;
      for (int i = 0; i < nv; ++i) {
        if (s_vals[i] > s_max) s_max = s_vals[i];
        if (u_vals[i] > u_max) u_max = u_vals[i];
      }
      if (s_max < 1e-10) s_max = 1e-10;
      if (u_max < 1e-10) u_max = 1e-10;

      std::vector<double> norm(nv);
      for (int i = 0; i < nv; ++i)
        norm[i] = s_vals[i] / s_max + u_vals[i] / u_max;

      // Compute percentile cutoff
      std::vector<double> norm_sorted = norm;
      std::sort(norm_sorted.begin(), norm_sorted.end());
      int cutoff_idx = (int)(nv * perc / 100.0);
      if (cutoff_idx >= nv) cutoff_idx = nv - 1;
      double cutoff = norm_sorted[cutoff_idx];

      // Keep only cells with norm >= cutoff
      std::vector<double> sk, uk;
      for (int i = 0; i < nv; ++i) {
        if (norm[i] >= cutoff) {
          sk.push_back(s_vals[i]);
          uk.push_back(u_vals[i]);
        }
      }
      // Fallback: if too few, use all
      if (sk.size() < 3) {
        sk = s_vals;
        uk = u_vals;
      }
      s_vals = sk;
      u_vals = uk;
      nv = (int)s_vals.size();
    }

    // OLS through origin (or with offset) on (potentially filtered) data
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int used = nv;
    for (int i = 0; i < nv; ++i) {
      double x = s_vals[i];
      double y = u_vals[i];
      sx += x; sy += y; sxx += x*x; sxy += x*y;
    }
    if (used < 2 || sxx < 1e-12) { gamma[g] = 0.0; gamma_r2[g] = 0.0; velocity_genes[g] = 0; continue; }

    if (fit_offset) {
      double mx = sx / used, my = sy / used;
      double var_x = sxx / used - mx * mx;
      if (var_x <= 0) { gamma[g] = 0.0; gamma_r2[g] = 0.0; velocity_genes[g] = 0; continue; }
      gamma[g] = (sxy / used - mx * my) / var_x;
    } else {
      gamma[g] = sxy / sxx;
    }
    if (!std::isfinite(gamma[g]) || gamma[g] < 0.0) { gamma[g] = 0.0; }

    // R² on all original data
    double ss_res = 0.0, ss_tot = 0.0;
    double mean_u = sy / used;
    for (int i = 0; i < nv; ++i) {
      double pred = gamma[g] * s_vals[i];
      ss_res += (u_vals[i] - pred) * (u_vals[i] - pred);
      ss_tot += (u_vals[i] - mean_u) * (u_vals[i] - mean_u);
    }
    gamma_r2[g] = ss_tot > 1e-12 ? 1.0 - ss_res / ss_tot : 0.0;
  }

  // Residual velocity
  NumericMatrix residual(n_genes, n_cells);
  for (int c = 0; c < n_cells; ++c) {
    for (int g = 0; g < n_genes; ++g) {
      residual(g, c) = Mu(g, c) - gamma[g] * Ms(g, c);
    }
  }

  // Velocity embedding (cosine projection)
  NumericMatrix velocity_embedding(n_cells, n_dims);
  NumericVector confidence(n_cells);
  NumericVector velocity_length(n_cells);
  scop_util::cosine_projection_embedding(residual, Ms, knn_idx, embedding, velocity_embedding, confidence, velocity_length);
  scop_util::velocity_confidence_py(residual, knn_idx, confidence, velocity_length);

  return List::create(
    _["velocity_embedding"] = velocity_embedding,
    _["confidence"] = confidence,
    _["velocity_length"] = velocity_length,
    _["gamma"] = gamma,
    _["r2"] = gamma_r2,
    _["velocity_genes"] = velocity_genes,
    _["residual"] = residual
  );
}


// ── 5. Stochastic velocity + embedding (refactored, takes Ms/Mu as input) ──

// [[Rcpp::export]]
List scvelo_stochastic_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    NumericMatrix Mss,
    NumericMatrix Mus,
    IntegerMatrix knn_idx,
    NumericMatrix embedding)
{
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  const int n_dims = embedding.ncol();

  if (Mu.nrow() != n_genes || Mu.ncol() != n_cells)
    stop("Ms and Mu must have identical dimensions");
  if (Mss.nrow() != n_genes || Mss.ncol() != n_cells)
    stop("Mss dimensions must match Ms");
  if (Mus.nrow() != n_genes || Mus.ncol() != n_cells)
    stop("Mus dimensions must match Ms");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");
  if (embedding.nrow() != n_cells)
    stop("embedding rows must match number of cells");

  const double* Ms_ptr = REAL(Ms);
  const double* Mu_ptr = REAL(Mu);
  const double* Mss_ptr = REAL(Mss);
  const double* Mus_ptr = REAL(Mus);

  NumericVector gamma(n_genes);
  NumericVector gamma_r2(n_genes);
  IntegerVector velocity_genes(n_genes, 1);
  std::vector<double> deterministic_res_std(n_genes, 1.0);
  std::vector<char> stochastic_update_gene(n_genes, 1);
  for (int g = 0; g < n_genes; ++g) {
    double num_det = 0.0, den_det = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double s = Ms_ptr[idx];
      double u = Mu_ptr[idx];
      if (!std::isfinite(s) || !std::isfinite(u)) continue;
      num_det += s * u;
      den_det += s * s;
    }
    const double gamma_det = den_det > 1e-12 ? num_det / den_det : 0.0;
    gamma[g] = std::isfinite(gamma_det) && gamma_det > 0.0 ? gamma_det : 0.0;

    double res_mean = 0.0;
    int n_res = 0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double r = Mu_ptr[idx] - gamma_det * Ms_ptr[idx];
      if (!std::isfinite(r)) continue;
      res_mean += r;
      ++n_res;
    }
    res_mean = n_res > 0 ? res_mean / static_cast<double>(n_res) : 0.0;
    double res_var = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double r = Mu_ptr[idx] - gamma_det * Ms_ptr[idx];
      if (!std::isfinite(r)) continue;
      res_var += (r - res_mean) * (r - res_mean);
    }
    double res_std = n_res > 0 ? std::sqrt(res_var / static_cast<double>(n_res)) : 1.0;
    if (std::isfinite(res_std) && res_std < 1e-8) {
      stochastic_update_gene[g] = 0;
      gamma_r2[g] = 1.0;
      velocity_genes[g] = 1;
      deterministic_res_std[g] = 1.0;
      continue;
    }
    if (!std::isfinite(res_std) || res_std < 1e-8) res_std = 1.0;
    deterministic_res_std[g] = res_std;

    double mean_u = 0.0;
    int n_u = 0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double u = Mu_ptr[idx];
      if (!std::isfinite(u)) continue;
      mean_u += u;
      ++n_u;
    }
    mean_u = n_u > 0 ? mean_u / static_cast<double>(n_u) : 0.0;
    double ss_res = 0.0, ss_tot = 0.0;
    bool has_ms = false, has_mu = false;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double s = Ms_ptr[idx];
      double u = Mu_ptr[idx];
      if (std::isfinite(s) && s > 0.0) has_ms = true;
      if (std::isfinite(u) && u > 0.0) has_mu = true;
      if (!std::isfinite(s) || !std::isfinite(u)) continue;
      double r = u - gamma[g] * s;
      ss_res += r * r;
      ss_tot += (u - mean_u) * (u - mean_u);
    }
    gamma_r2[g] = ss_tot > 1e-12 ? 1.0 - ss_res / ss_tot : 0.0;
    if (!std::isfinite(gamma_r2[g])) gamma_r2[g] = 0.0;
    velocity_genes[g] = (gamma_r2[g] > 0.01 && gamma[g] > 0.01 && has_ms && has_mu) ? 1 : 0;
  }

  int n_velocity_genes = 0;
  for (int g = 0; g < n_genes; ++g) n_velocity_genes += velocity_genes[g] > 0 ? 1 : 0;
  if (n_velocity_genes < 2 && n_genes > 0) {
    std::vector<double> r2_sorted(n_genes);
    for (int g = 0; g < n_genes; ++g) r2_sorted[g] = gamma_r2[g];
    std::sort(r2_sorted.begin(), r2_sorted.end());
    int idx80 = static_cast<int>(std::floor(0.80 * static_cast<double>(n_genes - 1)));
    if (idx80 < 0) idx80 = 0;
    if (idx80 >= n_genes) idx80 = n_genes - 1;
    const double min_r2 = r2_sorted[idx80];
    n_velocity_genes = 0;
    for (int g = 0; g < n_genes; ++g) {
      velocity_genes[g] = gamma_r2[g] > min_r2 ? 1 : 0;
      n_velocity_genes += velocity_genes[g] > 0 ? 1 : 0;
    }
  }

  // Match scvelo: stochastic generalized fit updates gamma only for
  // deterministic velocity genes; non-selected genes keep deterministic gamma.
  for (int g = 0; g < n_genes; ++g) {
    if (velocity_genes[g] == 0 || !stochastic_update_gene[g]) continue;

    double num_g2 = 0.0, den_g2 = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double x = 2.0 * Mss_ptr[idx] - Ms_ptr[idx];
      double y = 2.0 * Mus_ptr[idx] + Mu_ptr[idx];
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      num_g2 += x * y;
      den_g2 += x * x;
    }
    const double gamma2 = den_g2 > 1e-12 ? num_g2 / den_g2 : 0.0;

    double res2_mean = 0.0;
    int n_res2 = 0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double var_ss = 2.0 * Mss_ptr[idx] - Ms_ptr[idx];
      double cov_us = 2.0 * Mus_ptr[idx] + Mu_ptr[idx];
      double r2 = cov_us - gamma2 * var_ss;
      if (!std::isfinite(r2)) continue;
      res2_mean += r2;
      ++n_res2;
    }
    res2_mean = n_res2 > 0 ? res2_mean / static_cast<double>(n_res2) : 0.0;
    double res2_var = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      double var_ss = 2.0 * Mss_ptr[idx] - Ms_ptr[idx];
      double cov_us = 2.0 * Mus_ptr[idx] + Mu_ptr[idx];
      double r2 = cov_us - gamma2 * var_ss;
      if (!std::isfinite(r2)) continue;
      res2_var += (r2 - res2_mean) * (r2 - res2_mean);
    }
    double res2_std = n_res2 > 0 ? std::sqrt(res2_var / static_cast<double>(n_res2)) : 1.0;
    if (!std::isfinite(res2_std) || res2_std < 1e-8) res2_std = 1.0;
    const double res_std = deterministic_res_std[g];

    double num = 0.0, den = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      const double var_ss = 2.0 * Mss_ptr[idx] - Ms_ptr[idx];
      const double cov_us = 2.0 * Mus_ptr[idx] + Mu_ptr[idx];
      const double x1 = Ms_ptr[idx] / res_std;
      const double y1 = Mu_ptr[idx] / res_std;
      const double x2 = var_ss / res2_std;
      const double y2 = cov_us / res2_std;
      if (std::isfinite(x1) && std::isfinite(y1)) {
        num += x1 * y1;
        den += x1 * x1;
      }
      if (std::isfinite(x2) && std::isfinite(y2)) {
        num += x2 * y2;
        den += x2 * x2;
      }
    }
    gamma[g] = den > 1e-12 ? num / den : gamma[g];
    if (!std::isfinite(gamma[g]) || gamma[g] < 0.0) gamma[g] = 0.0;
  }

  NumericMatrix residual(n_genes, n_cells);
  double* residual_ptr = REAL(residual);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int g = 0; g < n_genes; ++g) {
    for (int c = 0; c < n_cells; ++c) {
      const int idx = c * n_genes + g;
      residual_ptr[idx] = Mu_ptr[idx] - gamma[g] * Ms_ptr[idx];
    }
  }

  // Velocity embedding (cosine projection)
  NumericMatrix velocity_embedding(n_cells, n_dims);
  NumericVector confidence(n_cells);
  NumericVector velocity_length(n_cells);
  scop_util::cosine_projection_embedding(residual, Ms, knn_idx, embedding, velocity_embedding, confidence, velocity_length);
  scop_util::velocity_confidence_py(residual, knn_idx, confidence, velocity_length);

  return List::create(
    _["velocity_embedding"] = velocity_embedding,
    _["confidence"] = confidence,
    _["velocity_length"] = velocity_length,
    _["gamma"] = gamma,
    _["r2"] = gamma_r2,
    _["velocity_genes"] = velocity_genes,
    _["residual"] = residual
  );
}


// ── 6. Velocity graph (cosine similarity on gene space) ───────────────────────

// [[Rcpp::export]]
List scvelo_velocity_graph_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    NumericMatrix residual,    // gene × cell velocity residuals
    IntegerMatrix knn_idx,     // cells × k (1-based)
    int n_neighbors_velo = -1,
    double softmax_scale = 4.0,
    bool sqrt_transform = false,
    int n_recurse_neighbors = 1)
{
  (void)softmax_scale;
  (void)Mu;
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;
  if (n_recurse_neighbors < 1) n_recurse_neighbors = 1;

  std::vector<int> rows, cols, rows_neg, cols_neg;
  std::vector<double> vals, vals_neg;
  NumericMatrix velocity(n_genes, n_cells);
  std::vector<double> velocity_norm(n_cells, 0.0);

  for (int cell = 0; cell < n_cells; ++cell) {
    double mean_v = 0.0;
    for (int g = 0; g < n_genes; ++g) {
      double v = residual(g, cell);
      if (sqrt_transform) v = std::sqrt(std::abs(v)) * (v < 0.0 ? -1.0 : 1.0);
      velocity(g, cell) = v;
      mean_v += v;
    }
    mean_v /= static_cast<double>(n_genes);
    double norm = 0.0;
    for (int g = 0; g < n_genes; ++g) {
      velocity(g, cell) -= mean_v;
      norm += velocity(g, cell) * velocity(g, cell);
    }
    velocity_norm[cell] = std::sqrt(norm);
  }

  auto add_neighbor = [&](std::vector<int>& out, std::vector<char>& seen, int nb) {
    if (nb == NA_INTEGER) return;
    nb -= 1;
    if (nb < 0 || nb >= n_cells || seen[nb]) return;
    seen[nb] = 1;
    out.push_back(nb);
  };

  // Reusable buffers for neighbor collection to avoid per-cell allocations
  std::vector<char> seen_buf(n_cells, 0);
  std::vector<int> current_buf, all_buf;
  current_buf.reserve(n_cells);
  all_buf.reserve(n_cells);

  auto collect_neighbors = [&](int cell) -> const std::vector<int>& {
    current_buf.clear();
    all_buf.clear();
    std::fill(seen_buf.begin(), seen_buf.end(), 0);
    seen_buf[cell] = 1;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int before = static_cast<int>(all_buf.size());
      add_neighbor(all_buf, seen_buf, knn_idx(cell, col));
      if (static_cast<int>(all_buf.size()) > before) current_buf.push_back(all_buf.back());
    }
    for (int depth = 1; depth < n_recurse_neighbors; ++depth) {
      std::vector<int> next;
      for (int parent : current_buf) {
        for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
          int before = static_cast<int>(all_buf.size());
          add_neighbor(all_buf, seen_buf, knn_idx(parent, col));
          if (static_cast<int>(all_buf.size()) > before) next.push_back(all_buf.back());
        }
      }
      current_buf.swap(next);
      if (current_buf.empty()) break;
    }
    return all_buf;
  };

  for (int cell = 0; cell < n_cells; ++cell) {
    double vn = velocity_norm[cell];
    if (vn < 1e-10) continue;

    const std::vector<int>& neighs = collect_neighbors(cell);
    static thread_local std::vector<double> delta;
    if (delta.size() != static_cast<size_t>(n_genes)) delta.resize(n_genes);
    for (int nb : neighs) {
      double mean_delta = 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double d = Ms(g, nb) - Ms(g, cell);
        if (sqrt_transform) d = std::sqrt(std::abs(d)) * (d < 0.0 ? -1.0 : 1.0);
        delta[g] = d;
        mean_delta += d;
      }
      mean_delta /= static_cast<double>(n_genes);

      double dot = 0.0, dn = 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double d = delta[g] - mean_delta;
        dot += velocity(g, cell) * d;
        dn += d * d;
      }
      dn = std::sqrt(dn);
      if (dn < 1e-10) continue;

      double cosine = dot / (vn * dn);
      if (cosine > 0 && std::isfinite(cosine)) {
        double value = std::min(1.0, cosine);
        rows.push_back(cell);
        cols.push_back(nb);
        vals.push_back(value);
      } else if (cosine < 0 && std::isfinite(cosine)) {
        double value = std::max(-1.0, cosine);
        rows_neg.push_back(cell);
        cols_neg.push_back(nb);
        vals_neg.push_back(value);
      }
    }
  }

  return List::create(
    _["velocity_graph_rows"] = IntegerVector(rows.begin(), rows.end()),
    _["velocity_graph_cols"] = IntegerVector(cols.begin(), cols.end()),
    _["velocity_graph_vals"] = NumericVector(vals.begin(), vals.end()),
    _["velocity_graph_neg_rows"] = IntegerVector(rows_neg.begin(), rows_neg.end()),
    _["velocity_graph_neg_cols"] = IntegerVector(cols_neg.begin(), cols_neg.end()),
    _["velocity_graph_neg_vals"] = NumericVector(vals_neg.begin(), vals_neg.end())
  );
}

// [[Rcpp::export]]
NumericMatrix scvelo_project_velocity_embedding_cpp(
    IntegerVector graph_rows,
    IntegerVector graph_cols,
    NumericVector graph_vals,
    IntegerVector graph_neg_rows,
    IntegerVector graph_neg_cols,
    NumericVector graph_neg_vals,
    NumericMatrix embedding,
    double scale = 10.0,
    bool self_transitions = true,
    bool use_negative_cosines = true)
{
  const int n_cells = embedding.nrow();
  const int n_dims = embedding.ncol();
  std::vector<std::vector<std::pair<int, double> > > rows(n_cells);
  std::vector<double> max_conf(n_cells, 0.0);

  for (int k = 0; k < graph_rows.size(); ++k) {
    int i = graph_rows[k];
    int j = graph_cols[k];
    if (i < 0 || i >= n_cells || j < 0 || j >= n_cells) continue;
    double v = graph_vals[k];
    if (!std::isfinite(v)) continue;
    max_conf[i] = std::max(max_conf[i], v);
    rows[i].push_back(std::make_pair(j, std::expm1(v * scale)));
  }
  for (int k = 0; k < graph_neg_rows.size(); ++k) {
    int i = graph_neg_rows[k];
    int j = graph_neg_cols[k];
    if (i < 0 || i >= n_cells || j < 0 || j >= n_cells) continue;
    double v = graph_neg_vals[k];
    if (!std::isfinite(v)) continue;
    double tv = use_negative_cosines ? -std::expm1(-v * scale) : std::expm1(v * scale) + 1.0;
    rows[i].push_back(std::make_pair(j, tv));
  }

  if (self_transitions && n_cells > 0) {
    std::vector<double> conf_sorted = max_conf;
    std::sort(conf_sorted.begin(), conf_sorted.end());
    double pos = 0.98 * static_cast<double>(n_cells - 1);
    int lo = static_cast<int>(std::floor(pos));
    int hi = static_cast<int>(std::ceil(pos));
    double frac = pos - static_cast<double>(lo);
    double ub = conf_sorted[lo] * (1.0 - frac) + conf_sorted[hi] * frac;
    for (int i = 0; i < n_cells; ++i) {
      double self_prob = std::min(1.0, std::max(0.0, ub - max_conf[i]));
      if (self_prob != 0.0) rows[i].push_back(std::make_pair(i, self_prob));
    }
  }

  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (const auto& entry : rows[i]) row_sum += std::abs(entry.second);
    if (std::abs(row_sum) > 1e-15) {
      for (auto& entry : rows[i]) entry.second /= row_sum;
    }
  }

  NumericMatrix velocity_embedding(n_cells, n_dims);
  for (int i = 0; i < n_cells; ++i) {
    double prob_sum = 0.0;
    int used = 0;
    std::vector<double> dx_sum(n_dims, 0.0);
    for (const auto& entry : rows[i]) {
      int j = entry.first;
      double p = entry.second;
      if (j == i || !std::isfinite(p) || p == 0.0) continue;
      double norm = 0.0;
      for (int d = 0; d < n_dims; ++d) {
        double delta = embedding(j, d) - embedding(i, d);
        norm += delta * delta;
      }
      norm = std::sqrt(norm);
      if (norm <= 1e-15 || !std::isfinite(norm)) continue;
      for (int d = 0; d < n_dims; ++d) {
        double ndx = (embedding(j, d) - embedding(i, d)) / norm;
        velocity_embedding(i, d) += p * ndx;
        dx_sum[d] += ndx;
      }
      prob_sum += p;
      ++used;
    }
    if (used > 0) {
      double mean_prob = prob_sum / static_cast<double>(used);
      for (int d = 0; d < n_dims; ++d) {
        velocity_embedding(i, d) -= mean_prob * dx_sum[d];
      }
    }
  }

  return velocity_embedding;
}

// ── 7. Velocity confidence metrics ─────────────────────────────────────────────

// [[Rcpp::export]]
List scvelo_velocity_confidence_cpp(
    NumericMatrix Ms,
    NumericMatrix residual,   // gene × cell velocity residuals
    IntegerMatrix knn_idx)    // cells × k (1-based)
{
  const int n_cells = Ms.ncol();

  NumericVector confidence(n_cells);
  NumericVector confidence_diff(n_cells);
  NumericVector velocity_length(n_cells);

  scop_util::velocity_confidence_py(residual, knn_idx, confidence, velocity_length);
  for (int cell = 0; cell < n_cells; ++cell) confidence_diff[cell] = confidence[cell];

  return List::create(
    _["confidence"] = confidence,
    _["confidence_diff"] = confidence_diff,
    _["velocity_length"] = velocity_length
  );
}

// ── 8. Terminal states (root_cells, end_points via Markov eigenvectors) ───────

// [[Rcpp::export]]
NumericMatrix scvelo_velocity_transition_cpp(
    NumericMatrix Ms,
    NumericMatrix residual,
    IntegerMatrix knn_idx,
    int n_neighbors_velo = -1,
    double softmax_scale = 10.0)
{
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (residual.nrow() != n_genes || residual.ncol() != n_cells)
    stop("Ms and residual must have identical dimensions");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match number of cells");
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;

  NumericMatrix T(n_cells, n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    double vn = 0.0;
    for (int g = 0; g < n_genes; ++g)
      vn += residual(g, cell) * residual(g, cell);
    vn = std::sqrt(vn);
    if (vn < 1e-10) {
      T(cell, cell) = 1.0;
      continue;
    }

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      double dot = 0.0, dn = 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double delta = Ms(g, nb) - Ms(g, cell);
        dot += residual(g, cell) * delta;
        dn += delta * delta;
      }
      dn = std::sqrt(dn);
      if (dn < 1e-10) continue;

      double cosine = dot / (vn * dn);
      if (std::isfinite(cosine)) {
        double value = 0.0;
        if (cosine > 0.0) {
          value = softmax_scale > 0.0 ? std::expm1(cosine * softmax_scale) : cosine;
        } else if (cosine < 0.0) {
          value = softmax_scale > 0.0 ? std::exp(cosine * softmax_scale) : 0.0;
        }
        if (value <= 0.0) continue;
        T(cell, nb) = value;
        row_sum += value;
      }
    }
    if (row_sum > 1e-12) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }
  return T;
}

static NumericMatrix scvelo_backward_transition(const NumericMatrix& T_forward) {
  const int n_cells = T_forward.nrow();
  NumericMatrix T_backward(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) {
      double value = T_forward(j, i);
      T_backward(i, j) = value;
      row_sum += value;
    }
    if (row_sum > 1e-12) {
      for (int j = 0; j < n_cells; ++j) T_backward(i, j) /= row_sum;
    } else {
      T_backward(i, i) = 1.0;
    }
  }
  return T_backward;
}

static NumericVector scvelo_smooth_connectivities(
    NumericVector score,
    IntegerMatrix knn_idx)
{
  const int n_cells = knn_idx.nrow();
  std::vector<std::vector<int> > adj(n_cells);
  for (int i = 0; i < n_cells; ++i) adj[i].push_back(i);
  for (int i = 0; i < n_cells; ++i) {
    for (int col = 0; col < knn_idx.ncol(); ++col) {
      int nb = knn_idx(i, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == i) continue;
      adj[i].push_back(nb);
    }
  }
  NumericVector smoothed(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
    adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
    for (int nb : adj[i]) smoothed[i] += score[nb];
    if (!adj[i].empty()) smoothed[i] /= static_cast<double>(adj[i].size());
  }
  return smoothed;
}

static NumericVector scvelo_clip_scale(NumericVector x) {
  const int n_cells = x.size();
  std::vector<double> sorted;
  sorted.reserve(n_cells);
  for (int i = 0; i < n_cells; ++i)
    if (std::isfinite(x[i])) sorted.push_back(x[i]);
  if (sorted.empty()) return NumericVector(n_cells);
  std::sort(sorted.begin(), sorted.end());
  int idx = std::min(static_cast<int>(sorted.size()) - 1,
                     static_cast<int>(std::floor(0.98 * (sorted.size() - 1))));
  double upper = sorted[idx];
  NumericVector out(n_cells);
  double xmin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n_cells; ++i) {
    double value = std::isfinite(x[i]) ? std::max(0.0, std::min(x[i], upper)) : 0.0;
    out[i] = value;
    if (value < xmin) xmin = value;
    if (value > xmax) xmax = value;
  }
  double range = xmax - xmin;
  if (range > 1e-12) {
    for (int i = 0; i < n_cells; ++i) out[i] = (out[i] - xmin) / range;
  }
  return out;
}

// [[Rcpp::export]]
List scvelo_terminal_states_transition_cpp(
    NumericMatrix transition_matrix,
    IntegerMatrix knn_idx)
{
  const int n_cells = transition_matrix.nrow();
  if (transition_matrix.ncol() != n_cells)
    stop("transition_matrix must be square");
  if (knn_idx.nrow() != n_cells)
    stop("knn_idx rows must match transition_matrix");

  NumericMatrix T_forward(clone(transition_matrix));
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) row_sum += T_forward(i, j);
    if (row_sum > 1e-12) {
      for (int j = 0; j < n_cells; ++j) T_forward(i, j) /= row_sum;
    } else {
      T_forward(i, i) = 1.0;
    }
  }

  NumericMatrix T_backward = scvelo_backward_transition(T_forward);
  NumericVector roots_raw = scop_util::stationary_distribution(T_backward, 1000, 1e-10);
  NumericVector ends_raw = scop_util::stationary_distribution(T_forward, 1000, 1e-10);

  NumericVector root_cells = scvelo_clip_scale(scvelo_smooth_connectivities(roots_raw, knn_idx));
  NumericVector end_points = scvelo_clip_scale(scvelo_smooth_connectivities(ends_raw, knn_idx));

  int n_root = 0, n_end = 0;
  for (int i = 0; i < n_cells; ++i) {
    if (root_cells[i] >= 0.95) ++n_root;
    if (end_points[i] >= 0.95) ++n_end;
  }

  return List::create(
    _["root_cells"] = root_cells,
    _["end_points"] = end_points,
    _["n_root_regions"] = n_root,
    _["n_end_regions"] = n_end
  );
}

// [[Rcpp::export]]
List scvelo_terminal_states_cpp(
    NumericMatrix velocity_embedding,
    NumericMatrix embedding,
    IntegerMatrix knn_idx,
    int n_neighbors_velo = 10,
    int seed = 0)
{
  const int n_cells = embedding.nrow();
  const int n_dims = embedding.ncol();

  if (velocity_embedding.nrow() != n_cells || velocity_embedding.ncol() != n_dims)
    stop("velocity_embedding must have same dimensions as embedding");

  // Build velocity transition matrix (embedding-space cosine)
  std::vector<double> T;
  scop_util::build_velocity_transition(velocity_embedding, embedding, knn_idx, n_neighbors_velo, T);

  NumericMatrix T_forward(n_cells, n_cells);
  NumericMatrix T_backward(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) {
      double value = T[i + j * n_cells];
      T_forward(i, j) = value;
      row_sum += value;
    }
    if (row_sum <= 1e-12) {
      T_forward(i, i) = 1.0;
    }
  }
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) {
      double value = T[j + i * n_cells];
      T_backward(i, j) = value;
      row_sum += value;
    }
    if (row_sum > 1e-12) {
      for (int j = 0; j < n_cells; ++j) T_backward(i, j) /= row_sum;
    } else {
      T_backward(i, i) = 1.0;
    }
  }

  NumericVector roots_raw = scop_util::stationary_distribution(T_backward, 1000, 1e-10);
  NumericVector ends_raw = scop_util::stationary_distribution(T_forward, 1000, 1e-10);

  NumericVector root_cells = scvelo_clip_scale(scvelo_smooth_connectivities(roots_raw, knn_idx));
  NumericVector end_points = scvelo_clip_scale(scvelo_smooth_connectivities(ends_raw, knn_idx));

  int n_root = 0, n_end = 0;
  for (int i = 0; i < n_cells; ++i) {
    if (root_cells[i] >= 0.95) ++n_root;
    if (end_points[i] >= 0.95) ++n_end;
  }

  return List::create(
    _["root_cells"] = root_cells,
    _["end_points"] = end_points,
    _["n_root_regions"] = n_root,
    _["n_end_regions"] = n_end
  );
}

static NumericMatrix scvelo_graph_transition_matrix(
    IntegerVector graph_rows,
    IntegerVector graph_cols,
    NumericVector graph_vals,
    IntegerVector graph_neg_rows,
    IntegerVector graph_neg_cols,
    NumericVector graph_neg_vals,
    int n_cells,
    bool backward,
    double scale = 10.0)
{
  NumericMatrix T(n_cells, n_cells);
  for (int k = 0; k < graph_rows.size(); ++k) {
    int i = graph_rows[k], j = graph_cols[k];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      T(i, j) += std::expm1(graph_vals[k] * scale);
    }
  }
  for (int k = 0; k < graph_neg_rows.size(); ++k) {
    int i = graph_neg_rows[k], j = graph_neg_cols[k];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      T(i, j) += std::exp(graph_neg_vals[k] * scale);
    }
  }
  if (backward) {
    NumericMatrix Tb(n_cells, n_cells);
    for (int i = 0; i < n_cells; ++i)
      for (int j = 0; j < n_cells; ++j)
        Tb(i, j) = T(j, i);
    T = Tb;
  }
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) row_sum += T(i, j);
    if (row_sum > 1e-12 && std::isfinite(row_sum)) {
      for (int j = 0; j < n_cells; ++j) T(i, j) /= row_sum;
    }
  }
  return T;
}

static double scvelo_percentile(std::vector<double> x, double pct) {
  if (x.empty()) return 0.0;
  std::sort(x.begin(), x.end());
  double pos = (pct / 100.0) * static_cast<double>(x.size() - 1);
  int lo = static_cast<int>(std::floor(pos));
  int hi = static_cast<int>(std::ceil(pos));
  if (lo == hi) return x[lo];
  double w = pos - lo;
  return x[lo] * (1.0 - w) + x[hi] * w;
}

static NumericMatrix scvelo_terminal_eigvecs(NumericMatrix T, double eps) {
  const int n_cells = T.nrow();
  NumericMatrix TT(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_cells; ++j)
      TT(i, j) = T(j, i);

  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(TT, Named("symmetric", false));
  ComplexVector evals = eig["values"];
  ComplexMatrix evecs = eig["vectors"];

  std::vector<std::pair<double, int>> selected;
  for (int i = 0; i < evals.size(); ++i) {
    double real_eval = evals[i].r;
    if (std::isfinite(real_eval) && real_eval >= 1.0 - eps)
      selected.push_back({real_eval, i});
  }
  std::sort(selected.begin(), selected.end(), std::greater<std::pair<double,int>>());

  NumericMatrix out(n_cells, selected.size());
  for (int comp = 0; comp < static_cast<int>(selected.size()); ++comp) {
    int idx = selected[comp].second;
    std::vector<double> values;
    values.reserve(n_cells);
    for (int i = 0; i < n_cells; ++i)
      values.push_back(std::abs(evecs(i, idx).r));

    double lower = scvelo_percentile(values, 2.0);
    double upper = scvelo_percentile(values, 98.0);
    double vmax = 0.0;
    for (int i = 0; i < n_cells; ++i) {
      double value = values[i] < lower ? 0.0 : std::min(values[i], upper);
      out(i, comp) = value;
      if (value > vmax) vmax = value;
    }
    if (vmax > 0.0) {
      for (int i = 0; i < n_cells; ++i) out(i, comp) /= vmax;
    }
  }
  return out;
}

// [[Rcpp::export]]
List scvelo_terminal_states_graph_cpp(
    IntegerVector graph_rows,
    IntegerVector graph_cols,
    NumericVector graph_vals,
    IntegerVector graph_neg_rows,
    IntegerVector graph_neg_cols,
    NumericVector graph_neg_vals,
    IntegerMatrix knn_idx,
    double eps = 1e-3)
{
  const int n_cells = knn_idx.nrow();

  NumericMatrix T_backward = scvelo_graph_transition_matrix(
    graph_rows, graph_cols, graph_vals,
    graph_neg_rows, graph_neg_cols, graph_neg_vals,
    n_cells, true);
  NumericMatrix root_eig = scvelo_terminal_eigvecs(T_backward, eps);
  NumericVector roots_raw(n_cells);
  for (int comp = 0; comp < root_eig.ncol(); ++comp)
    for (int i = 0; i < n_cells; ++i)
      roots_raw[i] += root_eig(i, comp);
  NumericVector root_cells = scvelo_clip_scale(scvelo_smooth_connectivities(roots_raw, knn_idx));

  NumericMatrix T_forward = scvelo_graph_transition_matrix(
    graph_rows, graph_cols, graph_vals,
    graph_neg_rows, graph_neg_cols, graph_neg_vals,
    n_cells, false);
  NumericMatrix end_eig = scvelo_terminal_eigvecs(T_forward, eps);
  NumericVector ends_raw(n_cells);
  for (int comp = 0; comp < end_eig.ncol(); ++comp)
    for (int i = 0; i < n_cells; ++i)
      ends_raw[i] += end_eig(i, comp);
  NumericVector end_points = scvelo_clip_scale(scvelo_smooth_connectivities(ends_raw, knn_idx));

  return List::create(
    _["root_cells"] = root_cells,
    _["end_points"] = end_points,
    _["n_root_regions"] = root_eig.ncol(),
    _["n_end_regions"] = end_eig.ncol()
  );
}


// ── 9. Velocity pseudotime from transition matrix (DPT via eigendecomposition) ─

// [[Rcpp::export]]
List scvelo_pseudotime_transition_cpp(
    NumericMatrix transition_matrix,
    NumericVector root_cells,
    NumericVector end_points)
{
  const int n_cells = transition_matrix.nrow();
  if (transition_matrix.ncol() != n_cells)
    stop("transition_matrix must be square");
  if (root_cells.size() != n_cells)
    stop("root_cells length must match n_cells");
  if (end_points.size() != n_cells)
    stop("end_points length must match n_cells");

  NumericMatrix T(clone(transition_matrix));
  for (int i = 0; i < n_cells; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n_cells; ++j) row_sum += T(i, j);
    if (row_sum > 1e-12) {
      for (int j = 0; j < n_cells; ++j) T(i, j) /= row_sum;
    } else {
      T(i, i) = 1.0;
    }
  }

  int root = 0;
  double rv = root_cells[0];
  for (int i = 1; i < n_cells; ++i)
    if (root_cells[i] > rv) { rv = root_cells[i]; root = i; }
  int end = 0;
  double ev = end_points[0];
  for (int i = 1; i < n_cells; ++i)
    if (end_points[i] > ev) { ev = end_points[i]; end = i; }

  NumericMatrix D(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      D(i, j) = (T(i, j) + T(j, i)) / 2.0;
    }
  }

  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(D, Named("symmetric", true));
  NumericVector evals_c = eig["values"];
  NumericMatrix evecs_c = eig["vectors"];

  std::vector<std::pair<double, int>> pairs;
  for (int i = 0; i < n_cells; ++i)
    pairs.push_back({evals_c[i], i});
  std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

  int k = std::min(10, n_cells);
  NumericMatrix dc(n_cells, k);
  for (int comp = 0; comp < k; ++comp) {
    int idx = pairs[comp].second;
    for (int i = 0; i < n_cells; ++i)
      dc(i, comp) = evecs_c(i, idx);
  }

  auto dpt_from = [&](int source) {
    NumericVector pt(n_cells);
    for (int i = 0; i < n_cells; ++i) {
      double dist = 0.0;
      for (int d = 0; d < k; ++d)
        dist += (dc(i, d) - dc(source, d)) * (dc(i, d) - dc(source, d));
      pt[i] = std::sqrt(dist);
    }
    double pmin = pt[0], pmax = pt[0];
    for (int i = 1; i < n_cells; ++i) {
      if (pt[i] < pmin) pmin = pt[i];
      if (pt[i] > pmax) pmax = pt[i];
    }
    double range = pmax - pmin;
    if (range > 0)
      for (int i = 0; i < n_cells; ++i) pt[i] = (pt[i] - pmin) / range;
    return pt;
  };

  NumericVector pseudotime_root = dpt_from(root);
  NumericVector pseudotime_end_inverse = dpt_from(end);
  for (int i = 0; i < n_cells; ++i) pseudotime_end_inverse[i] = 1.0 - pseudotime_end_inverse[i];

  NumericVector pseudotime(n_cells);
  for (int i = 0; i < n_cells; ++i)
    pseudotime[i] = 0.5 * (pseudotime_root[i] + pseudotime_end_inverse[i]);

  double pmin = pseudotime[0], pmax = pseudotime[0];
  for (int i = 1; i < n_cells; ++i) {
    if (pseudotime[i] < pmin) pmin = pseudotime[i];
    if (pseudotime[i] > pmax) pmax = pseudotime[i];
  }
  double prange = pmax - pmin;
  if (prange > 0)
    for (int i = 0; i < n_cells; ++i) pseudotime[i] = (pseudotime[i] - pmin) / prange;

  return List::create(
    _["pseudotime"] = pseudotime,
    _["pseudotime_root"] = pseudotime_root,
    _["pseudotime_end_inverse"] = pseudotime_end_inverse,
    _["root_cell"] = root + 1,
    _["end_cell"] = end + 1,
    _["diffusion_components"] = dc
  );
}

// [[Rcpp::export]]
List scvelo_pseudotime_cpp(
    NumericMatrix velocity_embedding,
    NumericMatrix embedding,
    IntegerMatrix knn_idx,
    NumericVector root_cells,
    NumericVector end_points,
    int n_neighbors_velo = 10)
{
  const int n_cells = embedding.nrow();
  const int n_dims = embedding.ncol();

  if (velocity_embedding.nrow() != n_cells || velocity_embedding.ncol() != n_dims)
    stop("velocity_embedding dimensions mismatch");
  if (root_cells.size() != n_cells)
    stop("root_cells length must match n_cells");
  if (end_points.size() != n_cells)
    stop("end_points length must match n_cells");

  // Build velocity transition matrix (embedding-space cosine)
  std::vector<double> T;
  scop_util::build_velocity_transition(velocity_embedding, embedding, knn_idx, n_neighbors_velo, T);

  int root = 0;
  double rv = root_cells[0];
  for (int i = 1; i < n_cells; ++i)
    if (root_cells[i] > rv) { rv = root_cells[i]; root = i; }
  int end = 0;
  double ev = end_points[0];
  for (int i = 1; i < n_cells; ++i)
    if (end_points[i] > ev) { ev = end_points[i]; end = i; }

  // Symmetrize: D = (T + T^T) / 2 for diffusion components
  NumericMatrix D(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      D(i, j) = (T[i + j * n_cells] + T[j + i * n_cells]) / 2.0;
    }
  }

  // Use R's eigen() for accurate eigendecomposition
  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(D, Named("symmetric", true));
  NumericVector evals_c = eig["values"];
  NumericMatrix evecs_c = eig["vectors"];

  // Sort by eigenvalue descending
  std::vector<std::pair<double, int>> pairs;
  for (int i = 0; i < n_cells; ++i)
    pairs.push_back({evals_c[i], i});
  std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

  // Diffusion components for the current C++ DPT approximation.
  int k = std::min(10, n_cells);
  NumericMatrix dc(n_cells, k);
  for (int comp = 0; comp < k; ++comp) {
    int idx = pairs[comp].second;
    for (int i = 0; i < n_cells; ++i)
      dc(i, comp) = evecs_c(i, idx);
  }

  auto dpt_from = [&](int source) {
    NumericVector pt(n_cells);
    for (int i = 0; i < n_cells; ++i) {
      double dist = 0.0;
      for (int d = 0; d < k; ++d)
        dist += (dc(i, d) - dc(source, d)) * (dc(i, d) - dc(source, d));
      pt[i] = std::sqrt(dist);
    }
    double pmin = pt[0], pmax = pt[0];
    for (int i = 1; i < n_cells; ++i) {
      if (pt[i] < pmin) pmin = pt[i];
      if (pt[i] > pmax) pmax = pt[i];
    }
    double range = pmax - pmin;
    if (range > 0)
      for (int i = 0; i < n_cells; ++i) pt[i] = (pt[i] - pmin) / range;
    return pt;
  };

  NumericVector pseudotime_root = dpt_from(root);
  NumericVector pseudotime_end_inverse = dpt_from(end);
  for (int i = 0; i < n_cells; ++i) pseudotime_end_inverse[i] = 1.0 - pseudotime_end_inverse[i];

  NumericVector pseudotime(n_cells);
  for (int i = 0; i < n_cells; ++i)
    pseudotime[i] = 0.5 * (pseudotime_root[i] + pseudotime_end_inverse[i]);

  double pmin = pseudotime[0], pmax = pseudotime[0];
  for (int i = 1; i < n_cells; ++i) {
    if (pseudotime[i] < pmin) pmin = pseudotime[i];
    if (pseudotime[i] > pmax) pmax = pseudotime[i];
  }
  double prange = pmax - pmin;
  if (prange > 0)
    for (int i = 0; i < n_cells; ++i) pseudotime[i] = (pseudotime[i] - pmin) / prange;

  return List::create(
    _["pseudotime"] = pseudotime,
    _["pseudotime_root"] = pseudotime_root,
    _["pseudotime_end_inverse"] = pseudotime_end_inverse,
    _["root_cell"] = root + 1,
    _["end_cell"] = end + 1,
    _["diffusion_components"] = dc
  );
}

// [[Rcpp::export]]
List scvelo_pseudotime_graph_cpp(
    IntegerVector graph_rows,
    IntegerVector graph_cols,
    NumericVector graph_vals,
    IntegerVector graph_neg_rows,
    IntegerVector graph_neg_cols,
    NumericVector graph_neg_vals,
    IntegerMatrix knn_idx,
    NumericVector root_cells,
    NumericVector end_points,
    int n_dcs = 10)
{
  const int n_cells = knn_idx.nrow();
  if (root_cells.size() != n_cells)
    stop("root_cells length must match number of cells");
  if (end_points.size() != n_cells)
    stop("end_points length must match number of cells");

  NumericMatrix C(n_cells, n_cells);
  for (int k = 0; k < graph_rows.size(); ++k) {
    int i = graph_rows[k], j = graph_cols[k];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      C(i, j) += graph_vals[k];
      C(j, i) += graph_vals[k];
    }
  }
  for (int k = 0; k < graph_neg_rows.size(); ++k) {
    int i = graph_neg_rows[k], j = graph_neg_cols[k];
    if (i >= 0 && i < n_cells && j >= 0 && j < n_cells) {
      double value = std::abs(graph_neg_vals[k]);
      C(i, j) += value;
      C(j, i) += value;
    }
  }

  NumericVector q(n_cells);
  for (int j = 0; j < n_cells; ++j) {
    for (int i = 0; i < n_cells; ++i) q[j] += C(i, j);
    if (q[j] <= 0.0 || !std::isfinite(q[j])) q[j] = 1.0;
  }

  NumericMatrix K(n_cells, n_cells);
  NumericVector z(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_cells; ++j) {
      K(i, j) = C(i, j) / (q[i] * q[j]);
      z[j] += K(i, j);
    }
  }
  for (int j = 0; j < n_cells; ++j) {
    z[j] = (z[j] > 0.0 && std::isfinite(z[j])) ? std::sqrt(z[j]) : 1.0;
  }

  NumericMatrix D(n_cells, n_cells);
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_cells; ++j)
      D(i, j) = K(i, j) / (z[i] * z[j]);

  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(D, Named("symmetric", true));
  NumericVector evals_c = eig["values"];
  NumericMatrix evecs_c = eig["vectors"];

  std::vector<std::pair<double, int>> pairs;
  for (int i = 0; i < n_cells; ++i)
    pairs.push_back({evals_c[i], i});
  std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

  int k_use = std::min(std::max(1, n_dcs), n_cells);
  NumericVector evals(k_use);
  NumericMatrix dc(n_cells, k_use);
  for (int comp = 0; comp < k_use; ++comp) {
    int idx = pairs[comp].second;
    evals[comp] = evals_c[idx];
    for (int i = 0; i < n_cells; ++i)
      dc(i, comp) = evecs_c(i, idx);
  }

  auto smoothed_argmax = [&](NumericVector score) {
    NumericVector smoothed = scvelo_smooth_connectivities(score, knn_idx);
    int best = 0;
    double best_val = smoothed[0];
    for (int i = 1; i < n_cells; ++i) {
      if (smoothed[i] > best_val) {
        best_val = smoothed[i];
        best = i;
      }
    }
    return best;
  };

  int root = smoothed_argmax(root_cells);
  int end = smoothed_argmax(end_points);
  auto has_informative_score = [&](NumericVector score) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_cells; ++i) {
      if (!std::isfinite(score[i])) continue;
      if (score[i] < xmin) xmin = score[i];
      if (score[i] > xmax) xmax = score[i];
    }
    return std::isfinite(xmin) && std::isfinite(xmax) && xmax > 0.0 && (xmax - xmin) > 1e-12;
  };
  bool use_root = has_informative_score(root_cells);
  bool use_end = has_informative_score(end_points);

  auto dpt_from = [&](int source) {
    NumericVector pt(n_cells);
    for (int i = 0; i < n_cells; ++i) {
      double dist = 0.0;
      for (int d = 0; d < k_use; ++d) {
        double lambda = evals[d];
        double delta = dc(source, d) - dc(i, d);
        if (lambda < 0.9994) {
          double weight = lambda / (1.0 - lambda);
          dist += std::pow(weight * delta, 2.0);
        } else {
          dist += delta * delta;
        }
      }
      pt[i] = std::sqrt(std::max(0.0, dist));
    }
    double pmax = 0.0;
    for (int i = 0; i < n_cells; ++i)
      if (std::isfinite(pt[i]) && pt[i] > pmax) pmax = pt[i];
    if (pmax > 0.0)
      for (int i = 0; i < n_cells; ++i) pt[i] /= pmax;
    return pt;
  };

  NumericVector pseudotime_root = dpt_from(root);
  NumericVector pseudotime_end_inverse = dpt_from(end);
  for (int i = 0; i < n_cells; ++i)
    pseudotime_end_inverse[i] = 1.0 - pseudotime_end_inverse[i];

  NumericVector pseudotime(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    if (use_root && use_end) {
      pseudotime[i] = 0.5 * (pseudotime_root[i] + pseudotime_end_inverse[i]);
    } else if (use_end) {
      pseudotime[i] = pseudotime_end_inverse[i];
    } else {
      pseudotime[i] = pseudotime_root[i];
    }
  }

  double pmin = std::numeric_limits<double>::infinity();
  double pmax = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n_cells; ++i) {
    if (!std::isfinite(pseudotime[i])) continue;
    if (pseudotime[i] < pmin) pmin = pseudotime[i];
    if (pseudotime[i] > pmax) pmax = pseudotime[i];
  }
  double prange = pmax - pmin;
  if (std::isfinite(prange) && prange > 0.0)
    for (int i = 0; i < n_cells; ++i) pseudotime[i] = (pseudotime[i] - pmin) / prange;

  return List::create(
    _["pseudotime"] = pseudotime,
    _["pseudotime_root"] = pseudotime_root,
    _["pseudotime_end_inverse"] = pseudotime_end_inverse,
    _["root_cell"] = root + 1,
    _["end_cell"] = end + 1,
    _["diffusion_components"] = dc,
    _["diffusion_eigenvalues"] = evals
  );
}


// ── 10. Rank velocity genes (Spearman-like correlation) ──────────────────────

// [[Rcpp::export]]
NumericVector scvelo_velocity_genes_cpp(
    NumericMatrix Ms,
    NumericMatrix velocity)
{
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  if (velocity.nrow() != n_genes || velocity.ncol() != n_cells)
    stop("Ms and velocity must have identical dimensions");

  // For each gene, compute absolute correlation between velocity and Ms
  NumericVector scores(n_genes);
  for (int g = 0; g < n_genes; ++g) {
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    int nv = 0;
    for (int c = 0; c < n_cells; ++c) {
      double x = Ms(g, c);
      double y = velocity(g, c);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      sx += x; sy += y; sxx += x*x; syy += y*y; sxy += x*y;
      ++nv;
    }
    if (nv < 3) { scores[g] = 0.0; continue; }
    double den = std::sqrt((nv*sxx - sx*sx) * (nv*syy - sy*sy));
    double corr = den > 1e-12 ? (nv*sxy - sx*sy) / den : 0.0;
    scores[g] = std::abs(corr);
  }
  return scores;
}


// ── 11. Keep backward-compatible wrapper (same API as before) ─────────────────

// [[Rcpp::export]]
List scvelo_stochastic_embedding_cpp(
  NumericMatrix spliced,
  NumericMatrix unspliced,
  IntegerMatrix knn_idx,
  NumericMatrix embedding
) {
  // Compute moments
  List moments = scvelo_moments_cpp(spliced, unspliced, knn_idx);
  NumericMatrix Ms = moments["Ms"];
  NumericMatrix Mu = moments["Mu"];
  List second = scvelo_second_order_moments_cpp(spliced, unspliced, knn_idx);
  NumericMatrix Mss = second["Mss"];
  NumericMatrix Mus = second["Mus"];
  // Run stochastic embedding
  return scvelo_stochastic_cpp(Ms, Mu, Mss, Mus, knn_idx, embedding);
}
