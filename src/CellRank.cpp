// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "velocity_utils.h"

using namespace Rcpp;
using namespace arma;

// ── 1. Transition matrix validation & normalization ──────────────────────────

// [[Rcpp::export]]
List cellrank_validate_transition_matrix_cpp(
    NumericMatrix T_,
    double eps = 1e-10,
    double min_self_loop = 0.01)
{
  int n = T_.nrow();
  if (T_.ncol() != n) stop("Transition matrix must be square");

  NumericMatrix T = clone(T_);
  int nans = 0, negs = 0, zeros = 0, loops = 0;
  scop_util::validate_transition_matrix(T, eps, min_self_loop, &nans, &negs, &zeros, &loops);

  return List::create(
    _["transition_matrix"] = T,
    _["nans_fixed"] = nans, _["negs_clipped"] = negs,
    _["zero_rows_fixed"] = zeros, _["self_loops_added"] = loops
  );
}

// ── 2. Stationary distribution ───────────────────────────────────────────────

// [[Rcpp::export]]
NumericVector cellrank_stationary_distribution_cpp(
    NumericMatrix T_, int max_iter = 1000, double tol = 1e-10)
{
  return scop_util::stationary_distribution(T_, max_iter, tol);
}

// ── 3. GPCCA Schur decomposition (using R's eigen() for robustness) ────────

// [[Rcpp::export]]
List cellrank_schur_cpp(NumericMatrix T_, int n_components = 10)
{
  int n = T_.nrow();
  if (n_components < 2) n_components = 2;
  if (n_components > n) n_components = n;

  // Stationary distribution
  NumericVector pi = cellrank_stationary_distribution_cpp(T_, 200, 1e-8);

  // Build T_bar = diag(sqrt(pi)) * T * diag(1/sqrt(pi))
  NumericMatrix T_bar(n, n);
  for (int i = 0; i < n; ++i) {
    double si = std::sqrt(std::max(pi[i], 1e-15));
    for (int j = 0; j < n; ++j) {
      double sj = std::sqrt(std::max(pi[j], 1e-15));
      T_bar(i, j) = si * T_(i, j) / sj;
    }
  }

  // Use R's eigen() for robustness on all matrix sizes
  Environment base("package:base");
  Function eigen_fun = base["eigen"];
  List eig = eigen_fun(T_bar, Named("symmetric", false));
  ComplexVector evals_c = eig["values"];
  ComplexMatrix evecs_c = eig["vectors"];

  // Sort by eigenvalue magnitude (descending), take real parts
  std::vector<std::pair<double, int>> pairs;
  for (int i = 0; i < n; ++i)
    pairs.push_back({std::abs(evals_c[i].r), i});
  std::sort(pairs.begin(), pairs.end(), std::greater<std::pair<double,int>>());

  int nc = std::min(n_components, n);
  NumericVector eigenvalues(nc);
  NumericMatrix schur_vecs(n, nc);
  for (int j = 0; j < nc; ++j) {
    int orig_j = pairs[j].second;
    eigenvalues[j] = evals_c[orig_j].r;
    for (int i = 0; i < n; ++i)
      schur_vecs(i, j) = evecs_c(i, orig_j).r;
  }

  // Macrostate assignment: largest non-trivial component per row. The leading
  // Perron vector is close to constant and can otherwise collapse all cells into
  // a single macrostate on connectivity kernels.
  IntegerVector macro(n);
  int start_component = (nc > 1 && std::abs(eigenvalues[0] - 1.0) < 1e-6) ? 1 : 0;
  for (int i = 0; i < n; ++i) {
    int best = start_component; double best_val = -1;
    for (int c = start_component; c < nc; ++c) {
      double v = std::abs(schur_vecs(i, c));
      if (v > best_val) { best_val = v; best = c; }
    }
    macro[i] = best + 1;
  }

  return List::create(
    _["eigenvalues"] = eigenvalues, _["schur_vectors"] = schur_vecs,
    _["macrostate_assignment"] = macro, _["n_macrostates"] = nc,
    _["stationary_distribution"] = pi,
    _["method"] = "eigen_fallback"
  );
}

// ── 4. Auto-detect number of macrostates via eigengap ────────────────────────

// [[Rcpp::export]]
int cellrank_auto_n_states_cpp(NumericVector eigenvalues, int min_states = 2, int max_states = 20)
{
  int n = eigenvalues.size();
  if (n <= min_states) return min_states;

  double max_gap = 0.0;
  int best_k = min_states;
  for (int k = min_states; k < std::min(max_states, n - 1); ++k) {
    double gap = std::abs(eigenvalues[k] - eigenvalues[k + 1]);
    if (gap > max_gap) {
      max_gap = gap;
      best_k = k + 1;
    }
  }

  for (int k = min_states; k < std::min(max_states, n - 1); ++k) {
    double denom = std::abs(eigenvalues[k]);
    if (denom > 1e-10) {
      double rel_gap = std::abs(eigenvalues[k] - eigenvalues[k + 1]) / denom;
      if (rel_gap > 0.3) return k + 1;
    }
  }

  return best_k;
}

// ── 5. Velocity kernel ───────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_velocity_kernel_cpp(
    NumericMatrix velocity_embedding,   // cells × dims
    NumericMatrix embedding,            // cells × dims
    IntegerMatrix knn_idx,              // cells × k (1-based)
    bool backward = false,
    double softmax_scale = 4.0,
    int n_neighbors_velo = -1)
{
  const int n_cells = velocity_embedding.nrow();
  const int n_dims = velocity_embedding.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    double vn = 0.0;
    for (int d = 0; d < n_dims; ++d)
      vn += velocity_embedding(cell, d) * velocity_embedding(cell, d);
    vn = std::sqrt(vn);
    if (vn < 1e-10) { T(cell, cell) = 1.0; continue; }

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      double dot = 0.0, dn = 0.0;
      for (int d = 0; d < n_dims; ++d) {
        double delta = embedding(nb, d) - embedding(cell, d);
        if (backward) {
          dot += (-velocity_embedding(cell, d)) * delta;
        } else {
          dot += velocity_embedding(cell, d) * delta;
        }
        dn += delta * delta;
      }
      dn = std::sqrt(dn);
      if (dn < 1e-10) continue;

      double cosine = dot / (vn * dn);
      if (cosine > 0 && std::isfinite(cosine)) {
        double weight = std::exp(cosine * softmax_scale);
        T(cell, nb) = weight;
        row_sum += weight;
      }
    }
    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// ── 5b. Gene-space velocity kernel (matching Python CellRank) ────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_velocity_kernel_gene_cpp(
    NumericMatrix gene_velocity,         // genes × cells
    NumericMatrix expression,            // genes × cells (Ms or counts)
    IntegerMatrix knn_idx,               // cells × k (1-based)
    bool backward = false,
    double softmax_scale = 4.0,
    int n_neighbors_velo = -1)
{
  const int n_genes = gene_velocity.nrow();
  const int n_cells = gene_velocity.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;

  // Precompute velocity norms per cell
  std::vector<double> vn_cells(n_cells, 0.0);
  for (int c = 0; c < n_cells; ++c) {
    double sq = 0.0;
    for (int g = 0; g < n_genes; ++g)
      sq += gene_velocity(g, c) * gene_velocity(g, c);
    vn_cells[c] = std::sqrt(sq);
  }

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    if (vn_cells[cell] < 1e-10) { T(cell, cell) = 1.0; continue; }

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      // Cosine correlation in gene space:
      // cos(v_cell, expr_nb - expr_cell)
      double dot = 0.0, ndsq = 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double delta = expression(g, nb) - expression(g, cell);
        double v = backward ? -gene_velocity(g, cell) : gene_velocity(g, cell);
        dot += v * delta;
        ndsq += delta * delta;
      }
      double nd = std::sqrt(ndsq);
      if (nd < 1e-10) continue;

      double cosine = dot / (vn_cells[cell] * nd);
      if (!std::isfinite(cosine)) continue;
      // CellRank uses exp(cosine / scale) for ALL neighbors
      double weight = std::exp(cosine / softmax_scale);
      T(cell, nb) = weight;
      row_sum += weight;
    }
    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// ── 6. Pseudotime kernel ─────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_pseudotime_kernel_cpp(
    NumericVector pseudotime,
    IntegerMatrix knn_idx,
    NumericVector cell_weights,
    double bandwidth = 1.0,
    bool backward = false)
{
  const int n_cells = pseudotime.size();
  const int n_neighbors = knn_idx.ncol();

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    double pt_cell = pseudotime[cell];
    double w_cell = (cell_weights.size() >= n_cells) ? cell_weights[cell] : 1.0;

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      double pt_nb = pseudotime[nb];
      double w_nb = (cell_weights.size() >= n_cells) ? cell_weights[nb] : 1.0;

      double dt;
      if (backward) {
        dt = pt_cell - pt_nb;
      } else {
        dt = pt_nb - pt_cell;
      }

      if (dt <= 0) continue;

      double weight = w_nb * std::exp(-dt * dt / (2.0 * bandwidth * bandwidth));
      T(cell, nb) = weight;
      row_sum += weight;
    }

    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// ── 7. CytoTRACE kernel ───────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_cytotrace_kernel_cpp(
    NumericVector gene_counts,
    IntegerMatrix knn_idx,
    double bandwidth = 1.0,
    bool backward = false)
{
  const int n_cells = gene_counts.size();
  const int n_neighbors = knn_idx.ncol();

  double gc_min = gene_counts[0], gc_max = gene_counts[0];
  for (int i = 1; i < n_cells; ++i) {
    if (gene_counts[i] < gc_min) gc_min = gene_counts[i];
    if (gene_counts[i] > gc_max) gc_max = gene_counts[i];
  }
  NumericVector cytotrace_score(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    cytotrace_score[i] = (gc_max > gc_min)
      ? (gene_counts[i] - gc_min) / (gc_max - gc_min) : 0.5;
  }

  NumericVector weights(n_cells, 1.0);
  NumericMatrix T = cellrank_pseudotime_kernel_cpp(
    cytotrace_score, knn_idx, weights, bandwidth, backward);

  return T;
}

// ── 8. CFLARE estimator ──────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_cflare_cpp(
    NumericMatrix T_,
    int n_states = 5,
    int max_iter = 1000,
    double tol = 1e-6)
{
  int n = T_.nrow();
  if (n_states > n) n_states = n;

  List validated = cellrank_validate_transition_matrix_cpp(T_);
  NumericMatrix T = validated["transition_matrix"];

  NumericVector pi = cellrank_stationary_distribution_cpp(T, max_iter, tol);

  List schur = cellrank_schur_cpp(T, n_states);
  NumericVector eigenvalues = schur["eigenvalues"];
  NumericMatrix schur_vecs = schur["schur_vectors"];
  IntegerVector macro = schur["macrostate_assignment"];
  int M = schur["n_macrostates"];

  // Build coarse transition matrix
  NumericMatrix P_coarse(M, M);
  NumericVector sizes(M);
  for (int i = 0; i < n; ++i) sizes[macro[i] - 1] += 1.0;
  for (int i = 0; i < n; ++i) {
    int ai = macro[i] - 1;
    for (int j = 0; j < n; ++j) {
      int aj = macro[j] - 1;
      P_coarse(ai, aj) += T(i, j);
    }
  }
  for (int a = 0; a < M; ++a) {
    double rs = 0.0;
    for (int b = 0; b < M; ++b) rs += P_coarse(a, b);
    if (rs > 0) for (int b = 0; b < M; ++b) P_coarse(a, b) /= rs;
  }

  std::vector<int> is_terminal, term_idx, trans_idx;
  int n_terminal;
  scop_util::detect_terminal_states(P_coarse, M, is_terminal, term_idx, trans_idx, n_terminal);
  int nT = term_idx.size();
  int nQ = trans_idx.size();

  NumericMatrix abs_prob(n, std::max(nT, 1));
  IntegerVector term_state(n, NA_INTEGER);
  NumericVector fate_conf(n, 0.0);
  scop_util::compute_absorption_probabilities(P_coarse, macro, n, M, is_terminal, term_idx, trans_idx, nT, nQ, abs_prob, term_state, fate_conf);
  IntegerVector terminal_obs(n);
  for (int i = 0; i < n; ++i) {
    int m = macro[i] - 1;
    terminal_obs[i] = (m >= 0 && m < M && is_terminal[m]) ? (m + 1) : 0;
  }

  return List::create(
    _["transition_matrix"] = T,
    _["stationary_distribution"] = pi,
    _["eigenvalues"] = eigenvalues,
    _["schur_vectors"] = schur_vecs,
    _["macrostate_assignment"] = macro,
    _["n_macrostates"] = M,
    _["terminal_states"] = terminal_obs,
    _["lineage_assignment"] = term_state,
    _["fate_confidence"] = fate_conf,
    _["absorption_probabilities"] = abs_prob,
    _["n_terminal_states"] = nT,
    _["method"] = "CFLARE"
  );
}

// ── 9. GPCCA Estimator ────────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_gpcca_cpp(
    NumericMatrix T_,
    int n_states = 5,
    int n_cells_terminal = 10,
    bool skip_perron = false)
{
  int n = T_.nrow();

  List validated = cellrank_validate_transition_matrix_cpp(T_);
  NumericMatrix T = validated["transition_matrix"];

  NumericVector pi = cellrank_stationary_distribution_cpp(T, 200, 1e-8);

  int schur_components = skip_perron ? (n_states + 2) : n_states;
  List schur_result = cellrank_schur_cpp(T, schur_components);
  NumericVector eigenvalues_all = schur_result["eigenvalues"];
  NumericMatrix schur_vecs_all = schur_result["schur_vectors"];
  IntegerVector macro_all = schur_result["macrostate_assignment"];

  // If skip_perron, drop the first component (Perron vector, eigenvalue ≈ 1)
  int start_comp = skip_perron ? 1 : 0;
  int M = std::min(n_states, schur_vecs_all.ncol() - start_comp);
  NumericVector eigenvalues(M);
  NumericMatrix schur_vecs(n, M);
  for (int j = 0; j < M; ++j) {
    eigenvalues[j] = eigenvalues_all[j + start_comp];
    for (int i = 0; i < n; ++i)
      schur_vecs(i, j) = schur_vecs_all(i, j + start_comp);
  }
  // Recompute macro from the selected Schur vectors
  IntegerVector macro(n);
  for (int i = 0; i < n; ++i) {
    int best = 0; double best_val = -1;
    for (int c = 0; c < M; ++c) {
      double v = std::abs(schur_vecs(i, c));
      if (v > best_val) { best_val = v; best = c; }
    }
    macro[i] = best + 1;
  }

  // G-PCCA approximation: use Schur vectors as soft macrostate memberships.
  // CellRank uses the real Schur vectors directly (with sign convention), but
  // our Armadillo schur vectors may have different sign conventions per column.
  // Using abs() is more robust and empirically gives better parity with
  // CellRank's fate probabilities than max(v,0) with sign normalization.
  arma::mat chi(n, M);
  chi.zeros();
  for (int i = 0; i < n; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < M; ++j) {
      double v = std::abs(schur_vecs(i, j));
      if (!std::isfinite(v)) v = 0.0;
      chi(i, j) = v;
      row_sum += v;
    }
    if (row_sum <= 1e-10) {
      int m = macro[i] - 1;
      if (m >= 0 && m < M) {
        chi(i, m) = 1.0;
        row_sum = 1.0;
      }
    }
    if (row_sum > 1e-10) {
      for (int j = 0; j < M; ++j) chi(i, j) /= row_sum;
    } else {
      for (int j = 0; j < M; ++j) chi(i, j) = 1.0 / M;
    }

    int best = 0;
    double best_val = chi(i, 0);
    for (int j = 1; j < M; ++j) {
      if (chi(i, j) > best_val) {
        best_val = chi(i, j);
        best = j;
      }
    }
    macro[i] = best + 1;
  }

  // P_coarse = chi^T * diag(pi) * T * chi
  arma::vec pi_arma(n);
  for (int i = 0; i < n; ++i) pi_arma(i) = pi[i];
  arma::mat D_pi = arma::diagmat(pi_arma);
  arma::mat T_arma(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      T_arma(i, j) = T(i, j);

  arma::mat P_coarse_arma = chi.t() * D_pi * T_arma * chi;
  for (int i = 0; i < M; ++i) {
    double rs = arma::accu(P_coarse_arma.row(i));
    if (rs > 1e-12) {
      P_coarse_arma.row(i) /= rs;
    } else {
      P_coarse_arma(i, i) = 1.0;
    }
  }
  NumericMatrix P_coarse = Rcpp::wrap(P_coarse_arma);

  std::vector<int> is_terminal, term_idx, trans_idx;
  int n_terminal;
  scop_util::detect_terminal_states(P_coarse, M, is_terminal, term_idx, trans_idx, n_terminal);
  int nT = term_idx.size();
  int nQ = trans_idx.size();

  NumericMatrix abs_prob(n, std::max(nT, 1));
  IntegerVector term_state(n, NA_INTEGER);
  NumericVector fate_conf(n, 0.0);
  scop_util::compute_absorption_probabilities(P_coarse, macro, n, M, is_terminal, term_idx, trans_idx, nT, nQ, abs_prob, term_state, fate_conf);
  IntegerVector terminal_obs(n);
  for (int i = 0; i < n; ++i) {
    int m = macro[i] - 1;
    terminal_obs[i] = (m >= 0 && m < M && is_terminal[m]) ? (m + 1) : 0;
  }
  // CellRank fate_confidence = max absorption probability per cell.
  // Fate confidence is already set by compute_absorption_probabilities:
  // 1.0 for terminal cells, best_p/total for non-terminal cells.

  return List::create(
    _["transition_matrix"] = T,
    _["stationary_distribution"] = pi,
    _["eigenvalues"] = eigenvalues,
    _["schur_vectors"] = schur_vecs,
    _["chi"] = wrap(chi),
    _["coarse_transition"] = wrap(P_coarse),
    _["macrostate_assignment"] = macro,
    _["n_macrostates"] = M,
    _["terminal_states"] = terminal_obs,
    _["lineage_assignment"] = term_state,
    _["fate_confidence"] = fate_conf,
    _["absorption_probabilities"] = abs_prob,
    _["n_terminal_states"] = nT,
    _["method"] = "GPCCA"
  );
}

// ── 10. Lineage drivers ──────────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_lineage_drivers_cpp(
    NumericMatrix expression,
    NumericMatrix abs_probs,
    IntegerVector lineage_idx = IntegerVector())
{
  const int n_genes = expression.nrow();
  const int n_cells = expression.ncol();
  const int n_lineages = abs_probs.ncol();

  if (lineage_idx.size() == 0) {
    lineage_idx = IntegerVector(n_lineages);
    for (int i = 0; i < n_lineages; ++i) lineage_idx[i] = i + 1;
  }

  int n_out = lineage_idx.size();
  NumericMatrix corr(n_genes, n_out);
  NumericMatrix pval(n_genes, n_out);
  NumericVector means(n_genes);
  NumericVector vars(n_genes);

  for (int g = 0; g < n_genes; ++g) {
    double m = 0.0;
    for (int c = 0; c < n_cells; ++c) m += expression(g, c);
    m /= n_cells;
    means[g] = m;
    double v = 0.0;
    for (int c = 0; c < n_cells; ++c) v += (expression(g, c) - m) * (expression(g, c) - m);
    v /= (n_cells - 1);
    vars[g] = v;
  }

  for (int li = 0; li < n_out; ++li) {
    int l = lineage_idx[li] - 1;
    if (l < 0 || l >= n_lineages) continue;

    double mean_prob = 0.0;
    for (int c = 0; c < n_cells; ++c) mean_prob += abs_probs(c, l);
    mean_prob /= n_cells;
    double var_prob = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      double d = abs_probs(c, l) - mean_prob;
      var_prob += d * d;
    }
    var_prob /= (n_cells - 1);

    for (int g = 0; g < n_genes; ++g) {
      double cov = 0.0;
      for (int c = 0; c < n_cells; ++c)
        cov += (expression(g, c) - means[g]) * (abs_probs(c, l) - mean_prob);
      cov /= (n_cells - 1);

      double denom = std::sqrt(vars[g] * var_prob);
      double r = denom > 1e-12 ? cov / denom : 0.0;
      if (r > 1.0) r = 1.0;
      if (r < -1.0) r = -1.0;
      corr(g, li) = r;

      double t_stat = n_cells > 3 ? r * std::sqrt((n_cells - 2.0) / (1.0 - r * r + 1e-15)) : 0.0;
      pval(g, li) = t_stat > 0 ? 1.0 / (1.0 + t_stat * t_stat) : 1.0;
    }
  }

  return List::create(
    _["correlation"] = corr,
    _["pval"] = pval,
    _["lineage_idx"] = lineage_idx
  );
}
