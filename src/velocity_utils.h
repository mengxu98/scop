#ifndef SCOP_VELOCITY_UTILS_H
#define SCOP_VELOCITY_UTILS_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <limits>

using namespace Rcpp;

namespace scop_util {

// ═══════════════════════════════════════════════════════════════════════════
// 1. POWER ITERATION — top-k eigenvalues/vectors of symmetric matrix
// ═══════════════════════════════════════════════════════════════════════════

inline void power_iteration_topk(
    const std::vector<double>& mat,     // column-major n×n
    int n, int k, int max_iter,
    std::vector<double>& eigvals,
    std::vector<double>& eigvecs)
{
    eigvals.assign(k, 0.0);
    eigvecs.assign(n * k, 0.0);
    if (n <= 0 || k <= 0) return;

    for (int comp = 0; comp < k && comp < n; ++comp) {
        int off = comp * n;
        for (int i = 0; i < n; ++i)
            eigvecs[off + i] = (i * 1103515245L + 12345L + comp * 777L) % 10000 / 10000.0;
        double norm = 0.0;
        for (int i = 0; i < n; ++i) norm += eigvecs[off + i] * eigvecs[off + i];
        norm = std::sqrt(norm);
        if (norm > 0) for (int i = 0; i < n; ++i) eigvecs[off + i] /= norm;

        double lambda = 0.0;
        std::vector<double> Av(n);
        for (int iter = 0; iter < max_iter; ++iter) {
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int j = 0; j < n; ++j)
                    sum += mat[j * n + i] * eigvecs[off + j];
                Av[i] = sum;
            }
            for (int p = 0; p < comp; ++p) {
                int poff = p * n;
                double dot = 0.0;
                for (int i = 0; i < n; ++i) dot += eigvecs[poff + i] * Av[i];
                for (int i = 0; i < n; ++i) Av[i] -= dot * eigvecs[poff + i];
            }
            norm = 0.0;
            for (int i = 0; i < n; ++i) norm += Av[i] * Av[i];
            norm = std::sqrt(norm);
            if (norm < 1e-15) break;
            for (int i = 0; i < n; ++i) eigvecs[off + i] = Av[i] / norm;

            double num = 0.0, den = 0.0;
            for (int i = 0; i < n; ++i) {
                double s = 0.0;
                for (int j = 0; j < n; ++j)
                    s += mat[j * n + i] * eigvecs[off + j];
                num += eigvecs[off + i] * s;
                den += eigvecs[off + i] * eigvecs[off + i];
            }
            double nl = den > 0 ? num / den : 0.0;
            if (std::abs(nl - lambda) < 1e-10) break;
            lambda = nl;
        }
        eigvals[comp] = lambda;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 2. COSINE PROJECTION — gene-space velocity → embedding-space projection
// ═══════════════════════════════════════════════════════════════════════════

inline void cosine_projection_embedding(
    const NumericMatrix& gene_velocity,  // genes × cells
    const NumericMatrix& Ms,             // genes × cells
    const IntegerMatrix& knn_idx,        // cells × k (1-based)
    const NumericMatrix& embedding,      // cells × dims
    NumericMatrix& velo_embedding,       // cells × dims (output)
    NumericVector& confidence,           // cells (output)
    NumericVector& velo_length)          // cells (output)
{
    int n_genes = Ms.nrow();
    int n_cells = Ms.ncol();
    int n_neighbors = knn_idx.ncol();
    int n_dims = embedding.ncol();

    // Velocity norm per cell
    for (int c = 0; c < n_cells; ++c) {
        double sq = 0.0;
        for (int g = 0; g < n_genes; ++g) {
            double v = gene_velocity(g, c);
            sq += v * v;
        }
        velo_length[c] = std::sqrt(std::max(sq, 0.0));
    }

    // Cosine projection onto neighbors in gene space (Ms)
    for (int cell = 0; cell < n_cells; ++cell) {
        double vn = velo_length[cell];
        if (vn <= 0.0) continue;

        double wsum = 0.0;
        int used = 0;
        for (int col = 0; col < n_neighbors; ++col) {
            int nb = knn_idx(cell, col);
            if (nb == NA_INTEGER) continue;
            nb -= 1;
            if (nb < 0 || nb >= n_cells || nb == cell) continue;

            double dot = 0.0, ndsq = 0.0;
            for (int g = 0; g < n_genes; ++g) {
                double d = Ms(g, nb) - Ms(g, cell);
                dot += gene_velocity(g, cell) * d;
                ndsq += d * d;
            }
            double nd = std::sqrt(ndsq);
            if (nd <= 0.0) continue;
            double cosval = dot / (vn * nd);
            if (cosval <= 0.0 || !std::isfinite(cosval)) continue;

            for (int d = 0; d < n_dims; ++d)
                velo_embedding(cell, d) += cosval * (embedding(nb, d) - embedding(cell, d));
            wsum += cosval;
            ++used;
        }
        if (wsum > 0.0) {
            for (int d = 0; d < n_dims; ++d)
                velo_embedding(cell, d) /= wsum;
            confidence[cell] = (double)used / (double)n_neighbors;
        }
    }
}

inline void velocity_confidence_py(
    const NumericMatrix& gene_velocity,  // genes × cells
    const IntegerMatrix& knn_idx,        // cells × k (1-based)
    NumericVector& confidence,           // cells (output)
    NumericVector& velocity_length)      // cells (output)
{
    int n_genes = gene_velocity.nrow();
    int n_cells = gene_velocity.ncol();
    int n_neighbors = knn_idx.ncol();

    std::vector<double> centered(n_genes * n_cells, 0.0);
    std::vector<double> norms(n_cells, 0.0);

    for (int cell = 0; cell < n_cells; ++cell) {
        double mean = 0.0;
        int n_valid = 0;
        for (int g = 0; g < n_genes; ++g) {
            double v = gene_velocity(g, cell);
            if (!std::isfinite(v)) continue;
            mean += v;
            ++n_valid;
        }
        if (n_valid > 0) mean /= static_cast<double>(n_valid);

        double ss = 0.0;
        for (int g = 0; g < n_genes; ++g) {
            double v = gene_velocity(g, cell);
            double cv = std::isfinite(v) ? v - mean : 0.0;
            centered[cell * n_genes + g] = cv;
            ss += cv * cv;
        }
        norms[cell] = std::sqrt(std::max(0.0, ss));
        velocity_length[cell] = std::round(norms[cell] * 100.0) / 100.0;
    }

    for (int cell = 0; cell < n_cells; ++cell) {
        if (norms[cell] <= 0.0) {
            confidence[cell] = 0.0;
            continue;
        }
        double r_sum = 0.0;
        int used = 0;
        for (int col = 0; col < n_neighbors; ++col) {
            int nb = knn_idx(cell, col);
            if (nb == NA_INTEGER) continue;
            nb -= 1;
            if (nb < 0 || nb >= n_cells || nb == cell || norms[nb] <= 0.0) continue;

            double dot = 0.0;
            for (int g = 0; g < n_genes; ++g)
                dot += centered[nb * n_genes + g] * centered[cell * n_genes + g];
            double corr = dot / (norms[nb] * norms[cell]);
            if (std::isfinite(corr)) {
                r_sum += corr;
                ++used;
            }
        }
        double r = used > 0 ? r_sum / static_cast<double>(used) : 0.0;
        confidence[cell] = std::max(0.0, r);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 3. OLS GAMMA — spliced/unspliced ratio through origin
// ═══════════════════════════════════════════════════════════════════════════

inline void ols_gamma_origin(
    const NumericMatrix& Ms,
    const NumericMatrix& Mu,
    NumericVector& gamma)
{
    int n_genes = Ms.nrow();
    int n_cells = Ms.ncol();
    for (int g = 0; g < n_genes; ++g) {
        double num = 0.0, den = 0.0;
        for (int c = 0; c < n_cells; ++c) {
            double s = Ms(g, c);
            num += s * Mu(g, c);
            den += s * s;
        }
        gamma[g] = den > 1e-12 ? num / den : 0.0;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 4. DSU (UNION-FIND) + EDGE — for Kruskal MST
// ═══════════════════════════════════════════════════════════════════════════

struct DSU {
    std::vector<int> parent, rank;
    explicit DSU(int n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }
    bool unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        if (rank[ra] < rank[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) ++rank[ra];
        return true;
    }
};

struct Edge {
    int i; int j; double weight;
};

// Build MST from edge list (sorted descending by weight)
inline NumericMatrix build_mst_matrix(
    const std::vector<Edge>& edges,
    int n_nodes)
{
    NumericMatrix tree(n_nodes, n_nodes);
    DSU dsu(n_nodes);
    int used = 0;
    for (const auto& e : edges) {
        if (dsu.unite(e.i, e.j)) {
            tree(e.i, e.j) = e.weight;
            ++used;
            if (used == n_nodes - 1) break;
        }
    }
    return tree;
}

// ═══════════════════════════════════════════════════════════════════════════
// 5. GAUSSIAN ELIMINATION — solve (I - Q) * X = R
// ═══════════════════════════════════════════════════════════════════════════

inline bool gaussian_elimination_solve(
    const NumericMatrix& Q, int nQ,
    const NumericMatrix& R, int nT,
    NumericMatrix& X)
{
    // Augmented matrix: [Q | R] of size nQ × (nQ + nT)
    std::vector<std::vector<double>> aug(nQ, std::vector<double>(nQ + nT));
    for (int i = 0; i < nQ; ++i) {
        for (int j = 0; j < nQ; ++j) aug[i][j] = Q(i, j);
        for (int j = 0; j < nT; ++j) aug[i][nQ + j] = R(i, j);
    }

    bool success = true;
    for (int col = 0; col < nQ; ++col) {
        int pivot = -1;
        double max_val = 0.0;
        for (int row = col; row < nQ; ++row) {
            if (std::abs(aug[row][col]) > max_val) {
                max_val = std::abs(aug[row][col]);
                pivot = row;
            }
        }
        if (pivot < 0 || max_val < 1e-14) { success = false; break; }
        std::swap(aug[col], aug[pivot]);
        double scale = aug[col][col];
        for (int j = 0; j < nQ + nT; ++j) aug[col][j] /= scale;
        for (int row = 0; row < nQ; ++row) {
            if (row == col) continue;
            double factor = aug[row][col];
            for (int j = 0; j < nQ + nT; ++j) aug[row][j] -= factor * aug[col][j];
        }
    }

    if (success) {
        X = NumericMatrix(nQ, nT);
        for (int i = 0; i < nQ; ++i)
            for (int j = 0; j < nT; ++j)
                X(i, j) = aug[i][nQ + j];
    }
    return success;
}

// ═══════════════════════════════════════════════════════════════════════════
// 6. BUILD VELOCITY TRANSITION — embedding-space cosine transition matrix
// ═══════════════════════════════════════════════════════════════════════════

inline void build_velocity_transition(
    const NumericMatrix& velocity_embedding,   // cells × dims
    const NumericMatrix& embedding,            // cells × dims
    const IntegerMatrix& knn_idx,              // cells × k (1-based)
    int n_neighbors_velo,
    std::vector<double>& T)                    // n × n column-major (output)
{
    int n_cells = embedding.nrow();
    int n_dims = embedding.ncol();
    int n_k = knn_idx.ncol();
    int nk = std::min(n_neighbors_velo, n_k);

    T.assign(n_cells * n_cells, 0.0);

    for (int i = 0; i < n_cells; ++i) {
        double vn = 0.0;
        for (int d = 0; d < n_dims; ++d)
            vn += velocity_embedding(i, d) * velocity_embedding(i, d);
        vn = std::sqrt(vn);
        if (vn < 1e-10) continue;

        double row_sum = 0.0;
        std::vector<double> row(n_cells, 0.0);

        for (int col = 0; col < nk; ++col) {
            int nb = knn_idx(i, col);
            if (nb == NA_INTEGER) continue;
            nb -= 1;
            if (nb < 0 || nb >= n_cells || nb == i) continue;
            double dot = 0.0, dn = 0.0;
            for (int d = 0; d < n_dims; ++d) {
                double delta = embedding(nb, d) - embedding(i, d);
                dot += velocity_embedding(i, d) * delta;
                dn += delta * delta;
            }
            dn = std::sqrt(dn);
            if (dn < 1e-10) continue;
            double cosine = dot / (vn * dn);
            if (cosine > 0) { row[nb] = cosine; row_sum += cosine; }
        }

        // Backward transitions
        for (int col = 0; col < n_k; ++col) {
            int nb = knn_idx(i, col);
            if (nb == NA_INTEGER) continue;
            nb -= 1;
            if (nb < 0 || nb >= n_cells || nb == i) continue;
            double dot = 0.0, dn = 0.0;
            for (int d = 0; d < n_dims; ++d) {
                double delta = embedding(i, d) - embedding(nb, d);
                dot -= velocity_embedding(i, d) * delta;
                dn += delta * delta;
            }
            dn = std::sqrt(dn);
            if (dn < 1e-10) continue;
            double cosine = dot / (vn * dn);
            if (cosine > 0 && row[nb] == 0.0) {
                row[nb] = cosine * 0.5;
                row_sum += cosine * 0.5;
            }
        }

        if (row_sum > 0) {
            for (int j = 0; j < n_cells; ++j)
                T[i + j * n_cells] = row[j] / row_sum;
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 7. STATIONARY DISTRIBUTION — power iteration on stochasic matrix
// ═══════════════════════════════════════════════════════════════════════════

inline NumericVector stationary_distribution(
    const NumericMatrix& T_, int max_iter = 1000, double tol = 1e-10)
{
    int n = T_.nrow();
    NumericVector pi(n, 1.0 / n);

    for (int iter = 0; iter < max_iter; ++iter) {
        NumericVector next_pi(n);
        // Row-stochastic matrix: pi * T = pi
        // pi[j] = sum_i pi[i] * T[i][j]
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int i = 0; i < n; ++i)
                s += pi[i] * T_(i, j);
            next_pi[j] = s;
        }
        double sum = 0.0;
        for (int i = 0; i < n; ++i) sum += next_pi[i];
        if (sum < 1e-15) { next_pi = NumericVector(n, 1.0/n); sum = 1.0; }
        for (int i = 0; i < n; ++i) next_pi[i] /= sum;
        double delta = 0.0;
        for (int i = 0; i < n; ++i) delta = std::max(delta, std::abs(next_pi[i] - pi[i]));
        pi = next_pi;
        if (delta < tol) break;
    }
    return pi;
}

// ═══════════════════════════════════════════════════════════════════════════
// 8. TRANSITION MATRIX VALIDATION
// ═══════════════════════════════════════════════════════════════════════════

inline void validate_transition_matrix(
    NumericMatrix& T,
    double eps = 1e-10,
    double min_self_loop = 0.01,
    int* p_nans = nullptr,
    int* p_negs = nullptr,
    int* p_zeros = nullptr,
    int* p_loops = nullptr)
{
    int n = T.nrow();
    int nans = 0, negs = 0, zeros = 0, loops = 0;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (!std::isfinite(T(i, j))) { T(i, j) = 0.0; ++nans; }

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (T(i, j) < -1e-10) { T(i, j) = 0.0; ++negs; }
            else if (T(i, j) < 0) T(i, j) = 0.0;

    for (int i = 0; i < n; ++i)
        if (T(i, i) < min_self_loop) { T(i, i) = min_self_loop; ++loops; }

    for (int i = 0; i < n; ++i) {
        double rs = 0.0;
        for (int j = 0; j < n; ++j) rs += T(i, j);
        if (rs < eps) {
            for (int j = 0; j < n; ++j) T(i, j) = 0.0;
            T(i, i) = 1.0;
            ++zeros;
        }
    }
    for (int i = 0; i < n; ++i) {
        double rs = 0.0;
        for (int j = 0; j < n; ++j) rs += T(i, j);
        if (rs > eps) for (int j = 0; j < n; ++j) T(i, j) /= rs;
    }

    if (p_nans) *p_nans = nans;
    if (p_negs) *p_negs = negs;
    if (p_zeros) *p_zeros = zeros;
    if (p_loops) *p_loops = loops;
}

// ═══════════════════════════════════════════════════════════════════════════
// 9. TERMINAL STATE DETECTION — from coarse transition matrix
// ═══════════════════════════════════════════════════════════════════════════

inline void detect_terminal_states(
    const NumericMatrix& P_coarse,
    int M,
    std::vector<int>& is_terminal,
    std::vector<int>& term_idx,
    std::vector<int>& trans_idx,
    int& n_terminal)
{
    // Sort by self-transition probability (descending)
    std::vector<std::pair<double, int>> self_trans;
    for (int a = 0; a < M; ++a)
        self_trans.push_back({P_coarse(a, a), a});
    std::sort(self_trans.begin(), self_trans.end(),
        std::greater<std::pair<double, int>>());

    // Find the largest relative gap in self-transition probabilities
    n_terminal = 1;
    for (int a = 1; a < M; ++a) {
        double gap = self_trans[a - 1].first - self_trans[a].first;
        if (gap > 0.10) { n_terminal = a; break; }
        n_terminal = a + 1;
    }
    if (n_terminal >= M) n_terminal = M - 1;
    // Ensure sensible range
    if (n_terminal < 1) n_terminal = 1;
    if (M >= 3 && n_terminal < 2) n_terminal = 2;

    is_terminal.assign(M, 0);
    for (int a = 0; a < n_terminal; ++a)
        is_terminal[self_trans[a].second] = 1;

    term_idx.clear();
    trans_idx.clear();
    for (int a = 0; a < M; ++a) {
        if (is_terminal[a]) term_idx.push_back(a);
        else trans_idx.push_back(a);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 10. ABSORPTION PROBABILITY — from macrostate + coarse matrix
// ═══════════════════════════════════════════════════════════════════════════

inline void compute_absorption_probabilities(
    const NumericMatrix& P_coarse,
    const IntegerVector& macro,
    int n, int M,
    const std::vector<int>& is_terminal,
    const std::vector<int>& term_idx,
    const std::vector<int>& trans_idx,
    int nT, int nQ,
    NumericMatrix& abs_prob,
    IntegerVector& term_state,
    NumericVector& fate_conf)
{
    NumericMatrix B_r;
    bool solve_ok = false;

    if (nQ > 0 && nT > 0) {
        NumericMatrix Q_r(nQ, nQ), R_r(nQ, nT);
        for (int i = 0; i < nQ; ++i) {
            for (int j = 0; j < nQ; ++j)
                Q_r(i, j) = (i == j ? 1.0 : 0.0) - P_coarse(trans_idx[i], trans_idx[j]);
            for (int j = 0; j < nT; ++j)
                R_r(i, j) = P_coarse(trans_idx[i], term_idx[j]);
        }
        solve_ok = gaussian_elimination_solve(Q_r, nQ, R_r, nT, B_r);
        if (!solve_ok) {
            // Diagonal approximation fallback
            B_r = NumericMatrix(nQ, nT);
            for (int i = 0; i < nQ; ++i) {
                double d = Q_r(i, i);
                if (std::abs(d) > 1e-10)
                    for (int j = 0; j < nT; ++j) B_r(i, j) = R_r(i, j) / d;
                else
                    for (int j = 0; j < nT; ++j) B_r(i, j) = 1.0 / nT;
            }
            solve_ok = true;
        }
    }

    for (int i = 0; i < n; ++i) {
        int m = macro[i] - 1;
        if (is_terminal[m]) {
            for (int t = 0; t < nT; ++t)
                abs_prob(i, t) = (term_idx[t] == m) ? 1.0 : 0.0;
            term_state[i] = m + 1;
            fate_conf[i] = 1.0;
        } else if (solve_ok && nQ > 0 && nT > 0) {
            int qi = -1;
            for (int q = 0; q < nQ; ++q) if (trans_idx[q] == m) { qi = q; break; }
            double total = 0.0;
            int best_t = 0;
            double best_p = 0.0;
            for (int t = 0; t < nT; ++t) {
                double p = (qi >= 0) ? B_r(qi, t) : 1.0 / nT;
                if (p < 0) p = 0.0;
                abs_prob(i, t) = p;
                total += p;
                if (p > best_p) { best_p = p; best_t = t; }
            }
            if (total > 0) {
                for (int t = 0; t < nT; ++t) abs_prob(i, t) /= total;
            }
            term_state[i] = term_idx[best_t] + 1;
            fate_conf[i] = best_p / (total > 0 ? total : 1.0);
        } else if (nT > 0) {
            for (int t = 0; t < nT; ++t)
                abs_prob(i, t) = (term_idx[t] == m) ? 1.0 : 0.0;
            term_state[i] = m + 1;
            fate_conf[i] = 1.0;
        }
    }
}

} // namespace scop_util

#endif // SCOP_VELOCITY_UTILS_H
