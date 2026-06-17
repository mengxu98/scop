#include <Rcpp.h>
#include "velocity_utils.h"

using namespace Rcpp;

// ── Core ODE solution (closed-form) ──────────────────────────────────────────

// du/dt = alpha - beta * u   →   u(tau) = u0*exp(-beta*tau) + alpha/beta * (1 - exp(-beta*tau))
inline double unspliced(double tau, double u0, double alpha, double beta) {
  double e = std::exp(-beta * tau);
  // avoid inf when beta=0
  if (!std::isfinite(e) || std::abs(beta) < 1e-12) {
    return u0 + alpha * tau;
  }
  double steady = alpha / beta;
  return u0 * e + steady * (1.0 - e);
}

// ds/dt = beta*u - gamma*s   →   closed form with induction from (u0,s0)
inline double spliced(double tau, double s0, double u0,
                       double alpha, double beta, double gamma) {
  double eu = std::exp(-beta * tau);
  double es = std::exp(-gamma * tau);
  if (!std::isfinite(eu) || !std::isfinite(es))
    return s0; // fallback

  double alpha_g = alpha / gamma;
  double diff = gamma - beta;
  double c = (std::abs(diff) > 1e-10)
    ? (alpha - u0 * beta) / diff
    : 0.0;

  return s0 * es + alpha_g * (1.0 - es) + c * (es - eu);
}

// ── Invert ODE to get tau from (u,s) ─────────────────────────────────────────

inline double tau_inv_u(double u, double uinf, double u0, double beta) {
  if (std::abs(beta) < 1e-12 || std::abs(u - uinf) < 1e-12)
    return 0.0;
  double arg = std::max((u - uinf) / (u0 - uinf), 1e-12);
  return -std::log(arg) / beta;
}

// ── Nelder-Mead Simplex Optimizer ─────────────────────────────────────────────

struct NMState {
  int n;                       // dimension
  std::vector<double> x0;      // initial guess
  std::vector<double> lb;      // lower bounds
  std::vector<double> ub;      // upper bounds
  int max_iter;
  double tol;
  double alpha_r;              // reflection coefficient (1.0)
  double gamma_e;              // expansion coefficient (2.0)
  double rho_c;                // contraction coefficient (0.5)
  double sigma_s;              // shrink coefficient (0.5)
};

// Evaluate the objective function at point x
typedef std::function<double(const std::vector<double>&)> NMFun;

static std::vector<double> nm_unpack_simplex(
    const std::vector<double>& simplex, int n, int idx)
{
  return std::vector<double>(simplex.begin() + idx * n, simplex.begin() + (idx + 1) * n);
}

static double nm_norm(const std::vector<double>& v) {
  double s = 0.0;
  for (double vi : v) s += vi * vi;
  return std::sqrt(s);
}

std::vector<double> nelder_mead(
    NMState& state, NMFun fun,
    double* final_loss = nullptr)
{
  int n = state.n;
  int npts = n + 1;
  std::vector<double> simplex(npts * n); // simplex points, row-major
  std::vector<double> fvals(npts);

  // Initialize simplex: x0, x0 + ei for each i
  for (int j = 0; j < n; ++j)
    simplex[0 * n + j] = state.x0[j];
  fvals[0] = fun(state.x0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      simplex[(i + 1) * n + j] = state.x0[j];
    double delta = std::max(std::abs(state.x0[i]) * 0.05, 0.00025);
    simplex[(i + 1) * n + i] += delta;
    // Project to bounds
    for (int j = 0; j < n; ++j) {
      simplex[(i + 1) * n + j] = std::max(state.lb[j],
        std::min(state.ub[j], simplex[(i + 1) * n + j]));
    }
    fvals[i + 1] = fun(nm_unpack_simplex(simplex, n, i + 1));
  }

  for (int iter = 0; iter < state.max_iter; ++iter) {
    // Sort simplex by fval (ascending)
    std::vector<int> idx(npts);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
      return fvals[a] < fvals[b];
    });

    // Check convergence
    double f_std = 0.0, f_mean = 0.0;
    for (int i = 0; i < npts; ++i) f_mean += fvals[i];
    f_mean /= npts;
    for (int i = 0; i < npts; ++i) f_std += (fvals[i] - f_mean) * (fvals[i] - f_mean);
    f_std = std::sqrt(f_std / npts);
    if (f_std < state.tol) break;

    int best = idx[0], worst = idx[n], s_worst = idx[n - 1];
    auto xb = nm_unpack_simplex(simplex, n, best);
    auto xw = nm_unpack_simplex(simplex, n, worst);

    // Centroid of all points except worst
    std::vector<double> x0(n, 0.0);
    for (int i = 0; i < npts; ++i) {
      if (i == worst) continue;
      auto xi = nm_unpack_simplex(simplex, n, i);
      for (int j = 0; j < n; ++j) x0[j] += xi[j];
    }
    for (int j = 0; j < n; ++j) x0[j] /= n;

    // Reflection
    std::vector<double> xr(n);
    for (int j = 0; j < n; ++j) {
      xr[j] = x0[j] + state.alpha_r * (x0[j] - xw[j]);
      xr[j] = std::max(state.lb[j], std::min(state.ub[j], xr[j]));
    }
    double fr = fun(xr);

    if (fr < fvals[best]) {
      // Expansion
      std::vector<double> xe(n);
      for (int j = 0; j < n; ++j) {
        xe[j] = x0[j] + state.gamma_e * (xr[j] - x0[j]);
        xe[j] = std::max(state.lb[j], std::min(state.ub[j], xe[j]));
      }
      double fe = fun(xe);
      if (fe < fr) {
        for (int j = 0; j < n; ++j) simplex[worst * n + j] = xe[j];
        fvals[worst] = fe;
      } else {
        for (int j = 0; j < n; ++j) simplex[worst * n + j] = xr[j];
        fvals[worst] = fr;
      }
    } else if (fr < fvals[s_worst]) {
      for (int j = 0; j < n; ++j) simplex[worst * n + j] = xr[j];
      fvals[worst] = fr;
    } else {
      // Contraction
      std::vector<double> xc(n);
      if (fr < fvals[worst]) {
        // Outside contraction
        for (int j = 0; j < n; ++j) {
          xc[j] = x0[j] + state.rho_c * (xr[j] - x0[j]);
          xc[j] = std::max(state.lb[j], std::min(state.ub[j], xc[j]));
        }
      } else {
        // Inside contraction
        for (int j = 0; j < n; ++j) {
          xc[j] = x0[j] + state.rho_c * (xw[j] - x0[j]);
          xc[j] = std::max(state.lb[j], std::min(state.ub[j], xc[j]));
        }
      }
      double fc = fun(xc);
      if (fc < std::min(fr, fvals[worst])) {
        for (int j = 0; j < n; ++j) simplex[worst * n + j] = xc[j];
        fvals[worst] = fc;
      } else {
        // Shrink toward best point
        for (int i = 0; i < npts; ++i) {
          if (i == best) continue;
          for (int j = 0; j < n; ++j) {
            simplex[i * n + j] = xb[j] + state.sigma_s * (simplex[i * n + j] - xb[j]);
            simplex[i * n + j] = std::max(state.lb[j],
              std::min(state.ub[j], simplex[i * n + j]));
          }
          fvals[i] = fun(nm_unpack_simplex(simplex, n, i));
        }
      }
    }
  }

  // Return best point
  std::vector<int> idx2(npts);
  std::iota(idx2.begin(), idx2.end(), 0);
  std::sort(idx2.begin(), idx2.end(), [&](int a, int b) {
    return fvals[a] < fvals[b];
  });
  if (final_loss) *final_loss = fvals[idx2[0]];
  return nm_unpack_simplex(simplex, n, idx2[0]);
}


// ── Dynamical velocity fitting (per gene, single EM iteration) ───────────────

struct DynamicalPars {
  double alpha;   // transcription rate (post-induction)
  double beta;    // splicing rate
  double gamma;   // degradation rate
  double t_;      // switching time
  double scaling; // u/s scaling factor
  double loss;    // final MSE
};

static DynamicalPars fit_one_gene_dynamical(
    const double* u, const double* s,  // cell-level data (n_cells)
    int n_cells,
    double alpha_init, double beta_init,
    double gamma_init, double t_init,
    int max_iter_nm)
{
  // Determine active cells (u>0 and s>0)
  std::vector<int> active;
  for (int i = 0; i < n_cells; ++i)
    if (u[i] > 0 && s[i] > 0)
      active.push_back(i);
  int na = (int)active.size();
  if (na < 10) {  // need enough cells for fitting
    return {alpha_init, beta_init, gamma_init, 0.0, 1.0, 1e10};
  }

  // Extract active data
  std::vector<double> ua(na), sa(na);
  for (int i = 0; i < na; ++i) {
    ua[i] = u[active[i]];
    sa[i] = s[active[i]];
  }

  // Simple initialization from data
  double u_max = *std::max_element(ua.begin(), ua.end());
  double s_max = *std::max_element(sa.begin(), sa.end());
  double u_min_pos = u_max;
  for (double v : ua) if (v > 0 && v < u_min_pos) u_min_pos = v;
  if (u_min_pos >= u_max) u_min_pos = u_max * 0.1;

  // Estimate beta from decay: beta ~ -log(u_min/u_max) / tau_range
  // gamma from s/u ratio at steady state
  double beta_est = beta_init;
  double gamma_est = gamma_init;
  double alpha_est = u_max * beta_est;
  double t_est = t_init;

  if (beta_init <= 0) beta_est = 1.0;
  if (gamma_init <= 0) gamma_est = 0.5;
  if (t_init <= 0) t_est = 2.0;

  // Parameter bounds
  double alpha_lb = 0.01, alpha_ub = 100.0;
  double beta_lb  = 0.01, beta_ub  = 50.0;
  double gamma_lb = 0.01, gamma_ub = 50.0;
  double t_lb     = 0.01, t_ub    = 20.0;

  // Clip initial values to bounds
  alpha_est = std::max(alpha_lb, std::min(alpha_ub, alpha_est));
  beta_est  = std::max(beta_lb,  std::min(beta_ub,  beta_est));
  gamma_est = std::max(gamma_lb, std::min(gamma_ub, gamma_est));
  t_est     = std::max(t_lb,     std::min(t_ub,     t_est));

  // Nelder-Mead on (alpha, beta, gamma, t_)
  // Objective: for each cell, invert ODE to get tau, predict (u,s), compute MSE
  NMState state;
  state.n = 4;
  state.x0 = {alpha_est, beta_est, gamma_est, t_est};
  state.lb = {alpha_lb, beta_lb, gamma_lb, t_lb};
  state.ub = {alpha_ub, beta_ub, gamma_ub, t_ub};
  state.max_iter = max_iter_nm;
  state.tol = 1e-4;
  state.alpha_r = 1.0;
  state.gamma_e = 2.0;
  state.rho_c = 0.5;
  state.sigma_s = 0.5;

  auto obj_fun = [&](const std::vector<double>& x) -> double {
    double a = x[0], b = x[1], g = x[2], ts = x[3];
    double u_inf = a / b;
    double mse = 0.0;
    int used = 0;
    for (int i = 0; i < na; ++i) {
      // Invert tau from u (unspliced equation is monotonic)
      double ui = ua[i], si = sa[i];
      // tau from unspliced: u = u0 + (u_inf - u0)*(1 - exp(-b*tau))
      // if ui >= u_inf, cell is beyond steady state
      double tau_i;
      if (ui >= u_inf * 0.99) {
        tau_i = ts * 3.0; // far in the future
      } else {
        double arg = 1.0 - ui / u_inf;
        if (arg < 1e-12) arg = 1e-12;
        tau_i = -std::log(arg) / b;
      }
      // Predict spliced
      double sp = spliced(tau_i, 0.0, 0.0, a, b, g);
      double up = unspliced(tau_i, 0.0, a, b);
      double e = (ui - up) * (ui - up) + (si - sp) * (si - sp);
      if (std::isfinite(e)) {
        mse += e;
        ++used;
      }
    }
    return used > 0 ? mse / used : 1e10;
  };

  double final_loss;
  auto opt = nelder_mead(state, obj_fun, &final_loss);

  return {opt[0], opt[1], opt[2], opt[3], 1.0, final_loss};
}

// ── Rcpp export: fit dynamical model per gene (Nelder-Mead, legacy) ───────

// [[Rcpp::export]]
List scvelo_dynamical_nm_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    IntegerVector use_genes,  // 1-based indices of genes to fit
    int max_iter = 10,
    double init_alpha = -1.0,
    double init_beta  = -1.0,
    double init_gamma = -1.0)
{
  int n_genes = Ms.nrow();
  int n_cells = Ms.ncol();
  if (Mu.nrow() != n_genes || Mu.ncol() != n_cells)
    stop("Ms and Mu must have identical dimensions");

  // Which genes to fit
  std::vector<bool> do_fit(n_genes, true);
  if (use_genes.size() > 0) {
    do_fit.assign(n_genes, false);
    for (int i = 0; i < use_genes.size(); ++i) {
      int g = use_genes[i] - 1;
      if (g >= 0 && g < n_genes) do_fit[g] = true;
    }
  }

  // Output arrays
  NumericVector alpha_out(n_genes);
  NumericVector beta_out(n_genes);
  NumericVector gamma_out(n_genes);
  NumericVector t_out(n_genes);
  NumericVector loss_out(n_genes);
  IntegerVector converged(n_genes, 0);

  int n_fit = 0;
  for (int g = 0; g < n_genes; ++g) {
    if (!do_fit[g]) {
      alpha_out[g] = 0; beta_out[g] = 0; gamma_out[g] = 0;
      t_out[g] = 0; loss_out[g] = 1e10;
      continue;
    }

    // Copy per-gene data
    std::vector<double> u(n_cells), s(n_cells);
    int valid = 0;
    for (int c = 0; c < n_cells; ++c) {
      u[c] = Mu(g, c);
      s[c] = Ms(g, c);
      if (u[c] > 0 && s[c] > 0) ++valid;
    }
    if (valid < 10) {
      alpha_out[g] = 0; beta_out[g] = 0; gamma_out[g] = 0;
      t_out[g] = 0; loss_out[g] = 1e10;
      continue;
    }

    // Initialization: deterministic gamma as starting point
    double gamma_det = 0.0;
    {
      double num = 0.0, den = 0.0;
      for (int c = 0; c < n_cells; ++c) {
        double si = Ms(g, c), ui = Mu(g, c);
        num += si * ui;
        den += si * si;
      }
      gamma_det = den > 1e-12 ? num / den : 0.5;
    }
    double ms_mean = 0.0;
    for (int c = 0; c < n_cells; ++c) ms_mean += Ms(g, c);
    ms_mean /= n_cells;
    double alpha_init = init_alpha > 0
      ? init_alpha
      : gamma_det * (ms_mean > 0.1 ? ms_mean : 0.1);
    double beta_init  = init_beta > 0 ? init_beta : 1.0;
    double gamma_init = init_gamma > 0 ? init_gamma : std::max(0.1, gamma_det);
    double t_init = 2.0;

    int nm_iter = std::max(5, std::min(max_iter / 5, 50));
    DynamicalPars par = fit_one_gene_dynamical(
      u.data(), s.data(), n_cells,
      alpha_init, beta_init, gamma_init, t_init,
      nm_iter
    );

    alpha_out[g] = par.alpha;
    beta_out[g]  = par.beta;
    gamma_out[g] = par.gamma;
    t_out[g]     = par.t_;
    loss_out[g]  = par.loss;
    converged[g] = par.loss < 1e9 ? 1 : 0;
    ++n_fit;

    // Progress callback for R
    if (n_fit % 50 == 0) Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["alpha"] = alpha_out,
    _["beta"] = beta_out,
    _["gamma"] = gamma_out,
    _["t_"] = t_out,
    _["loss"] = loss_out,
    _["converged"] = converged,
    _["n_fitted"] = n_fit
  );
}

// [[Rcpp::export]]
List scvelo_dynamical_velocity_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    NumericVector alpha,
    NumericVector beta,
    NumericVector gamma,
    NumericVector t_,
    IntegerMatrix knn_idx,
    NumericMatrix embedding)
{
  int n_genes = Ms.nrow();
  int n_cells = Ms.ncol();
  int n_neighbors = knn_idx.ncol();
  int n_dims = embedding.ncol();

  if ((int)alpha.size() < n_genes || (int)beta.size() < n_genes ||
      (int)gamma.size() < n_genes || (int)t_.size() < n_genes)
    stop("parameter vectors must match n_genes");

  // Compute per-cell velocity from fitted parameters
  // v = du/dt = alpha - beta*u  (unspliced velocity)
  NumericMatrix velocity(n_genes, n_cells);
  for (int g = 0; g < n_genes; ++g) {
    double a = alpha[g], b = beta[g];
    for (int c = 0; c < n_cells; ++c) {
      velocity(g, c) = a - b * Mu(g, c);
    }
  }

  // Velocity embedding (cosine projection, same as stochastic/deterministic)
  NumericMatrix velocity_embedding(n_cells, n_dims);
  NumericVector confidence(n_cells);
  NumericVector velocity_length(n_cells);
  scop_util::cosine_projection_embedding(velocity, Ms, knn_idx, embedding, velocity_embedding, confidence, velocity_length);
  scop_util::velocity_confidence_py(velocity, knn_idx, confidence, velocity_length);

  return List::create(
    _["velocity"] = velocity,
    _["velocity_embedding"] = velocity_embedding,
    _["confidence"] = confidence,
    _["velocity_length"] = velocity_length
  );
}

// ── EM-based dynamical model fitting (per gene) ──────────────────────────────
// Matches scvelo's Expectation-Maximization approach for the dynamical model.
// The model: du/dt = alpha - beta*u (on phase), ds/dt = beta*u - gamma*s
// With switching time t_, cells before t_ follow (alpha, beta); after t_, alpha=0.

// [[Rcpp::export]]
List scvelo_dynamical_em_cpp(
    NumericMatrix Ms,
    NumericMatrix Mu,
    IntegerVector use_genes,
    int max_iter_em = 10,
    double conv_tol = 1e-6,
    int em_oversampling = 2,
    double init_alpha = -1.0,
    double init_beta  = -1.0,
    double init_gamma = -1.0)
{
  const int n_genes = Ms.nrow();
  const int n_cells = Ms.ncol();
  if (Mu.nrow() != n_genes || Mu.ncol() != n_cells)
    stop("Ms and Mu must have identical dimensions");

  std::vector<bool> do_fit(n_genes, true);
  if (use_genes.size() > 0) {
    do_fit.assign(n_genes, false);
    for (int i = 0; i < use_genes.size(); ++i) {
      int g = use_genes[i] - 1;
      if (g >= 0 && g < n_genes) do_fit[g] = true;
    }
  }

  NumericVector alpha_out(n_genes);
  NumericVector beta_out(n_genes);
  NumericVector gamma_out(n_genes);
  NumericVector t_out(n_genes);
  NumericVector scaling_out(n_genes);
  NumericVector loss_out(n_genes);
  IntegerVector converged(n_genes, 0);

  // Shared/u scaling: median of (u_max / s_max) across genes
  double global_scale = 1.0;
  {
    std::vector<double> ratios;
    for (int g = 0; g < n_genes; ++g) {
      if (!do_fit[g]) continue;
      double smax = 0.0, umax = 0.0;
      for (int c = 0; c < n_cells; ++c) {
        if (Ms(g, c) > smax) smax = Ms(g, c);
        if (Mu(g, c) > umax) umax = Mu(g, c);
      }
      if (smax > 0 && umax > 0) ratios.push_back(umax / smax);
    }
    if (!ratios.empty()) {
      std::sort(ratios.begin(), ratios.end());
      global_scale = ratios[ratios.size() / 2];
    }
  }

  int n_fit = 0;
  for (int g = 0; g < n_genes; ++g) {
    if (!do_fit[g]) {
      alpha_out[g] = 0; beta_out[g] = 0; gamma_out[g] = 0;
      t_out[g] = 0; scaling_out[g] = 1.0; loss_out[g] = 1e10;
      continue;
    }

    std::vector<double> u(n_cells), s(n_cells);
    int valid = 0;
    for (int c = 0; c < n_cells; ++c) {
      u[c] = Mu(g, c);
      s[c] = Ms(g, c);
      if (u[c] > 0 && s[c] > 0) ++valid;
    }
    if (valid < 10) {
      alpha_out[g] = 0; beta_out[g] = 0; gamma_out[g] = 0;
      t_out[g] = 0; scaling_out[g] = 1.0; loss_out[g] = 1e10;
      continue;
    }

    // Compute per-gene statistics for initialization
    double s_mean = 0, u_mean = 0;
    double ss_tot = 0, su_tot = 0, uu_tot = 0;
    for (int c = 0; c < n_cells; ++c) {
      s_mean += s[c]; u_mean += u[c];
      ss_tot += s[c] * s[c];
      su_tot += s[c] * u[c];
      uu_tot += u[c] * u[c];
    }
    s_mean /= n_cells; u_mean /= n_cells;

    // OLS gamma = <s,u>/<s,s>
    double gamma_det = (ss_tot > 1e-12) ? su_tot / ss_tot : 0.5;
    double alpha_init = init_alpha > 0 ? init_alpha : gamma_det * (s_mean > 0.1 ? s_mean : 0.1);
    double beta_init  = init_beta > 0 ? init_beta : 1.0;
    double gamma_init = init_gamma > 0 ? init_gamma : std::max(0.1, gamma_det);

    // Clip to bounds
    double a = alpha_init, b = beta_init, gm = gamma_init;
    a = std::max(0.01, std::min(100.0, a));
    b = std::max(0.1, std::min(50.0, b));
    gm = std::max(0.01, std::min(50.0, gm));
    double ts = 2.0;  // initial switching time
    double scaling = global_scale > 0 ? global_scale : 1.0;

    // EM iterations: E-step assigns latent time, M-step optimizes parameters
    double prev_loss = 1e10;
    for (int em_iter = 0; em_iter < max_iter_em; ++em_iter) {
      // E-step: assign each cell a latent time tau based on current params
      // Two regimes:
      //   Phase 1 (induction, t < t_): u(t) = alpha/beta * (1 - exp(-beta*t))
      //   Phase 2 (repression, t >= t_): u(t) = u_ * exp(-beta*(t - t_))
      //                                     where u_ = alpha/beta * (1 - exp(-beta*t_))
      double u_steady = a / b;
      double u_switch = u_steady * (1.0 - std::exp(-b * ts));
      double s_switch = 0.0;
      {
        double eu = std::exp(-b * ts), es = std::exp(-gm * ts);
        s_switch = (a / gm) * (1.0 - es) - (a / gm - u_switch) * (es - eu) * gm / (gm - b > 1e-10 ? (gm - b) : 1e-10);
        // Simplified: at switching time, s = alpha/gamma * (1 - exp(-gamma*t_))
        //   + correction from transient u
        s_switch = (a / gm) * (1.0 - es);
        if (std::abs(gm - b) > 1e-10) {
          double c_term = (a - u_switch * b) / (gm - b);
          s_switch += c_term * (eu - es);
        }
        if (!std::isfinite(s_switch)) s_switch = 0.0;
      }

      // Assign latent times
      std::vector<double> tau(n_cells);
      std::vector<int> regime(n_cells);  // 0=induction, 1=repression
      for (int c = 0; c < n_cells; ++c) {
        double uc = u[c];
        // If u < u_switch: induction phase
        if (uc <= u_switch * 1.01 && b > 1e-10) {
          regime[c] = 0;
          double arg = 1.0 - uc / u_steady;
          if (arg < 1e-12) arg = 1e-12;
          if (arg > 1.0) arg = 1.0;
          tau[c] = -std::log(arg) / b;
        } else if (b > 1e-10) {
          // Repression phase
          regime[c] = 1;
          if (u_steady > 1e-10 && uc < u_steady) {
            double arg = uc / u_steady;
            if (arg > 1.0) arg = 1.0;
            tau[c] = ts + (-std::log(arg)) / b;
          } else {
            tau[c] = ts * 3.0;
          }
        } else {
          regime[c] = 0;
          tau[c] = uc / (a > 1e-10 ? a : 1e-10);
        }
        if (!std::isfinite(tau[c]) || tau[c] < 0) tau[c] = 0.0;
      }

      // M-step: optimize alpha, beta, gamma given tau assignments
      // Objective: minimize sum of (u_obs - u_pred)^2 + (s_obs - s_pred)^2
      // Use gradient descent with line search
      double best_loss = 1e10;
      double best_a = a, best_b = b, best_g = gm, best_t = ts;

      // Multiple restarts for robustness
      for (int restart = 0; restart < 3; ++restart) {
        double ra = a, rb = b, rg = gm, rt = ts;
        if (restart == 1) { rb = 2.0; rg = gamma_det > 0 ? gamma_det * 2 : 1.0; }
        if (restart == 2) { rb = 0.5; rg = gamma_det > 0 ? gamma_det * 0.5 : 0.5; }

        double lr = 0.01;
        for (int gd_iter = 0; gd_iter < 200; ++gd_iter) {
          double loss = 0.0;
          double grad_a = 0.0, grad_b = 0.0, grad_g = 0.0;
          int used = 0;

          for (int c = 0; c < n_cells; ++c) {
            double ui = u[c], si = s[c], ti = tau[c];
            double u_pred, s_pred;

            if (ti <= rt || regime[c] == 0) {
              // Induction: u = alpha/beta * (1 - exp(-beta*tau))
              double eu = std::exp(-rb * ti);
              u_pred = ra / rb * (1.0 - eu);
              if (std::abs(rg - rb) > 1e-10) {
                double es = std::exp(-rg * ti);
                s_pred = (ra / rg) * (1.0 - es) + (ra / rb * eu - ra / rg) * (es - eu) * rg / (rg - rb);
              } else {
                s_pred = ra / rg * (1.0 - std::exp(-rg * ti));
              }
            } else {
              // Repression phase
              double u_s = ra / rb * (1.0 - std::exp(-rb * rt));
              double eu_t = std::exp(-rb * (ti - rt));
              u_pred = u_s * eu_t;
              double es_t = std::exp(-rg * (ti - rt));
              double s_at_t = (ra / rg) * (1.0 - std::exp(-rg * rt));
              s_pred = s_at_t + (u_s * rb / (rg - rb > 1e-10 ? (rg - rb) : 1e-10))
                        * (es_t - eu_t) + s_at_t * (1.0 - es_t);
              s_pred = s_at_t * es_t;  // simplified: s decays after switching
            }

            if (!std::isfinite(u_pred)) u_pred = 0.0;
            if (!std::isfinite(s_pred)) s_pred = 0.0;

            double err_u = ui - u_pred;
            double err_s = si - s_pred;
            loss += err_u * err_u + err_s * err_s;
            ++used;
          }

          // Numerical gradients
          double eps = 1e-5;
          auto compute_loss_at = [&](double pa, double pb, double pg, double pt) -> double {
            double l = 0.0; int n = 0;
            for (int c = 0; c < n_cells; ++c) {
              double ui = u[c], si = s[c], ti = tau[c];
              double up, sp;
              if (ti <= pt) {
                double eu = std::exp(-pb * ti);
                up = pa / pb * (1.0 - eu);
                sp = pa / pg * (1.0 - std::exp(-pg * ti));
              } else {
                double us = pa / pb * (1.0 - std::exp(-pb * pt));
                up = us * std::exp(-pb * (ti - pt));
                sp = pa / pg * (1.0 - std::exp(-pg * pt)) * std::exp(-pg * (ti - pt));
              }
              if (!std::isfinite(up)) up = 0.0;
              if (!std::isfinite(sp)) sp = 0.0;
              l += (ui - up) * (ui - up) + (si - sp) * (si - sp);
              ++n;
            }
            return n > 0 ? l / n : 1e10;
          };

          double ga = (compute_loss_at(ra + eps, rb, rg, rt) - compute_loss_at(ra - eps, rb, rg, rt)) / (2 * eps);
          double gb = (compute_loss_at(ra, rb + eps, rg, rt) - compute_loss_at(ra, rb - eps, rg, rt)) / (2 * eps);
          double gg = (compute_loss_at(ra, rb, rg + eps, rt) - compute_loss_at(ra, rb, rg - eps, rt)) / (2 * eps);

          ra -= lr * ga; rb -= lr * gb; rg -= lr * gg;
          ra = std::max(0.01, std::min(100.0, ra));
          rb = std::max(0.1, std::min(50.0, rb));
          rg = std::max(0.01, std::min(50.0, rg));

          double curr = compute_loss_at(ra, rb, rg, rt);
          if (curr < best_loss) {
            best_loss = curr;
            best_a = ra; best_b = rb; best_g = rg; best_t = rt;
          }
          if (gd_iter > 10 && std::abs(ga) + std::abs(gb) + std::abs(gg) < conv_tol) break;
        }
      }

      a = best_a; b = best_b; gm = best_g; ts = best_t;

      // Refine switching time
      double best_ts_loss = 1e10;
      for (double t_try = 0.5; t_try <= 10.0; t_try += 0.5) {
        double l = 0.0;
        for (int c = 0; c < n_cells; ++c) {
          double ui = u[c], si = s[c];
          double up, sp;
          double u_ss = a / b * (1.0 - std::exp(-b * t_try));
          if (ui <= u_ss * 1.01 && b > 1e-10) {
            double arg = 1.0 - ui / (a / b);
            if (arg < 1e-12) arg = 1e-12;
            if (arg > 1.0) arg = 1.0;
            up = ui;
            sp = a / gm * (1.0 - std::exp(-gm * (-std::log(arg) / b)));
          } else if (b > 1e-10) {
            up = u_ss * std::exp(-b * 0.5);  // approximate
            sp = a / gm * (1.0 - std::exp(-gm * t_try));
          } else { up = 0; sp = 0; }
          if (!std::isfinite(up)) up = 0.0;
          if (!std::isfinite(sp)) sp = 0.0;
          l += (ui - up) * (ui - up) + (si - sp) * (si - sp);
        }
        if (l < best_ts_loss) { best_ts_loss = l; ts = t_try; }
      }

      // Compute final loss
      double final_loss = 0.0;
      for (int c = 0; c < n_cells; ++c) {
        double ui = u[c], si = s[c];
        double u_st = a / b;
        double u_sw = u_st * (1.0 - std::exp(-b * ts));
        double up, sp;
        if (ui <= u_sw * 1.01) {
          up = ui;  // exact since we invert tau from u
          double tau_i = b > 1e-10 ? -std::log(std::max(1.0 - ui / u_st, 1e-12)) / b : 0.0;
          if (!std::isfinite(tau_i) || tau_i < 0) tau_i = 0.0;
          double es = std::exp(-gm * tau_i);
          sp = a / gm * (1.0 - es);
          if (std::abs(gm - b) > 1e-10) {
            sp += (a / b * std::exp(-b * tau_i) - a / gm) * (gm / (gm - b)) * (es - std::exp(-b * tau_i));
          }
        } else {
          double s_at_sw = a / gm * (1.0 - std::exp(-gm * ts));
          sp = s_at_sw;  // simplified
          up = u_sw;
        }
        if (!std::isfinite(up)) up = 0.0;
        if (!std::isfinite(sp)) sp = 0.0;
        final_loss += (ui - up) * (ui - up) + (si - sp) * (si - sp);
      }
      final_loss /= n_cells;

      if (std::abs(prev_loss - final_loss) < conv_tol) break;
      prev_loss = final_loss;
    }

    alpha_out[g] = a;
    beta_out[g] = b;
    gamma_out[g] = gm;
    t_out[g] = ts;
    scaling_out[g] = scaling;
    loss_out[g] = prev_loss;
    converged[g] = (prev_loss < 1e9) ? 1 : 0;
    ++n_fit;

    if (n_fit % 50 == 0) Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["alpha"] = alpha_out,
    _["beta"] = beta_out,
    _["gamma"] = gamma_out,
    _["t_"] = t_out,
    _["scaling"] = scaling_out,
    _["loss"] = loss_out,
    _["converged"] = converged,
    _["n_fitted"] = n_fit
  );
}
