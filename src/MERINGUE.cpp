#include <Rcpp.h>
#include <R_ext/Random.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

double rbits(int bits) {
  int_least64_t v = 0;
  for (int n = 0; n <= bits; n += 16) {
    int v1 = (int) std::floor(unif_rand() * 65536.0);
    v = 65536 * v + v1;
  }
  const int_least64_t one64 = 1L;
  return static_cast<double>(v & ((one64 << bits) - 1));
}

double r_unif_index(double dn, bool rounding) {
  if (rounding) {
    return std::floor(dn * unif_rand());
  }
  if (dn <= 0.0) {
    return 0.0;
  }
  int bits = (int) std::ceil(std::log2(dn));
  double dv;
  do {
    dv = rbits(bits);
  } while (dn <= dv);
  return dv;
}

double moran_i_centered(
    const double* values,
    double mean_z,
    int n,
    const std::vector<double>& row_std
) {
  std::vector<double> centered(n);
  long double v_l = 0.0L;
  for (int i = 0; i < n; ++i) {
    centered[i] = values[i] - mean_z;
    v_l += static_cast<long double>(centered[i]) * centered[i];
  }
  double v = static_cast<double>(v_l);
  if (v <= 0.0) {
    return NA_REAL;
  }
  long double S0_l = 0.0L;
  long double cv_l = 0.0L;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double wn = row_std[static_cast<std::size_t>(i) * n + j];
      S0_l += static_cast<long double>(wn);
      cv_l += static_cast<long double>(wn) * centered[i] * centered[j];
    }
  }
  double S0 = static_cast<double>(S0_l);
  double cv = static_cast<double>(cv_l);
  return (static_cast<double>(n) / S0) * (cv / v);
}

double moran_i_centered_edges(
    const double* values,
    double mean_z,
    int n,
    const std::vector<double>& row_std,
    const std::vector<int>& edges,
    double S0
) {
  std::vector<double> centered(n);
  long double v_l = 0.0L;
  for (int i = 0; i < n; ++i) {
    centered[i] = values[i] - mean_z;
    v_l += static_cast<long double>(centered[i]) * centered[i];
  }
  double v = static_cast<double>(v_l);
  if (v <= 0.0) {
    return NA_REAL;
  }
  long double cv_l = 0.0L;
  for (std::size_t e = 0; e < edges.size(); e += 2) {
    const int i = edges[e];
    const int j = edges[e + 1];
    double wn = row_std[static_cast<std::size_t>(i) * n + j];
    cv_l += static_cast<long double>(wn) * centered[i] * centered[j];
  }
  double cv = static_cast<double>(cv_l);
  return (static_cast<double>(n) / S0) * (cv / v);
}

double moran_i_centered_edges_fast(
    const double* values,
    double mean_z,
    int n,
    const std::vector<double>& row_std,
    const std::vector<int>& edges,
    double S0,
    std::vector<double>& centered,
    double* cv_out
) {
  double v = 0.0;
  for (int i = 0; i < n; ++i) {
    const double centered_i = values[i] - mean_z;
    centered[i] = centered_i;
    v += centered_i * centered_i;
  }
  if (v <= 0.0) {
    *cv_out = 0.0;
    return NA_REAL;
  }
  double cv = 0.0;
  const std::size_t n_edges = edges.size();
  for (std::size_t e = 0; e < n_edges; e += 2) {
    const int i = edges[e];
    const int j = edges[e + 1];
    cv += row_std[static_cast<std::size_t>(i) * n + j] *
      centered[i] * centered[j];
  }
  *cv_out = cv;
  return (static_cast<double>(n) / S0) * (cv / v);
}

std::vector<int> nonzero_edges(
    const std::vector<double>& row_std,
    int n
) {
  std::vector<int> edges;
  edges.reserve(static_cast<std::size_t>(n) * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (row_std[static_cast<std::size_t>(i) * n + j] != 0.0) {
        edges.push_back(i);
        edges.push_back(j);
      }
    }
  }
  return edges;
}

}

std::vector<double> row_standardize_weight(const NumericMatrix& weight, int n) {
  std::vector<double> row_std(static_cast<std::size_t>(n) * n);
  std::vector<double> rs(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n; ++j) {
      row_sum += weight(i, j);
    }
    rs[i] = (row_sum == 0.0) ? 1.0 : row_sum;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      row_std[static_cast<std::size_t>(i) * n + j] = weight(i, j) / rs[i];
    }
  }
  return row_std;
}

struct MoranWeightStats {
  double S0;
  double s1;
  double s2;
};

inline MoranWeightStats compute_moran_weight_stats(const std::vector<double>& row_std, int n) {
  long double s1_l = 0.0L;
  long double s2_l = 0.0L;
  long double S0_l = 0.0L;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double wn = row_std[static_cast<std::size_t>(i) * n + j];
      S0_l += static_cast<long double>(wn);
      double wn_t = row_std[static_cast<std::size_t>(j) * n + i];
      double t = wn + wn_t;
      s1_l += static_cast<long double>(t) * t;
    }
  }
  for (int i = 0; i < n; ++i) {
    double rs_i = 0.0;
    double cs_i = 0.0;
    for (int j = 0; j < n; ++j) {
      rs_i += row_std[static_cast<std::size_t>(i) * n + j];
      cs_i += row_std[static_cast<std::size_t>(j) * n + i];
    }
    double t = rs_i + cs_i;
    s2_l += static_cast<long double>(t) * t;
  }
  return MoranWeightStats{
    static_cast<double>(S0_l),
    0.5 * static_cast<double>(s1_l),
    static_cast<double>(s2_l)
  };
}

inline void moran_test_statistics_fast(
    const double* z,
    const std::vector<double>& row_std,
    int n,
    const std::string& alt_norm,
    const MoranWeightStats& ws,
    double& observed,
    double& expected,
    double& sd,
    double& p_value
) {
  std::vector<double> zc(n);
  long double sum_l = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum_l += static_cast<long double>(z[i]);
  }
  const double mean_z = static_cast<double>(sum_l / static_cast<long double>(n));
  for (int i = 0; i < n; ++i) {
    zc[i] = z[i] - mean_z;
  }

  observed = moran_i_centered(zc.data(), 0.0, n, row_std);
  expected = -1.0 / static_cast<double>(n - 1);

  long double z2_l = 0.0L;
  long double z4_l = 0.0L;
  for (int i = 0; i < n; ++i) {
    z2_l += static_cast<long double>(zc[i]) * zc[i];
    long double z4 = static_cast<long double>(zc[i]) * zc[i];
    z4_l += z4 * z4;
  }

  const double S0 = ws.S0;
  const double s1 = ws.s1;
  const double s2 = ws.s2;
  const double n_d = static_cast<double>(n);
  const double k = (static_cast<double>(z4_l) / n_d) /
    std::pow(static_cast<double>(z2_l) / n_d, 2.0);
  const double vi =
    (n_d * ((std::pow(n_d, 2.0) - 3.0 * n_d + 3.0) * s1 - n_d * s2 + 3.0 * std::pow(S0, 2.0)) -
     k * (n_d * (n_d - 1.0) * s1 - 2.0 * n_d * s2 + 6.0 * std::pow(S0, 2.0))) /
      ((n_d - 1.0) * (n_d - 2.0) * (n_d - 3.0) * std::pow(S0, 2.0)) -
    std::pow(expected, 2.0);
  sd = std::sqrt(vi);

  double pv = R::pnorm(observed, expected, sd, TRUE, FALSE);
  if (alt_norm == "two.sided") {
    p_value = (observed <= expected) ? 2.0 * pv : 2.0 * (1.0 - pv);
  } else if (alt_norm == "greater") {
    p_value = 1.0 - pv;
  } else {
    p_value = pv;
  }
}

void moran_test_statistics(
    const double* z,
    const std::vector<double>& row_std,
    int n,
    const std::string& alt_norm,
    double& observed,
    double& expected,
    double& sd,
    double& p_value
) {
  MoranWeightStats ws = compute_moran_weight_stats(row_std, n);
  moran_test_statistics_fast(z, row_std, n, alt_norm, ws, observed, expected, sd, p_value);
}

// [[Rcpp::export]]
NumericMatrix meringue_moran_matrix_cpp(
    NumericMatrix expr,
    NumericMatrix weight,
    String alternative,
    bool rounding_sample = false,
    int n_perm = 0,
    int n_threads = 1
) {
  Rcpp::RNGScope rng_scope;
  const int n_genes = expr.nrow();
  const int n = expr.ncol();
  const std::string alt = alternative;

  if (weight.nrow() != n || weight.ncol() != n) {
    stop("weight must be a square matrix matching the expression columns");
  }
  if (n_perm < 0) {
    stop("n_perm must be non-negative");
  }
  n_threads = std::max(1, std::min(n_threads, std::max(1, n_genes)));

  std::string alt_norm;
  if (alt == "greater" || alt == "less" || alt == "two.sided") {
    alt_norm = alt;
  } else {
    alt_norm = "greater";
  }

  const std::vector<double> row_std = row_standardize_weight(weight, n);
  if (n_perm <= 0) {
    NumericMatrix out(n_genes, 4);
    colnames(out) = CharacterVector::create("observed", "expected", "sd", "p_value");
    rownames(out) = rownames(expr);

    MoranWeightStats ws = compute_moran_weight_stats(row_std, n);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
    {
      std::vector<double> z(n);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int g = 0; g < n_genes; ++g) {
        for (int j = 0; j < n; ++j) {
          z[j] = expr(g, j);
        }
        double observed, expected, sd, p_value;
        moran_test_statistics_fast(z.data(), row_std, n, alt_norm, ws, observed, expected, sd, p_value);
        out(g, 0) = observed;
        out(g, 1) = expected;
        out(g, 2) = sd;
        out(g, 3) = p_value;
      }
    }
    return out;
  }

  long double S0_l = 0.0L;
  for (std::size_t k = 0; k < row_std.size(); ++k) {
    S0_l += static_cast<long double>(row_std[k]);
  }
  const double S0 = static_cast<double>(S0_l);
  const std::vector<int> edges = nonzero_edges(row_std, n);

  const double n_d = static_cast<double>(n);
  std::vector<int> sample_indices(static_cast<std::size_t>(n_perm) * n);
  for (int p = 0; p < n_perm; ++p) {
    for (int i = 0; i < n; ++i) {
      sample_indices[static_cast<std::size_t>(p) * n + i] =
        static_cast<int>(r_unif_index(n_d, rounding_sample));
    }
  }

  std::vector<double> expr_rows(static_cast<std::size_t>(n_genes) * n);
  for (int g = 0; g < n_genes; ++g) {
    for (int j = 0; j < n; ++j) {
      expr_rows[static_cast<std::size_t>(g) * n + j] = expr(g, j);
    }
  }
  const int n_out = 5;
  std::vector<double> results(
    static_cast<std::size_t>(n_genes) * n_out,
    NA_REAL
  );

  #ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
  #endif
  for (int g = 0; g < n_genes; ++g) {
    const double* z = &expr_rows[static_cast<std::size_t>(g) * n];
    std::vector<double> zc(n);
    long double sum_l = 0.0L;
    for (int i = 0; i < n; ++i) {
      sum_l += static_cast<long double>(z[i]);
    }
    const double mean_z = static_cast<double>(sum_l / static_cast<long double>(n));
    for (int i = 0; i < n; ++i) {
      zc[i] = z[i] - mean_z;
    }
    const double observed = moran_i_centered(zc.data(), 0.0, n, row_std);

    std::vector<double> sim;
    sim.reserve(static_cast<std::size_t>(n_perm));
    std::vector<double> foo(n);
    std::vector<double> centered(n);
    double cv = 0.0;
    for (int p = 0; p < n_perm; ++p) {
      long double foo_sum = 0.0L;
      const std::size_t offset = static_cast<std::size_t>(p) * n;
      for (int i = 0; i < n; ++i) {
        foo[i] = z[sample_indices[offset + i]];
        foo_sum += static_cast<long double>(foo[i]);
      }
      const double foo_mean = static_cast<double>(
        foo_sum / static_cast<long double>(n)
      );
      const double I_foo = moran_i_centered_edges_fast(
        foo.data(), foo_mean, n, row_std, edges, S0, centered, &cv
      );
      if (!std::isnan(I_foo)) {
        sim.push_back(I_foo);
      }
    }

    const int n_sim = static_cast<int>(sim.size());
    double expected = NA_REAL;
    double sd = NA_REAL;
    double p_value = NA_REAL;
    if (n_sim > 0) {
      long double sim_sum = 0.0L;
      for (int i = 0; i < n_sim; ++i) {
        sim_sum += static_cast<long double>(sim[i]);
      }
      expected = static_cast<double>(sim_sum / static_cast<long double>(n_sim));

      if (n_sim > 1) {
        long double dev_sum = 0.0L;
        for (int i = 0; i < n_sim; ++i) {
          const double d = sim[i] - expected;
          dev_sum += static_cast<long double>(d) * d;
        }
        sd = std::sqrt(
          static_cast<double>(dev_sum) / static_cast<double>(n_sim - 1)
        );
      }

      int count = 0;
      for (int i = 0; i < n_sim; ++i) {
        if (alt_norm == "two.sided") {
          if (std::abs(sim[i]) >= std::abs(observed)) {
            ++count;
          }
        } else if (alt_norm == "greater") {
          if (sim[i] >= observed) {
            ++count;
          }
        } else if (sim[i] <= observed) {
          ++count;
        }
      }
      p_value = static_cast<double>(count + 1) /
        static_cast<double>(n_sim + 1);
    }

    const std::size_t out_offset = static_cast<std::size_t>(g) * n_out;
    results[out_offset] = observed;
    results[out_offset + 1] = expected;
    results[out_offset + 2] = sd;
    results[out_offset + 3] = p_value;
    results[out_offset + 4] = static_cast<double>(n_perm);
  }

  NumericMatrix out(n_genes, n_out);
  for (int g = 0; g < n_genes; ++g) {
    for (int j = 0; j < n_out; ++j) {
      out(g, j) = results[static_cast<std::size_t>(g) * n_out + j];
    }
  }
  colnames(out) = CharacterVector::create(
    "observed", "expected", "sd", "p_value", "N"
  );
  rownames(out) = rownames(expr);
  return out;
}

// [[Rcpp::export]]
NumericVector meringue_moran_cpp(
    NumericVector z,
    NumericMatrix weight,
    int n_perm,
    String alternative,
    bool rounding_sample = false
) {
  Rcpp::RNGScope rng_scope;
  const int n = z.size();
  const std::string alt = alternative;

  std::string alt_norm;
  if (alt == "greater" || alt == "less" || alt == "two.sided") {
    alt_norm = alt;
  } else {
    alt_norm = "greater";
  }

  const std::vector<double> row_std = row_standardize_weight(weight, n);

  std::vector<double> zc(n);
  long double sum_l = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum_l += static_cast<long double>(z[i]);
  }
  const double mean_z = static_cast<double>(sum_l / static_cast<long double>(n));
  for (int i = 0; i < n; ++i) {
    zc[i] = z[i] - mean_z;
  }

  double observed = moran_i_centered(zc.data(), 0.0, n, row_std);

  if (n_perm <= 0) {
    double expected, sd, p_value;
    moran_test_statistics(&z[0], row_std, n, alt_norm, observed, expected, sd, p_value);
    NumericVector out(4);
    out.attr("names") = CharacterVector::create("observed", "expected", "sd", "p_value");
    out[0] = observed;
    out[1] = expected;
    out[2] = sd;
    out[3] = p_value;
    return out;
  }

  const double n_d = static_cast<double>(n);
  std::vector<double> sim;
  sim.reserve(static_cast<std::size_t>(n_perm));
  std::vector<double> foo(n);
  std::vector<double> centered(n);

  long double S0_l = 0.0L;
  for (std::size_t k = 0; k < row_std.size(); ++k) {
    S0_l += static_cast<long double>(row_std[k]);
  }
  const double S0 = static_cast<double>(S0_l);
  const std::vector<int> edges = nonzero_edges(row_std, n);

  double cv = 0.0;
  for (int p = 0; p < n_perm; ++p) {
    for (int i = 0; i < n; ++i) {
      foo[i] = z[static_cast<int>(r_unif_index(n_d, rounding_sample))];
    }
    long double foo_sum = 0.0L;
    for (int i = 0; i < n; ++i) {
      foo_sum += static_cast<long double>(foo[i]);
    }
    const double foo_mean = static_cast<double>(foo_sum / static_cast<long double>(n));
    double I_foo = moran_i_centered_edges_fast(
      foo.data(), foo_mean, n, row_std, edges, S0, centered, &cv
    );
    if (!R_IsNaN(I_foo)) {
      sim.push_back(I_foo);
    }
  }

  const int n_sim = static_cast<int>(sim.size());
  long double sim_sum = 0.0L;
  for (int i = 0; i < n_sim; ++i) {
    sim_sum += static_cast<long double>(sim[i]);
  }
  const double expected = static_cast<double>(sim_sum / static_cast<long double>(n_sim));

  long double dev_sum = 0.0L;
  for (int i = 0; i < n_sim; ++i) {
    double d = sim[i] - expected;
    dev_sum += static_cast<long double>(d) * d;
  }
  const double sd = std::sqrt(static_cast<double>(dev_sum) / static_cast<double>(n_sim - 1));

  int count = 0;
  for (int i = 0; i < n_sim; ++i) {
    if (alt_norm == "two.sided") {
      if (std::abs(sim[i]) >= std::abs(observed)) {
        ++count;
      }
    } else if (alt_norm == "greater") {
      if (sim[i] >= observed) {
        ++count;
      }
    } else {
      if (sim[i] <= observed) {
        ++count;
      }
    }
  }
  const double p_value = static_cast<double>(count + 1) / static_cast<double>(n_sim + 1);

  NumericVector out(5);
  out.attr("names") = CharacterVector::create("observed", "expected", "sd", "p_value", "N");
  out[0] = observed;
  out[1] = expected;
  out[2] = sd;
  out[3] = p_value;
  out[4] = static_cast<double>(n_perm);
  return out;
}
