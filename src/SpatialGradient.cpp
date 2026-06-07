// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

struct DgCSlots {
  int n_rows;
  int n_cols;
  IntegerVector row_idx;
  IntegerVector col_ptr;
  NumericVector values;
};

DgCSlots get_dgc_slots_spatial_gradient(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  IntegerVector row_idx = mat.slot("i");
  IntegerVector col_ptr = mat.slot("p");
  NumericVector values = mat.slot("x");
  if (dims.size() != 2 || col_ptr.size() != dims[1] + 1 || row_idx.size() != values.size()) {
    stop("expr must be a valid dgCMatrix.");
  }
  return DgCSlots{dims[0], dims[1], row_idx, col_ptr, values};
}

struct Entry {
  int col;
  double value;
};

struct FitStats {
  double intercept;
  double slope;
  double r;
  double r2;
  double mae;
  double rmse;
  double p_value;
  int n;
  int n_nonzero;
};

struct GradientStats {
  double tot_var;
  double norm_var;
  double rel_var;
  double p_value;
};

inline bool finite2(double x, double y) {
  return std::isfinite(x) && std::isfinite(y);
}

double pearson_r(const std::vector<double>& x, const std::vector<double>& y) {
  const int n = static_cast<int>(x.size());
  if (n < 3) {
    return NA_REAL;
  }
  double sx = 0.0, sy = 0.0, sx2 = 0.0, sy2 = 0.0, sxy = 0.0;
  for (int i = 0; i < n; ++i) {
    sx += x[i];
    sy += y[i];
    sx2 += x[i] * x[i];
    sy2 += y[i] * y[i];
    sxy += x[i] * y[i];
  }
  const double num = n * sxy - sx * sy;
  const double den_x = n * sx2 - sx * sx;
  const double den_y = n * sy2 - sy * sy;
  if (!(den_x > 0.0) || !(den_y > 0.0)) {
    return NA_REAL;
  }
  double r = num / std::sqrt(den_x * den_y);
  if (r > 1.0) {
    r = 1.0;
  } else if (r < -1.0) {
    r = -1.0;
  }
  return r;
}

FitStats fit_linear(
  const std::vector<double>& x,
  const std::vector<double>& y,
  int n_random,
  int seed,
  int feature_index
) {
  FitStats out;
  out.intercept = NA_REAL;
  out.slope = NA_REAL;
  out.r = NA_REAL;
  out.r2 = NA_REAL;
  out.mae = NA_REAL;
  out.rmse = NA_REAL;
  out.p_value = NA_REAL;
  out.n = static_cast<int>(x.size());
  out.n_nonzero = 0;

  const int n = out.n;
  if (n < 3) {
    return out;
  }

  double sx = 0.0, sy = 0.0, sx2 = 0.0, sxy = 0.0;
  for (int i = 0; i < n; ++i) {
    sx += x[i];
    sy += y[i];
    sx2 += x[i] * x[i];
    sxy += x[i] * y[i];
    if (std::isfinite(y[i]) && y[i] > 0.0) {
      ++out.n_nonzero;
    }
  }
  const double den = n * sx2 - sx * sx;
  if (!(den > 0.0)) {
    return out;
  }
  out.slope = (n * sxy - sx * sy) / den;
  out.intercept = (sy - out.slope * sx) / n;
  out.r = pearson_r(x, y);
  if (std::isfinite(out.r)) {
    out.r2 = out.r * out.r;
  }

  double abs_sum = 0.0;
  double sq_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double pred = out.intercept + out.slope * x[i];
    const double err = y[i] - pred;
    abs_sum += std::abs(err);
    sq_sum += err * err;
  }
  out.mae = abs_sum / n;
  out.rmse = std::sqrt(sq_sum / n);

  if (std::isfinite(out.r)) {
    if (n_random > 0) {
      std::vector<double> y_perm = y;
      std::mt19937 rng(static_cast<unsigned int>(seed + feature_index * 104729));
      const double obs = std::abs(out.r);
      int ge = 0;
      int used = 0;
      for (int iter = 0; iter < n_random; ++iter) {
        std::shuffle(y_perm.begin(), y_perm.end(), rng);
        const double rp = pearson_r(x, y_perm);
        if (std::isfinite(rp)) {
          ++used;
          if (std::abs(rp) >= obs) {
            ++ge;
          }
        }
      }
      if (used > 0) {
        out.p_value = (static_cast<double>(ge) + 1.0) / (static_cast<double>(used) + 1.0);
      }
    } else if (n > 2 && out.r2 < 1.0) {
      const double t = std::abs(out.r) * std::sqrt((n - 2.0) / std::max(1e-15, 1.0 - out.r2));
      out.p_value = 2.0 * R::pt(-t, static_cast<double>(n - 2), 1, 0);
    } else if (out.r2 >= 1.0) {
      out.p_value = 0.0;
    }
  }
  return out;
}

std::vector<double> normalize_minmax(const std::vector<double>& values) {
  double v_min = std::numeric_limits<double>::infinity();
  double v_max = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (std::isfinite(values[i])) {
      v_min = std::min(v_min, values[i]);
      v_max = std::max(v_max, values[i]);
    }
  }
  std::vector<double> out(values.size(), NA_REAL);
  if (!std::isfinite(v_min) || !std::isfinite(v_max) || !(v_max > v_min)) {
    return out;
  }
  const double denom = v_max - v_min;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (std::isfinite(values[i])) {
      out[i] = (values[i] - v_min) / denom;
    }
  }
  return out;
}

double total_variation_rescaled(const std::vector<double>& gradient) {
  std::vector<double> finite;
  finite.reserve(gradient.size());
  for (std::size_t i = 0; i < gradient.size(); ++i) {
    if (std::isfinite(gradient[i])) {
      finite.push_back(gradient[i]);
    }
  }
  if (finite.size() < 2) {
    return NA_REAL;
  }
  std::vector<double> scaled = normalize_minmax(finite);
  bool any_finite = false;
  for (std::size_t i = 0; i < scaled.size(); ++i) {
    if (std::isfinite(scaled[i])) {
      any_finite = true;
      break;
    }
  }
  if (!any_finite) {
    return 0.0;
  }
  double out = 0.0;
  for (std::size_t i = 1; i < scaled.size(); ++i) {
    if (std::isfinite(scaled[i]) && std::isfinite(scaled[i - 1])) {
      out += std::abs(scaled[i] - scaled[i - 1]);
    }
  }
  return out;
}

double relative_variation(const std::vector<double>& gradient) {
  double tv = total_variation_rescaled(gradient);
  if (!std::isfinite(tv)) {
    return NA_REAL;
  }
  double first = NA_REAL;
  double last = NA_REAL;
  for (std::size_t i = 0; i < gradient.size(); ++i) {
    if (std::isfinite(gradient[i])) {
      if (!std::isfinite(first)) {
        first = gradient[i];
      }
      last = gradient[i];
    }
  }
  if (!std::isfinite(first) || !std::isfinite(last)) {
    return NA_REAL;
  }
  const double net = std::abs(last - first);
  if (tv <= 0.0) {
    return 0.0;
  }
  return net / tv;
}

std::vector<double> local_linear_gradient(
  const std::vector<double>& x,
  const std::vector<double>& y,
  const NumericVector& centers,
  int n_bins
) {
  const int n = static_cast<int>(x.size());
  std::vector<double> out(static_cast<std::size_t>(centers.size()), NA_REAL);
  if (n < 3) {
    return out;
  }
  int k = static_cast<int>(std::ceil(static_cast<double>(n) / std::max(1, n_bins) * 3.0));
  k = std::max(3, std::min(n, k));
  std::vector<double> distances(static_cast<std::size_t>(n));

  for (int c = 0; c < centers.size(); ++c) {
    if (!std::isfinite(centers[c])) {
      continue;
    }
    const double center = centers[c];
    for (int i = 0; i < n; ++i) {
      distances[static_cast<std::size_t>(i)] = std::abs(x[static_cast<std::size_t>(i)] - center);
    }
    std::vector<double> sorted = distances;
    std::nth_element(sorted.begin(), sorted.begin() + k - 1, sorted.end());
    double bandwidth = sorted[static_cast<std::size_t>(k - 1)];
    if (!(bandwidth > 0.0)) {
      bandwidth = sorted.back();
    }
    if (!(bandwidth > 0.0)) {
      continue;
    }

    double sw = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < n; ++i) {
      const double d = distances[static_cast<std::size_t>(i)] / bandwidth;
      if (d > 1.0) {
        continue;
      }
      const double u = 1.0 - d * d * d;
      const double w = u * u * u;
      const double xc = x[static_cast<std::size_t>(i)] - center;
      const double yi = y[static_cast<std::size_t>(i)];
      sw += w;
      sx += w * xc;
      sy += w * yi;
      sxx += w * xc * xc;
      sxy += w * xc * yi;
    }
    if (!(sw > 0.0)) {
      continue;
    }
    const double den = sw * sxx - sx * sx;
    if (den > 0.0) {
      const double slope = (sw * sxy - sx * sy) / den;
      const double intercept = (sy - slope * sx) / sw;
      out[static_cast<std::size_t>(c)] = intercept;
    } else {
      out[static_cast<std::size_t>(c)] = sy / sw;
    }
  }
  return out;
}

std::vector<double> binned_means(
  const std::vector<double>& y_all,
  const IntegerVector& bins,
  int n_bins
) {
  std::vector<double> sums(static_cast<std::size_t>(n_bins), 0.0);
  std::vector<int> counts(static_cast<std::size_t>(n_bins), 0);
  for (int col = 0; col < bins.size(); ++col) {
    if (bins[col] == NA_INTEGER || !std::isfinite(y_all[static_cast<std::size_t>(col)])) {
      continue;
    }
    const int bin = bins[col];
    sums[static_cast<std::size_t>(bin)] += y_all[static_cast<std::size_t>(col)];
    counts[static_cast<std::size_t>(bin)] += 1;
  }
  std::vector<double> out(static_cast<std::size_t>(n_bins), NA_REAL);
  for (int bin = 0; bin < n_bins; ++bin) {
    if (counts[static_cast<std::size_t>(bin)] > 0) {
      out[static_cast<std::size_t>(bin)] =
        sums[static_cast<std::size_t>(bin)] / counts[static_cast<std::size_t>(bin)];
    }
  }
  return out;
}

GradientStats gradient_test(
  const std::vector<double>& x_valid,
  const std::vector<double>& gradient,
  const NumericVector& bin_centers,
  int n_bins,
  int n_random,
  int seed,
  int feature_index
) {
  GradientStats out;
  out.tot_var = total_variation_rescaled(gradient);
  out.norm_var = std::isfinite(out.tot_var) && gradient.size() > 0 ?
    out.tot_var / static_cast<double>(gradient.size()) : NA_REAL;
  out.rel_var = relative_variation(gradient);
  out.p_value = NA_REAL;
  if (!std::isfinite(out.tot_var) || n_random <= 0 || x_valid.size() < 3) {
    return out;
  }

  std::mt19937 rng(static_cast<unsigned int>(seed + feature_index * 104729));
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  int lt = 0;
  int used = 0;
  std::vector<double> random_y(x_valid.size(), 0.0);
  for (int iter = 0; iter < n_random; ++iter) {
    for (std::size_t i = 0; i < random_y.size(); ++i) {
      random_y[i] = unif(rng);
    }
    std::vector<double> random_gradient = local_linear_gradient(
      x_valid,
      random_y,
      bin_centers,
      n_bins
    );
    const double random_tv = total_variation_rescaled(random_gradient);
    if (std::isfinite(random_tv)) {
      ++used;
      if (random_tv < out.tot_var) {
        ++lt;
      }
    }
  }
  if (used > 0) {
    out.p_value = static_cast<double>(lt) / static_cast<double>(used);
  }
  return out;
}

void add_model_fit(
  const std::string& variable,
  const std::string& model,
  const std::vector<double>& gradient,
  const std::vector<double>& model_values,
  std::vector<std::string>& fit_variable,
  std::vector<std::string>& fit_model,
  std::vector<double>& fit_mae,
  std::vector<double>& fit_rmse,
  std::vector<double>& fit_r2,
  std::vector<double>& fit_slope,
  std::vector<double>& fit_intercept
) {
  double abs_sum = 0.0;
  double sq_sum = 0.0;
  std::vector<double> g;
  std::vector<double> m;
  const std::size_t n = std::min(gradient.size(), model_values.size());
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(gradient[i]) && std::isfinite(model_values[i])) {
      const double err = gradient[i] - model_values[i];
      abs_sum += std::abs(err);
      sq_sum += err * err;
      g.push_back(gradient[i]);
      m.push_back(model_values[i]);
    }
  }
  if (g.empty()) {
    return;
  }
  const double r = pearson_r(g, m);
  fit_variable.push_back(variable);
  fit_model.push_back(model);
  fit_mae.push_back(abs_sum / static_cast<double>(g.size()));
  fit_rmse.push_back(std::sqrt(sq_sum / static_cast<double>(g.size())));
  fit_r2.push_back(std::isfinite(r) ? r * r : NA_REAL);
  fit_slope.push_back(NA_REAL);
  fit_intercept.push_back(NA_REAL);
}

std::vector<double> annotation_distances(const NumericMatrix& coords, const LogicalVector& reference_spots) {
  const int n = coords.nrow();
  std::vector<int> refs;
  refs.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (reference_spots[i] == TRUE && finite2(coords(i, 0), coords(i, 1))) {
      refs.push_back(i);
    }
  }
  if (refs.empty()) {
    stop("The annotation reference contains no finite spots.");
  }
  std::vector<double> dist(n, NA_REAL);
  for (int i = 0; i < n; ++i) {
    if ((i & 255) == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (!finite2(coords(i, 0), coords(i, 1))) {
      continue;
    }
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < refs.size(); ++k) {
      const int j = refs[k];
      const double dx = coords(i, 0) - coords(j, 0);
      const double dy = coords(i, 1) - coords(j, 1);
      const double d2 = dx * dx + dy * dy;
      if (d2 < best) {
        best = d2;
      }
    }
    dist[i] = std::sqrt(best);
  }
  return dist;
}

std::vector<double> trajectory_distances(const NumericMatrix& coords, const NumericMatrix& trajectory) {
  const int n = coords.nrow();
  const int k = trajectory.nrow();
  if (k < 2 || trajectory.ncol() != 2) {
    stop("trajectory must contain at least two x/y coordinates.");
  }
  std::vector<double> seg_len(k - 1, 0.0);
  std::vector<double> cum(k - 1, 0.0);
  double total = 0.0;
  for (int s = 0; s < k - 1; ++s) {
    const double dx = trajectory(s + 1, 0) - trajectory(s, 0);
    const double dy = trajectory(s + 1, 1) - trajectory(s, 1);
    seg_len[s] = std::sqrt(dx * dx + dy * dy);
    cum[s] = total;
    total += seg_len[s];
  }
  if (!(total > 0.0)) {
    stop("trajectory length must be greater than zero.");
  }

  std::vector<double> dist(n, NA_REAL);
  for (int i = 0; i < n; ++i) {
    if ((i & 255) == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (!finite2(coords(i, 0), coords(i, 1))) {
      continue;
    }
    double best_d2 = std::numeric_limits<double>::infinity();
    double best_along = NA_REAL;
    for (int s = 0; s < k - 1; ++s) {
      if (!(seg_len[s] > 0.0)) {
        continue;
      }
      const double ax = trajectory(s, 0);
      const double ay = trajectory(s, 1);
      const double bx = trajectory(s + 1, 0);
      const double by = trajectory(s + 1, 1);
      const double vx = bx - ax;
      const double vy = by - ay;
      const double wx = coords(i, 0) - ax;
      const double wy = coords(i, 1) - ay;
      double t = (wx * vx + wy * vy) / (seg_len[s] * seg_len[s]);
      if (t < 0.0) {
        t = 0.0;
      } else if (t > 1.0) {
        t = 1.0;
      }
      const double px = ax + t * vx;
      const double py = ay + t * vy;
      const double dx = coords(i, 0) - px;
      const double dy = coords(i, 1) - py;
      const double d2 = dx * dx + dy * dy;
      if (d2 < best_d2) {
        best_d2 = d2;
        best_along = cum[s] + t * seg_len[s];
      }
    }
    dist[i] = best_along;
  }
  return dist;
}

IntegerVector make_bins(const std::vector<double>& dist, int n_bins, NumericVector& centers) {
  double d_min = std::numeric_limits<double>::infinity();
  double d_max = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < dist.size(); ++i) {
    if (std::isfinite(dist[i])) {
      d_min = std::min(d_min, dist[i]);
      d_max = std::max(d_max, dist[i]);
    }
  }
  if (!std::isfinite(d_min) || !std::isfinite(d_max)) {
    stop("No finite spatial distances are available.");
  }
  if (!(d_max > d_min)) {
    n_bins = 1;
  }
  IntegerVector bins(dist.size(), NA_INTEGER);
  NumericVector sums(n_bins);
  IntegerVector counts(n_bins);
  for (std::size_t i = 0; i < dist.size(); ++i) {
    if (!std::isfinite(dist[i])) {
      continue;
    }
    int bin = 0;
    if (n_bins > 1) {
      bin = static_cast<int>(std::floor((dist[i] - d_min) / (d_max - d_min) * n_bins));
      if (bin >= n_bins) {
        bin = n_bins - 1;
      }
      if (bin < 0) {
        bin = 0;
      }
    }
    bins[i] = bin;
    sums[bin] += dist[i];
    counts[bin] += 1;
  }
  centers = NumericVector(n_bins, NA_REAL);
  for (int bin = 0; bin < n_bins; ++bin) {
    if (counts[bin] > 0) {
      centers[bin] = sums[bin] / counts[bin];
    }
  }
  return bins;
}

} // namespace

// [[Rcpp::export]]
List spatial_gradient_screening_cpp(
  S4 expr,
  NumericMatrix coords,
  LogicalVector reference_spots,
  NumericMatrix trajectory,
  CharacterVector variables,
  std::string mode,
  int n_bins = 50,
  int n_random = 0,
  int seed = 123,
  int min_spots = 3
) {
  DgCSlots mat = get_dgc_slots_spatial_gradient(expr);
  if (coords.nrow() != mat.n_cols || coords.ncol() != 2) {
    stop("coords must have one x/y row per expression column.");
  }
  if (variables.size() != mat.n_rows) {
    stop("variables must contain one name per expression row.");
  }
  if (n_bins < 1) {
    stop("n_bins must be positive.");
  }
  if (n_random < 0) {
    stop("n_random must be non-negative.");
  }
  if (min_spots < 1) {
    stop("min_spots must be positive.");
  }

  std::vector<std::vector<Entry> > rows(static_cast<std::size_t>(mat.n_rows));
  for (int col = 0; col < mat.n_cols; ++col) {
    for (int ptr = mat.col_ptr[col]; ptr < mat.col_ptr[col + 1]; ++ptr) {
      const int row = mat.row_idx[ptr];
      if (row < 0 || row >= mat.n_rows) {
        stop("expr row index is out of bounds.");
      }
      const double value = mat.values[ptr];
      if (std::isfinite(value) && value != 0.0) {
        rows[static_cast<std::size_t>(row)].push_back(Entry{col, value});
      }
    }
  }

  std::vector<double> dist = mode == "annotation" ?
    annotation_distances(coords, reference_spots) :
    trajectory_distances(coords, trajectory);

  NumericVector bin_centers;
  IntegerVector bins = make_bins(dist, n_bins, bin_centers);
  int n_nonempty_bins = 0;
  for (int b = 0; b < bin_centers.size(); ++b) {
    if (std::isfinite(bin_centers[b])) {
      ++n_nonempty_bins;
    }
  }

  std::vector<std::string> sig_variable;
  std::vector<double> sig_tot_var;
  std::vector<double> sig_norm_var;
  std::vector<double> sig_rel_var;
  std::vector<double> sig_p;
  std::vector<int> sig_n;
  std::vector<int> sig_n_nonzero;
  std::vector<double> sig_linear_r2;
  std::vector<double> sig_linear_slope;
  std::vector<double> sig_linear_intercept;

  std::vector<std::string> fit_variable;
  std::vector<std::string> fit_model;
  std::vector<double> fit_mae;
  std::vector<double> fit_rmse;
  std::vector<double> fit_r2;
  std::vector<double> fit_slope;
  std::vector<double> fit_intercept;

  std::vector<std::string> sc_variable;
  std::vector<double> sc_distance;
  std::vector<double> sc_value;
  std::vector<double> sc_estimate;
  std::vector<std::string> sc_reference;
  std::vector<std::string> sc_mode;
  sc_variable.reserve(static_cast<std::size_t>(mat.n_rows) * std::max(1, n_nonempty_bins));

  std::vector<double> y_all(static_cast<std::size_t>(mat.n_cols), 0.0);
  std::vector<double> y_all_norm(static_cast<std::size_t>(mat.n_cols), NA_REAL);
  std::vector<double> x_valid;
  std::vector<double> y_valid;

  for (int row = 0; row < mat.n_rows; ++row) {
    if ((row & 31) == 0) {
      Rcpp::checkUserInterrupt();
    }
    std::fill(y_all.begin(), y_all.end(), 0.0);
    const std::vector<Entry>& entries = rows[static_cast<std::size_t>(row)];
    for (std::size_t k = 0; k < entries.size(); ++k) {
      y_all[static_cast<std::size_t>(entries[k].col)] = entries[k].value;
    }

    x_valid.clear();
    y_valid.clear();
    x_valid.reserve(mat.n_cols);
    y_valid.reserve(mat.n_cols);
    int n_nonzero = 0;
    double y_min = std::numeric_limits<double>::infinity();
    double y_max = -std::numeric_limits<double>::infinity();
    for (int col = 0; col < mat.n_cols; ++col) {
      if (std::isfinite(dist[col]) && std::isfinite(y_all[static_cast<std::size_t>(col)])) {
        const double y_raw = y_all[static_cast<std::size_t>(col)];
        y_min = std::min(y_min, y_raw);
        y_max = std::max(y_max, y_raw);
        if (y_all[static_cast<std::size_t>(col)] > 0.0) {
          ++n_nonzero;
        }
      }
    }
    if (n_nonzero < min_spots || !std::isfinite(y_min) || !std::isfinite(y_max) || !(y_max > y_min)) {
      continue;
    }
    const double y_den = y_max - y_min;
    std::fill(y_all_norm.begin(), y_all_norm.end(), NA_REAL);
    for (int col = 0; col < mat.n_cols; ++col) {
      if (std::isfinite(dist[col]) && std::isfinite(y_all[static_cast<std::size_t>(col)])) {
        const double y_norm = (y_all[static_cast<std::size_t>(col)] - y_min) / y_den;
        y_all_norm[static_cast<std::size_t>(col)] = y_norm;
        x_valid.push_back(dist[col]);
        y_valid.push_back(y_norm);
      }
    }
    if (x_valid.size() < 3) {
      continue;
    }

    FitStats fit = fit_linear(x_valid, y_valid, 0, seed, row);
    fit.n_nonzero = n_nonzero;
    std::vector<double> gradient = local_linear_gradient(
      x_valid,
      y_valid,
      bin_centers,
      n_bins
    );
    std::vector<double> gradient_scaled = normalize_minmax(gradient);
    GradientStats grad_stats = gradient_test(
      x_valid,
      gradient,
      bin_centers,
      n_bins,
      n_random,
      seed,
      row
    );
    if (!std::isfinite(grad_stats.tot_var)) {
      continue;
    }
    const std::string variable = as<std::string>(variables[row]);
    sig_variable.push_back(variable);
    sig_tot_var.push_back(grad_stats.tot_var);
    sig_norm_var.push_back(grad_stats.norm_var);
    sig_rel_var.push_back(grad_stats.rel_var);
    sig_p.push_back(grad_stats.p_value);
    sig_n.push_back(fit.n);
    sig_n_nonzero.push_back(fit.n_nonzero);
    sig_linear_r2.push_back(fit.r2);
    sig_linear_slope.push_back(fit.slope);
    sig_linear_intercept.push_back(fit.intercept);

    double d_min = std::numeric_limits<double>::infinity();
    double d_max = -std::numeric_limits<double>::infinity();
    for (int bin = 0; bin < bin_centers.size(); ++bin) {
      if (std::isfinite(bin_centers[bin])) {
        d_min = std::min(d_min, static_cast<double>(bin_centers[bin]));
        d_max = std::max(d_max, static_cast<double>(bin_centers[bin]));
      }
    }
    std::vector<double> model_ascending(static_cast<std::size_t>(bin_centers.size()), NA_REAL);
    std::vector<double> model_descending(static_cast<std::size_t>(bin_centers.size()), NA_REAL);
    std::vector<double> model_peak(static_cast<std::size_t>(bin_centers.size()), NA_REAL);
    std::vector<double> model_valley(static_cast<std::size_t>(bin_centers.size()), NA_REAL);
    std::vector<double> model_linear(static_cast<std::size_t>(bin_centers.size()), NA_REAL);
    if (std::isfinite(d_min) && std::isfinite(d_max) && d_max > d_min) {
      for (int bin = 0; bin < bin_centers.size(); ++bin) {
        if (!std::isfinite(bin_centers[bin])) {
          continue;
        }
        const double t = (bin_centers[bin] - d_min) / (d_max - d_min);
        model_ascending[static_cast<std::size_t>(bin)] = t;
        model_descending[static_cast<std::size_t>(bin)] = 1.0 - t;
        model_peak[static_cast<std::size_t>(bin)] = 1.0 - std::abs(2.0 * t - 1.0);
        model_valley[static_cast<std::size_t>(bin)] = std::abs(2.0 * t - 1.0);
        if (std::isfinite(fit.intercept) && std::isfinite(fit.slope)) {
          model_linear[static_cast<std::size_t>(bin)] = fit.intercept + fit.slope * bin_centers[bin];
        }
      }
      model_linear = normalize_minmax(model_linear);
    }
    add_model_fit(variable, "ascending", gradient_scaled, model_ascending, fit_variable, fit_model, fit_mae, fit_rmse, fit_r2, fit_slope, fit_intercept);
    add_model_fit(variable, "descending", gradient_scaled, model_descending, fit_variable, fit_model, fit_mae, fit_rmse, fit_r2, fit_slope, fit_intercept);
    add_model_fit(variable, "peak", gradient_scaled, model_peak, fit_variable, fit_model, fit_mae, fit_rmse, fit_r2, fit_slope, fit_intercept);
    add_model_fit(variable, "valley", gradient_scaled, model_valley, fit_variable, fit_model, fit_mae, fit_rmse, fit_r2, fit_slope, fit_intercept);
    add_model_fit(variable, "linear", gradient_scaled, model_linear, fit_variable, fit_model, fit_mae, fit_rmse, fit_r2, fit_slope, fit_intercept);

    std::vector<double> bin_values = binned_means(y_all_norm, bins, bin_centers.size());
    for (int bin = 0; bin < bin_centers.size(); ++bin) {
      if (!std::isfinite(bin_centers[bin]) || !std::isfinite(bin_values[static_cast<std::size_t>(bin)])) {
        continue;
      }
      sc_variable.push_back(variable);
      sc_distance.push_back(bin_centers[bin]);
      sc_value.push_back(bin_values[static_cast<std::size_t>(bin)]);
      sc_estimate.push_back(gradient_scaled[static_cast<std::size_t>(bin)]);
      sc_reference.push_back(mode);
      sc_mode.push_back(mode);
    }
  }

  DataFrame significance = DataFrame::create(
    _["variable"] = sig_variable,
    _["tot_var"] = sig_tot_var,
    _["p_value"] = sig_p,
    _["fdr"] = NumericVector(sig_p.size(), NA_REAL),
    _["norm_var"] = sig_norm_var,
    _["rel_var"] = sig_rel_var,
    _["n_spots"] = sig_n,
    _["n_nonzero"] = sig_n_nonzero,
    _["linear_r2"] = sig_linear_r2,
    _["linear_slope"] = sig_linear_slope,
    _["linear_intercept"] = sig_linear_intercept,
    _["stringsAsFactors"] = false
  );
  DataFrame model_fits = DataFrame::create(
    _["variable"] = fit_variable,
    _["model"] = fit_model,
    _["mae"] = fit_mae,
    _["rmse"] = fit_rmse,
    _["r2"] = fit_r2,
    _["slope"] = fit_slope,
    _["intercept"] = fit_intercept,
    _["stringsAsFactors"] = false
  );
  DataFrame screening = DataFrame::create(
    _["variable"] = sc_variable,
    _["distance"] = sc_distance,
    _["value"] = sc_value,
    _["estimate"] = sc_estimate,
    _["reference"] = sc_reference,
    _["mode"] = sc_mode,
    _["stringsAsFactors"] = false
  );

  return List::create(
    _["screening"] = screening,
    _["significance"] = significance,
    _["model_fits"] = model_fits
  );
}
