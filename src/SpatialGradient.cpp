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
  std::vector<double> sig_p;
  std::vector<int> sig_n;
  std::vector<int> sig_n_nonzero;

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
    for (int col = 0; col < mat.n_cols; ++col) {
      if (std::isfinite(dist[col]) && std::isfinite(y_all[static_cast<std::size_t>(col)])) {
        x_valid.push_back(dist[col]);
        y_valid.push_back(y_all[static_cast<std::size_t>(col)]);
        if (y_all[static_cast<std::size_t>(col)] > 0.0) {
          ++n_nonzero;
        }
      }
    }
    if (n_nonzero < min_spots || x_valid.size() < 3) {
      continue;
    }

    FitStats fit = fit_linear(x_valid, y_valid, n_random, seed, row);
    fit.n_nonzero = n_nonzero;
    if (!std::isfinite(fit.r2)) {
      continue;
    }
    const std::string variable = as<std::string>(variables[row]);
    sig_variable.push_back(variable);
    sig_tot_var.push_back(fit.r2);
    sig_p.push_back(fit.p_value);
    sig_n.push_back(fit.n);
    sig_n_nonzero.push_back(fit.n_nonzero);

    fit_variable.push_back(variable);
    fit_model.push_back("linear");
    fit_mae.push_back(fit.mae);
    fit_rmse.push_back(fit.rmse);
    fit_r2.push_back(fit.r2);
    fit_slope.push_back(fit.slope);
    fit_intercept.push_back(fit.intercept);

    NumericVector bin_value_sum(bin_centers.size());
    IntegerVector bin_value_count(bin_centers.size());
    for (int col = 0; col < mat.n_cols; ++col) {
      if (bins[col] == NA_INTEGER) {
        continue;
      }
      const int bin = bins[col];
      bin_value_sum[bin] += y_all[static_cast<std::size_t>(col)];
      bin_value_count[bin] += 1;
    }
    for (int bin = 0; bin < bin_centers.size(); ++bin) {
      if (bin_value_count[bin] <= 0 || !std::isfinite(bin_centers[bin])) {
        continue;
      }
      sc_variable.push_back(variable);
      sc_distance.push_back(bin_centers[bin]);
      sc_value.push_back(bin_value_sum[bin] / bin_value_count[bin]);
      sc_estimate.push_back(fit.intercept + fit.slope * bin_centers[bin]);
      sc_reference.push_back(mode);
      sc_mode.push_back(mode);
    }
  }

  DataFrame significance = DataFrame::create(
    _["variable"] = sig_variable,
    _["tot_var"] = sig_tot_var,
    _["p_value"] = sig_p,
    _["fdr"] = NumericVector(sig_p.size(), NA_REAL),
    _["n_spots"] = sig_n,
    _["n_nonzero"] = sig_n_nonzero,
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
