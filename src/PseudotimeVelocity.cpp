// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

inline bool velocity_value_missing(double v) {
  return R_IsNA(v) || R_IsNaN(v) || !std::isfinite(v);
}

// Compute velocity for a single cell using KNN method.
// Returns a vector of length n_dims, or all zeros if computation fails.
arma::rowvec velocity_knn_single(
    int i,
    const arma::mat& x_emb,
    const arma::vec& pseudotime,
    const IntegerMatrix& neighbors,
    int k_use
) {
  const int n_dims = x_emb.n_cols;
  arma::rowvec velocity(n_dims, arma::fill::zeros);

  // Count valid neighbors (non-NA)
  int n_valid = 0;
  for (int j = 0; j < k_use; ++j) {
    int nb = neighbors(i, j);
    if (nb != NA_INTEGER && nb >= 1 && nb <= static_cast<int>(x_emb.n_rows)) {
      ++n_valid;
    } else {
      break;  // NA values fill the rest of the column
    }
  }

  if (n_valid == 0) {
    return velocity;
  }

  // Collect valid neighbor indices (0-based for Armadillo)
  std::vector<int> nb_idx;
  nb_idx.reserve(n_valid);
  for (int j = 0; j < n_valid; ++j) {
    nb_idx.push_back(neighbors(i, j) - 1);  // 1-based R -> 0-based C++
  }

  // Compute pos_diff and time_diff
  arma::mat pos_diff(n_valid, n_dims);
  arma::vec time_diff(n_valid);
  arma::vec dists(n_valid);

  double self_time = pseudotime(i);

  for (int j = 0; j < n_valid; ++j) {
    int nb = nb_idx[j];
    time_diff(j) = pseudotime(nb) - self_time;

    double d2 = 0.0;
    for (int dim = 0; dim < n_dims; ++dim) {
      double pd = x_emb(nb, dim) - x_emb(i, dim);
      pos_diff(j, dim) = pd;
      d2 += pd * pd;
    }
    dists(j) = std::sqrt(std::max(d2, 0.0));
  }

  // Clamp zero distances
  for (int j = 0; j < n_valid; ++j) {
    if (dists(j) < 1e-15) {
      dists(j) = 1e-10;
    }
  }

  // Weights = time_diff / dists
  arma::vec weights(n_valid);
  for (int j = 0; j < n_valid; ++j) {
    double w = time_diff(j) / dists(j);
    if (velocity_value_missing(w)) {
      w = 0.0;
    }
    weights(j) = w;
  }

  // velocity = colSums(pos_diff * weights) / n_valid
  for (int j = 0; j < n_valid; ++j) {
    double w = weights(j);
    for (int dim = 0; dim < n_dims; ++dim) {
      velocity(dim) += pos_diff(j, dim) * w;
    }
  }
  velocity /= static_cast<double>(n_valid);

  // Check for invalid values
  for (int dim = 0; dim < n_dims; ++dim) {
    if (velocity_value_missing(velocity(dim))) {
      velocity.zeros();
      return velocity;
    }
  }

  return velocity;
}

// Compute gradient for a single cell using Gaussian-weighted gradient method.
arma::rowvec velocity_gradient_single(
    int i,
    const arma::mat& x_emb,
    const arma::vec& pseudotime,
    const IntegerMatrix& neighbors,
    int k_use,
    double smooth
) {
  const int n_dims = x_emb.n_cols;
  arma::rowvec gradient(n_dims, arma::fill::zeros);

  // Count valid neighbors
  int n_valid = 0;
  for (int j = 0; j < k_use; ++j) {
    int nb = neighbors(i, j);
    if (nb != NA_INTEGER && nb >= 1 && nb <= static_cast<int>(x_emb.n_rows)) {
      ++n_valid;
    } else {
      break;
    }
  }

  if (n_valid == 0) {
    return gradient;
  }

  std::vector<int> nb_idx;
  nb_idx.reserve(n_valid);
  for (int j = 0; j < n_valid; ++j) {
    nb_idx.push_back(neighbors(i, j) - 1);
  }

  arma::mat pos_diff(n_valid, n_dims);
  arma::vec time_diff(n_valid);
  arma::vec dists(n_valid);

  double self_time = pseudotime(i);

  for (int j = 0; j < n_valid; ++j) {
    int nb = nb_idx[j];
    time_diff(j) = pseudotime(nb) - self_time;

    double d2 = 0.0;
    for (int dim = 0; dim < n_dims; ++dim) {
      double pd = x_emb(nb, dim) - x_emb(i, dim);
      pos_diff(j, dim) = pd;
      d2 += pd * pd;
    }
    dists(j) = std::sqrt(std::max(d2, 0.0));
  }

  // Clamp zero distances
  for (int j = 0; j < n_valid; ++j) {
    if (dists(j) < 1e-15) {
      dists(j) = 1e-10;
    }
  }

  // sigma = mean(dists) * smooth
  double mean_dist = arma::mean(dists);
  double sigma = mean_dist * smooth;
  if (sigma < 1e-15) {
    sigma = 1e-10;
  }

  // weights = exp(-dists^2 / (2 * sigma^2)) * time_diff
  arma::vec weights(n_valid);
  double neg_half_sigma2_inv = -0.5 / (sigma * sigma);
  for (int j = 0; j < n_valid; ++j) {
    double gaussian = std::exp(dists(j) * dists(j) * neg_half_sigma2_inv);
    weights(j) = gaussian * time_diff(j);
  }

  // gradient = colSums(pos_diff * weights) / sum(abs(weights) + 1e-10)
  double abs_sum = 0.0;
  for (int j = 0; j < n_valid; ++j) {
    abs_sum += std::abs(weights(j));
  }
  abs_sum += 1e-10;

  for (int j = 0; j < n_valid; ++j) {
    double w = weights(j);
    for (int dim = 0; dim < n_dims; ++dim) {
      gradient(dim) += pos_diff(j, dim) * w;
    }
  }
  gradient /= abs_sum;

  // Check for invalid values
  for (int dim = 0; dim < n_dims; ++dim) {
    if (velocity_value_missing(gradient(dim))) {
      gradient.zeros();
      return gradient;
    }
  }

  return gradient;
}

// Normalize velocity matrix: each row is divided by its norm and multiplied by mean norm.
// Rows with zero norm are left as zero.
void normalize_velocity_rows(arma::mat& v_emb) {
  const int n_cells = v_emb.n_rows;
  const int n_dims = v_emb.n_cols;

  arma::vec norms(n_cells);
  double sum_norm = 0.0;
  int n_nonzero = 0;

  for (int i = 0; i < n_cells; ++i) {
    double n2 = 0.0;
    for (int dim = 0; dim < n_dims; ++dim) {
      double val = v_emb(i, dim);
      n2 += val * val;
    }
    double n = std::sqrt(n2);
    norms(i) = n;
    if (n > 0) {
      sum_norm += n;
      ++n_nonzero;
    }
  }

  double mean_norm = (n_nonzero > 0) ? sum_norm / static_cast<double>(n_nonzero) : 1.0;

  for (int i = 0; i < n_cells; ++i) {
    double n = norms(i);
    if (n < 1e-15) {
      n = 1.0;
    }
    double scale = mean_norm / n;
    for (int dim = 0; dim < n_dims; ++dim) {
      v_emb(i, dim) *= scale;
    }
  }
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix pseudotime_velocity_knn(
    NumericMatrix x_emb,
    NumericVector pseudotime,
    IntegerMatrix neighbors,
    bool normalize = true
) {
  const int n_cells = x_emb.nrow();
  const int n_dims = x_emb.ncol();
  const int k_use = neighbors.ncol();

  if (pseudotime.size() != n_cells) {
    stop("pseudotime length must equal number of rows in x_emb");
  }
  if (neighbors.nrow() != n_cells) {
    stop("neighbors must have the same number of rows as x_emb");
  }

  // Convert to Armadillo
  arma::mat emb(x_emb.begin(), n_cells, n_dims, false);
  arma::vec pt(pseudotime.begin(), n_cells, false);

  arma::mat v_emb(n_cells, n_dims, arma::fill::zeros);

  for (int i = 0; i < n_cells; ++i) {
    v_emb.row(i) = velocity_knn_single(i, emb, pt, neighbors, k_use);
  }

  if (normalize) {
    normalize_velocity_rows(v_emb);
  }

  // Copy back to NumericMatrix
  NumericMatrix result(n_cells, n_dims);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_dims; ++j) {
      result(i, j) = v_emb(i, j);
    }
  }

  List dimnames = x_emb.attr("dimnames");
  if (dimnames.size() > 0 && !Rf_isNull(dimnames[0])) {
    result.attr("dimnames") = List::create(dimnames[0], R_NilValue);
  }

  return result;
}

// [[Rcpp::export]]
NumericMatrix pseudotime_velocity_gradient(
    NumericMatrix x_emb,
    NumericVector pseudotime,
    IntegerMatrix neighbors,
    double smooth = 0.5,
    bool normalize = true
) {
  const int n_cells = x_emb.nrow();
  const int n_dims = x_emb.ncol();
  const int k_use = neighbors.ncol();

  if (pseudotime.size() != n_cells) {
    stop("pseudotime length must equal number of rows in x_emb");
  }
  if (neighbors.nrow() != n_cells) {
    stop("neighbors must have the same number of rows as x_emb");
  }

  arma::mat emb(x_emb.begin(), n_cells, n_dims, false);
  arma::vec pt(pseudotime.begin(), n_cells, false);

  arma::mat v_emb(n_cells, n_dims, arma::fill::zeros);

  for (int i = 0; i < n_cells; ++i) {
    v_emb.row(i) = velocity_gradient_single(i, emb, pt, neighbors, k_use, smooth);
  }

  if (normalize) {
    normalize_velocity_rows(v_emb);
  }

  NumericMatrix result(n_cells, n_dims);
  for (int i = 0; i < n_cells; ++i) {
    for (int j = 0; j < n_dims; ++j) {
      result(i, j) = v_emb(i, j);
    }
  }

  List dimnames = x_emb.attr("dimnames");
  if (dimnames.size() > 0 && !Rf_isNull(dimnames[0])) {
    result.attr("dimnames") = List::create(dimnames[0], R_NilValue);
  }

  return result;
}
