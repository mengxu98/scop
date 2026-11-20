// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <thisutils/log_message.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

inline bool velocity_value_missing(double v) {
  return R_IsNA(v) || R_IsNaN(v) || !std::isfinite(v);
}

struct VelocityNeighborDifferences {
  arma::mat positions;
  arma::vec times;
  arma::vec distances;
};

VelocityNeighborDifferences pseudotime_neighbor_differences(
    int i,
    const arma::mat& x_emb,
    const arma::vec& pseudotime,
    const IntegerMatrix& neighbors,
    int k_use
) {
  const int n_dims = x_emb.n_cols;
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
    return {};
  }

  VelocityNeighborDifferences differences{
    arma::mat(n_valid, n_dims),
    arma::vec(n_valid),
    arma::vec(n_valid)
  };

  double self_time = pseudotime(i);

  for (int j = 0; j < n_valid; ++j) {
    int nb = neighbors(i, j) - 1;
    differences.times(j) = pseudotime(nb) - self_time;

    double d2 = 0.0;
    for (int dim = 0; dim < n_dims; ++dim) {
      double pd = x_emb(nb, dim) - x_emb(i, dim);
      differences.positions(j, dim) = pd;
      d2 += pd * pd;
    }
    differences.distances(j) = std::sqrt(std::max(d2, 0.0));
  }

  for (int j = 0; j < n_valid; ++j) {
    if (differences.distances(j) < 1e-15) {
      differences.distances(j) = 1e-10;
    }
  }

  return differences;
}

arma::rowvec velocity_knn_single(
    int i,
    const arma::mat& x_emb,
    const arma::vec& pseudotime,
    const IntegerMatrix& neighbors,
    int k_use
) {
  const int n_dims = x_emb.n_cols;
  arma::rowvec velocity(n_dims, arma::fill::zeros);

  const VelocityNeighborDifferences differences = pseudotime_neighbor_differences(
    i, x_emb, pseudotime, neighbors, k_use
  );
  const int n_valid = differences.times.n_elem;
  if (n_valid == 0) {
    return velocity;
  }
  const arma::mat& pos_diff = differences.positions;
  const arma::vec& time_diff = differences.times;
  const arma::vec& dists = differences.distances;

  arma::vec weights(n_valid);
  for (int j = 0; j < n_valid; ++j) {
    double w = time_diff(j) / dists(j);
    if (velocity_value_missing(w)) {
      w = 0.0;
    }
    weights(j) = w;
  }

  for (int j = 0; j < n_valid; ++j) {
    double w = weights(j);
    for (int dim = 0; dim < n_dims; ++dim) {
      velocity(dim) += pos_diff(j, dim) * w;
    }
  }
  velocity /= static_cast<double>(n_valid);

  for (int dim = 0; dim < n_dims; ++dim) {
    if (velocity_value_missing(velocity(dim))) {
      velocity.zeros();
      return velocity;
    }
  }

  return velocity;
}

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

  const VelocityNeighborDifferences differences = pseudotime_neighbor_differences(
    i, x_emb, pseudotime, neighbors, k_use
  );
  const int n_valid = differences.times.n_elem;
  if (n_valid == 0) {
    return gradient;
  }
  const arma::mat& pos_diff = differences.positions;
  const arma::vec& time_diff = differences.times;
  const arma::vec& dists = differences.distances;

  double mean_dist = arma::mean(dists);
  double sigma = mean_dist * smooth;
  if (sigma < 1e-15) {
    sigma = 1e-10;
  }

  arma::vec weights(n_valid);
  double neg_half_sigma2_inv = -0.5 / (sigma * sigma);
  for (int j = 0; j < n_valid; ++j) {
    double gaussian = std::exp(dists(j) * dists(j) * neg_half_sigma2_inv);
    weights(j) = gaussian * time_diff(j);
  }

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

  for (int dim = 0; dim < n_dims; ++dim) {
    if (velocity_value_missing(gradient(dim))) {
      gradient.zeros();
      return gradient;
    }
  }

  return gradient;
}

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

template <typename Velocity>
NumericMatrix pseudotime_velocity_matrix(
    NumericMatrix& x_emb,
    NumericVector& pseudotime,
    const IntegerMatrix& neighbors,
    bool normalize,
    Velocity velocity
) {
  const int n_cells = x_emb.nrow();
  const int n_dims = x_emb.ncol();
  const int k_use = neighbors.ncol();

  if (pseudotime.size() != n_cells) {
    thisutils::log_message("pseudotime length must equal number of rows in x_emb", "error");
  }
  if (neighbors.nrow() != n_cells) {
    thisutils::log_message("neighbors must have the same number of rows as x_emb", "error");
  }

  arma::mat emb(x_emb.begin(), n_cells, n_dims, false);
  arma::vec pt(pseudotime.begin(), n_cells, false);

  arma::mat v_emb(n_cells, n_dims, arma::fill::zeros);

  for (int i = 0; i < n_cells; ++i) {
    v_emb.row(i) = velocity(i, emb, pt, neighbors, k_use);
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

}

// [[Rcpp::export]]
NumericMatrix pseudotime_velocity_knn(
    NumericMatrix x_emb,
    NumericVector pseudotime,
    IntegerMatrix neighbors,
    bool normalize = true
) {
  return pseudotime_velocity_matrix(x_emb, pseudotime, neighbors, normalize, velocity_knn_single);
}

// [[Rcpp::export]]
NumericMatrix pseudotime_velocity_gradient(
    NumericMatrix x_emb,
    NumericVector pseudotime,
    IntegerMatrix neighbors,
    double smooth = 0.5,
    bool normalize = true
) {
  return pseudotime_velocity_matrix(
    x_emb, pseudotime, neighbors, normalize,
    [smooth](int i, const arma::mat& emb, const arma::vec& pt, const IntegerMatrix& neighbors, int k_use) {
      return velocity_gradient_single(i, emb, pt, neighbors, k_use, smooth);
    }
  );
}
