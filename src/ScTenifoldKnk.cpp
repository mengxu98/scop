#include <RcppArmadillo.h>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RSpectra)]]

class SctenifoldDowndatedSymMatProd {
private:
  typedef Eigen::Index Index;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Map<const Vector> MapConstVec;
  typedef Eigen::Map<Vector> MapVec;

  const Matrix& covariance_;
  const Vector& y_;

public:
  SctenifoldDowndatedSymMatProd(const Matrix& covariance, const Vector& y) :
    covariance_(covariance), y_(y) {}

  Index rows() const { return covariance_.rows(); }
  Index cols() const { return covariance_.cols(); }

  void perform_op(const double* x_in, double* y_out) const {
    MapConstVec x(x_in, covariance_.cols());
    MapVec out(y_out, covariance_.rows());
    out.noalias() = covariance_ * x;
    out.noalias() -= y_ * y_.dot(x);
  }
};

class SctenifoldEigenThreadGuard {
private:
  int old_threads_;
  bool active_;

public:
  explicit SctenifoldEigenThreadGuard(int n_threads) :
    old_threads_(Eigen::nbThreads()), active_(n_threads > 0) {
    if (active_) {
      Eigen::setNbThreads(n_threads);
    }
  }

  ~SctenifoldEigenThreadGuard() {
    if (active_) {
      Eigen::setNbThreads(old_threads_);
    }
  }
};

static int sctenifold_worker_count(int requested, int n_tasks) {
  if (requested <= 1 || n_tasks <= 1) {
    return 1;
  }

  int n_threads = std::min(requested, n_tasks);
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    n_threads = std::min(n_threads, static_cast<int>(hardware));
  }
  return std::max(1, n_threads);
}

static arma::mat sctenifold_crossprod(const arma::mat& x) {
  return x.t() * x;
}

static arma::mat sctenifold_safe_inv(const arma::mat& x) {
  arma::mat out;
  if (!arma::inv(out, x)) {
    out = arma::pinv(x);
  }
  return out;
}

static double sctenifold_tensor_residual(
  const std::vector<arma::mat>& tensors,
  const arma::mat& u1,
  const arma::mat& u2,
  const arma::mat& u3,
  const arma::mat& u4,
  const arma::rowvec& lambdas
) {
  const arma::uword n1 = u1.n_rows;
  const arma::uword n2 = u2.n_rows;
  const arma::uword n4 = u4.n_rows;
  const arma::uword k = u1.n_cols;
  double ss = 0.0;

  for (arma::uword l = 0; l < n4; ++l) {
    const arma::mat& x = tensors[l];
    for (arma::uword j = 0; j < n2; ++j) {
      for (arma::uword i = 0; i < n1; ++i) {
        double est = 0.0;
        for (arma::uword r = 0; r < k; ++r) {
          est += lambdas[r] * u1(i, r) * u2(j, r) * u3(0, r) * u4(l, r);
        }
        const double diff = est - x(i, j);
        ss += diff * diff;
      }
    }
  }

  return std::sqrt(ss);
}

static arma::mat sctenifold_tensor_mttkrp(
  const std::vector<arma::mat>& tensors,
  const arma::mat& u1,
  const arma::mat& u2,
  const arma::mat& u3,
  const arma::mat& u4,
  const int mode
) {
  const arma::uword n_genes = u1.n_rows;
  const arma::uword n_net = u4.n_rows;
  const arma::uword k = u1.n_cols;

  if (mode == 1) {
    arma::mat out(n_genes, k, arma::fill::zeros);
    for (arma::uword l = 0; l < n_net; ++l) {
      arma::mat prod = tensors[l] * u2;
      prod.each_row() %= (u3.row(0) % u4.row(l));
      out += prod;
    }
    return out;
  }

  if (mode == 2) {
    arma::mat out(n_genes, k, arma::fill::zeros);
    for (arma::uword l = 0; l < n_net; ++l) {
      arma::mat prod = tensors[l].t() * u1;
      prod.each_row() %= (u3.row(0) % u4.row(l));
      out += prod;
    }
    return out;
  }

  if (mode == 3) {
    arma::mat out(1, k, arma::fill::zeros);
    for (arma::uword l = 0; l < n_net; ++l) {
      arma::rowvec diag_vals = arma::diagvec(u1.t() * tensors[l] * u2).t();
      out.row(0) += u4.row(l) % diag_vals;
    }
    return out;
  }

  arma::mat out(n_net, k, arma::fill::zeros);
  for (arma::uword l = 0; l < n_net; ++l) {
    arma::rowvec diag_vals = arma::diagvec(u1.t() * tensors[l] * u2).t();
    out.row(l) = u3.row(0) % diag_vals;
  }
  return out;
}

static arma::mat sctenifold_tensor_average(
  const arma::uword n_genes,
  const arma::uword n_net,
  const arma::mat& u1,
  const arma::mat& u2,
  const arma::mat& u3,
  const arma::mat& u4,
  const arma::rowvec& lambdas
) {
  const arma::uword k = u1.n_cols;
  arma::mat out(n_genes, n_genes, arma::fill::zeros);

  for (arma::uword l = 0; l < n_net; ++l) {
    for (arma::uword r = 0; r < k; ++r) {
      out += (lambdas[r] * u3(0, r) * u4(l, r)) *
        (u1.col(r) * u2.col(r).t());
    }
  }
  out /= static_cast<double>(n_net);

  const double max_abs = arma::abs(out).max();
  if (max_abs > 0.0 && std::isfinite(max_abs)) {
    out /= max_abs;
  }
  return out;
}

static Eigen::MatrixXd sctenifold_pcnet_covariance_raw_eigen(
  NumericMatrix x,
  int n_comp,
  int ncv,
  int maxit,
  double tol,
  int n_threads
) {
  const int n_genes = x.nrow();
  const int n_cells = x.ncol();
  if (n_cells < 3) {
    stop("x must contain at least three cells");
  }
  if (n_comp < 2 || n_comp >= n_genes || n_comp >= n_cells) {
    stop("n_comp must be >= 2 and lower than the number of genes and cells");
  }

  Eigen::MatrixXd scaled(n_cells, n_genes);
  for (int g = 0; g < n_genes; ++g) {
    double mean = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      mean += x(g, c);
    }
    mean /= static_cast<double>(n_cells);

    double ss = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      const double centered = x(g, c) - mean;
      ss += centered * centered;
    }
    const double scale = std::sqrt(ss / static_cast<double>(n_cells - 1));
    if (scale == 0.0 || !std::isfinite(scale)) {
      stop("all retained genes must have finite non-zero variance");
    }
    for (int c = 0; c < n_cells; ++c) {
      scaled(c, g) = (x(g, c) - mean) / scale;
    }
  }

  const Eigen::MatrixXd covariance = scaled * scaled.transpose();
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_genes, n_genes);
  const int ncv_eff = ncv > n_comp ? std::min(ncv, n_cells) :
    std::min(n_cells, std::max(2 * n_comp + 1, 20));

  auto compute_gene = [&] (int g) {
    const Eigen::VectorXd y = scaled.col(g);
    SctenifoldDowndatedSymMatProd op(covariance, y);
    Spectra::SymEigsSolver<
      double,
      Spectra::LARGEST_ALGE,
      SctenifoldDowndatedSymMatProd
    > eigs(&op, n_comp, ncv_eff);
    eigs.init();
    const int nconv = eigs.compute(maxit, tol, Spectra::LARGEST_ALGE);
    if (eigs.info() != Spectra::SUCCESSFUL || nconv < n_comp) {
      throw std::runtime_error(
        "Spectra failed to converge while constructing pcNet at gene index " +
          std::to_string(g + 1)
      );
    }

    const Eigen::VectorXd values = eigs.eigenvalues();
    const Eigen::MatrixXd vectors = eigs.eigenvectors();
    Eigen::VectorXd projection = vectors.transpose() * y;
    projection.array() /= values.array();
    const Eigen::VectorXd score_projection = vectors * projection;
    const Eigen::VectorXd beta = scaled.transpose() * score_projection;
    for (int j = 0; j < n_genes; ++j) {
      if (j != g) {
        out(g, j) = beta(j);
      }
    }
  };

  const int n_workers = sctenifold_worker_count(n_threads, n_genes);
  if (n_workers == 1) {
    for (int g = 0; g < n_genes; ++g) {
      compute_gene(g);
    }
  } else {
    SctenifoldEigenThreadGuard eigen_threads(1);
    std::atomic<int> next_gene(0);
    std::atomic<int> failed(0);
    std::vector<std::string> errors(static_cast<std::size_t>(n_workers));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(n_workers));

    for (int worker = 0; worker < n_workers; ++worker) {
      workers.emplace_back([&, worker] {
        while (failed.load(std::memory_order_relaxed) == 0) {
          const int g = next_gene.fetch_add(1, std::memory_order_relaxed);
          if (g >= n_genes) {
            break;
          }
          try {
            compute_gene(g);
          } catch (const std::exception& e) {
            errors[static_cast<std::size_t>(worker)] = e.what();
            failed.store(1, std::memory_order_relaxed);
            break;
          } catch (...) {
            errors[static_cast<std::size_t>(worker)] =
              "Unknown error while constructing pcNet";
            failed.store(1, std::memory_order_relaxed);
            break;
          }
        }
      });
    }

    for (std::thread& worker : workers) {
      worker.join();
    }

    if (failed.load(std::memory_order_relaxed) != 0) {
      for (const std::string& error : errors) {
        if (!error.empty()) {
          stop(error);
        }
      }
      stop("Unknown error while constructing pcNet");
    }
  }

  return out;
}

static double sctenifold_order_statistic(std::vector<double>& values, std::size_t k) {
  std::nth_element(values.begin(), values.begin() + k, values.end());
  return values[k];
}

static double sctenifold_quantile_type7(std::vector<double>& values, double prob) {
  if (values.empty()) {
    stop("cannot compute quantile of an empty vector");
  }
  if (prob <= 0.0) {
    return sctenifold_order_statistic(values, 0);
  }
  if (prob >= 1.0) {
    return sctenifold_order_statistic(values, values.size() - 1);
  }
  const double index = 1.0 + (static_cast<double>(values.size()) - 1.0) * prob;
  const std::size_t lo = static_cast<std::size_t>(std::floor(index));
  const double gamma = index - static_cast<double>(lo);
  if (lo >= values.size()) {
    return sctenifold_order_statistic(values, values.size() - 1);
  }

  const std::size_t left_index = lo - 1;
  const double left = sctenifold_order_statistic(values, left_index);
  if (gamma == 0.0 || left_index + 1 >= values.size()) {
    return left;
  }
  const double right = sctenifold_order_statistic(values, left_index + 1);
  return left + gamma * (right - left);
}

static S4 sctenifold_dgCMatrix_from_dense_threshold(
  const Eigen::MatrixXd& coefficients,
  bool scale_scores,
  bool symmetric,
  double q
) {
  const int n = coefficients.rows();
  if (n != coefficients.cols()) {
    stop("coefficients must be a square matrix");
  }
  if (q < 0.0 || q > 1.0 || !std::isfinite(q)) {
    stop("q must be a finite number between 0 and 1");
  }
  Eigen::MatrixXd work = symmetric ?
    ((coefficients + coefficients.transpose()) / 2.0).eval() :
    coefficients;

  std::vector<double> abs_values;
  abs_values.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
  double max_abs = 0.0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const double abs_value = std::abs(work(i, j));
      abs_values.push_back(abs_value);
      if (abs_value > max_abs) {
        max_abs = abs_value;
      }
    }
  }
  const double cutoff = sctenifold_quantile_type7(abs_values, q);
  const double scale = scale_scores ? max_abs : 1.0;

  std::vector<int> i_slot;
  std::vector<int> p_slot;
  std::vector<double> x_slot;
  p_slot.reserve(static_cast<std::size_t>(n) + 1);
  p_slot.push_back(0);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (i == j) {
        continue;
      }
      const double raw = work(i, j);
      if (std::abs(raw) < cutoff || raw == 0.0) {
        continue;
      }
      i_slot.push_back(i);
      x_slot.push_back(scale_scores ? raw / scale : raw);
    }
    p_slot.push_back(static_cast<int>(i_slot.size()));
  }

  S4 out("dgCMatrix");
  out.slot("i") = wrap(i_slot);
  out.slot("p") = wrap(p_slot);
  out.slot("x") = wrap(x_slot);
  out.slot("Dim") = IntegerVector::create(n, n);
  out.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
  out.slot("factors") = List::create();
  return out;
}

// [[Rcpp::export]]
NumericMatrix sctenifold_pcnet_covariance_raw_cpp(
  NumericMatrix x,
  int n_comp = 3,
  int ncv = 0,
  int maxit = 1000,
  double tol = 1e-10,
  int n_threads = 1
) {
  Eigen::MatrixXd coefficients = sctenifold_pcnet_covariance_raw_eigen(
    x,
    n_comp,
    ncv,
    maxit,
    tol,
    n_threads
  );
  NumericMatrix out(coefficients.rows(), coefficients.cols());
  for (int j = 0; j < coefficients.cols(); ++j) {
    for (int i = 0; i < coefficients.rows(); ++i) {
      out(i, j) = coefficients(i, j);
    }
  }
  return out;
}

// [[Rcpp::export]]
S4 sctenifold_pcnet_covariance_sparse_cpp(
  NumericMatrix x,
  int n_comp = 3,
  bool scale_scores = true,
  bool symmetric = false,
  double q = 0.0,
  int ncv = 0,
  int maxit = 1000,
  double tol = 1e-10,
  int n_threads = 1
) {
  Eigen::MatrixXd coefficients = sctenifold_pcnet_covariance_raw_eigen(
    x,
    n_comp,
    ncv,
    maxit,
    tol,
    n_threads
  );
  return sctenifold_dgCMatrix_from_dense_threshold(
    coefficients,
    scale_scores,
    symmetric,
    q
  );
}

// [[Rcpp::export]]
NumericMatrix sctenifold_tensor_decomposition_cpp(
  List x_list,
  List init_u,
  int max_iter = 1000,
  double tol = 1e-5
) {
  const int n_net_int = x_list.size();
  if (n_net_int < 1) {
    stop("x_list must contain at least one network");
  }
  const arma::uword n_net = static_cast<arma::uword>(n_net_int);

  std::vector<arma::mat> tensors;
  tensors.reserve(n_net);
  NumericMatrix first_mat(x_list[0]);
  const arma::uword n_genes = first_mat.nrow();
  if (n_genes != static_cast<arma::uword>(first_mat.ncol())) {
    stop("all networks must be square matrices");
  }

  double tensor_ss = 0.0;
  for (arma::uword l = 0; l < n_net; ++l) {
    NumericMatrix xm(x_list[l]);
    if (static_cast<arma::uword>(xm.nrow()) != n_genes ||
        static_cast<arma::uword>(xm.ncol()) != n_genes) {
      stop("all networks must have the same dimensions");
    }
    arma::mat x(xm.begin(), n_genes, n_genes, false);
    tensors.push_back(x);
    tensor_ss += arma::accu(x % x);
  }
  const double tensor_norm = std::sqrt(tensor_ss);
  if (tensor_norm == 0.0 || !std::isfinite(tensor_norm)) {
    stop("tensor norm must be finite and non-zero");
  }

  arma::mat u1 = as<arma::mat>(init_u[0]);
  arma::mat u2 = as<arma::mat>(init_u[1]);
  arma::mat u3 = as<arma::mat>(init_u[2]);
  arma::mat u4 = as<arma::mat>(init_u[3]);
  if (u1.n_rows != n_genes || u2.n_rows != n_genes ||
      u3.n_rows != 1 || u4.n_rows != n_net ||
      u1.n_cols != u2.n_cols || u1.n_cols != u3.n_cols ||
      u1.n_cols != u4.n_cols) {
    stop("init_u dimensions do not match tensor dimensions");
  }

  const arma::uword k = u1.n_cols;
  arma::rowvec lambdas(k, arma::fill::ones);
  double prev_resid = NA_REAL;
  bool converged = false;
  int curr_iter = 1;

  while (curr_iter < max_iter && !converged) {
    for (int mode = 1; mode <= 4; ++mode) {
      arma::mat v(k, k, arma::fill::ones);
      if (mode != 1) v %= sctenifold_crossprod(u1);
      if (mode != 2) v %= sctenifold_crossprod(u2);
      if (mode != 3) v %= sctenifold_crossprod(u3);
      if (mode != 4) v %= sctenifold_crossprod(u4);

      arma::mat tmp = sctenifold_tensor_mttkrp(
        tensors,
        u1,
        u2,
        u3,
        u4,
        mode
      ) * sctenifold_safe_inv(v);

      for (arma::uword r = 0; r < k; ++r) {
        lambdas[r] = arma::norm(tmp.col(r), 2);
        tmp.col(r) /= lambdas[r];
      }

      if (mode == 1) {
        u1 = tmp;
      } else if (mode == 2) {
        u2 = tmp;
      } else if (mode == 3) {
        u3 = tmp;
      } else {
        u4 = tmp;
      }
    }

    const double resid = sctenifold_tensor_residual(tensors, u1, u2, u3, u4, lambdas);
    if (curr_iter > 1 && std::abs(resid - prev_resid) / tensor_norm < tol) {
      converged = true;
    } else {
      prev_resid = resid;
      ++curr_iter;
    }
  }

  arma::mat out = sctenifold_tensor_average(n_genes, n_net, u1, u2, u3, u4, lambdas);
  NumericMatrix result = wrap(out);
  return result;
}

// [[Rcpp::export]]
NumericMatrix sctenifold_strict_direction_cpp(NumericMatrix x, double lambda = 1.0) {
  const R_xlen_t nr = x.nrow();
  const R_xlen_t nc = x.ncol();
  if (nr != nc) {
    stop("x must be a square matrix");
  }
  if (lambda == 0.0) {
    NumericMatrix out = clone(x);
    out.attr("dimnames") = x.attr("dimnames");
    return out;
  }
  NumericMatrix out(nr, nc);
  for (R_xlen_t j = 0; j < nc; ++j) {
    for (R_xlen_t i = 0; i < nr; ++i) {
      const double value = x(i, j);
      double directed = value;
      if (std::abs(value) < std::abs(x(j, i))) {
        directed = 0.0;
      }
      out(i, j) = ((1.0 - lambda) * value) + (lambda * directed);
    }
  }
  out.attr("dimnames") = x.attr("dimnames");
  return out;
}

// [[Rcpp::export]]
NumericMatrix sctenifold_manifold_matrix_cpp(NumericMatrix x, NumericMatrix y) {
  const R_xlen_t n = x.nrow();
  if (n != x.ncol() || n != y.nrow() || n != y.ncol()) {
    stop("x and y must be square matrices with matching dimensions");
  }

  double sum_x = 0.0;
  double sum_y = 0.0;
  for (R_xlen_t j = 0; j < n; ++j) {
    for (R_xlen_t i = 0; i < n; ++i) {
      sum_x += x(i, j) + 1.0;
      sum_y += y(i, j) + 1.0;
    }
  }
  const double bridge = 0.9 * (sum_x + sum_y) / (2.0 * static_cast<double>(n));
  const R_xlen_t n2 = n * 2;
  NumericMatrix w(n2, n2);

  for (R_xlen_t j = 0; j < n; ++j) {
    for (R_xlen_t i = 0; i < n; ++i) {
      w(i, j) = -(x(i, j) + 1.0);
      w(i + n, j + n) = -(y(i, j) + 1.0);
    }
    w(j + n, j) = -bridge;
    w(j, j + n) = -bridge;
  }

  for (R_xlen_t i = 0; i < n2; ++i) {
    w(i, i) = 0.0;
  }
  for (R_xlen_t j = 0; j < n2; ++j) {
    double col_sum = 0.0;
    for (R_xlen_t i = 0; i < n2; ++i) {
      col_sum += w(i, j);
    }
    w(j, j) = -col_sum;
  }

  return w;
}

// [[Rcpp::export]]
NumericVector sctenifold_pair_distances_cpp(NumericMatrix aligned) {
  const R_xlen_t nr = aligned.nrow();
  const R_xlen_t nc = aligned.ncol();
  if (nr % 2 != 0) {
    stop("aligned must contain X genes followed by Y genes");
  }
  const R_xlen_t n_genes = nr / 2;
  NumericVector distances(n_genes);
  for (R_xlen_t i = 0; i < n_genes; ++i) {
    double ss = 0.0;
    for (R_xlen_t j = 0; j < nc; ++j) {
      const double diff = aligned(i, j) - aligned(i + n_genes, j);
      ss += diff * diff;
    }
    distances[i] = std::sqrt(ss);
  }
  return distances;
}
