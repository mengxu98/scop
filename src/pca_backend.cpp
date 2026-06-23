#include <RcppArmadillo.h>

static int component_count(int requested, int features, int cells) {
  int limit = std::min(features, cells - 1);
  if (limit < 1) return 0;
  return std::max(1, std::min(requested, limit));
}

// [[Rcpp::export]]
Rcpp::List pca_backend_run(const arma::mat& X,
                           int npcs,
                           bool weight_by_var = true) {
  const int features = X.n_rows;
  const int cells = X.n_cols;
  const int keep = component_count(npcs, features, cells);
  if (keep < 1) {
    Rcpp::stop("pca_backend_run: matrix is too small");
  }

  arma::mat gram = X * X.t();
  arma::vec values;
  arma::mat vectors;
  if (!arma::eig_sym(values, vectors, gram)) {
    Rcpp::stop("pca_backend_run: eigensolver failed");
  }

  arma::mat loadings(features, keep);
  arma::vec eigvals(keep);
  for (int j = 0; j < keep; ++j) {
    const arma::uword source = static_cast<arma::uword>(values.n_elem - 1 - j);
    eigvals[j] = values[source];
    loadings.col(j) = vectors.col(source);
  }

  arma::mat embeddings = X.t() * loadings;
  arma::vec singular = arma::sqrt(arma::clamp(eigvals, 0.0, arma::datum::inf));
  arma::vec sdev = singular / std::sqrt(static_cast<double>(std::max(1, cells - 1)));
  if (!weight_by_var) {
    for (int j = 0; j < keep; ++j) {
      if (singular[j] > 0.0) embeddings.col(j) /= singular[j];
    }
  }

  return Rcpp::List::create(
    Rcpp::_["embeddings"] = embeddings,
    Rcpp::_["loadings"] = loadings,
    Rcpp::_["sdev"] = sdev,
    Rcpp::_["eigvals"] = eigvals
  );
}

arma::mat cca_crossprod_matrix(const arma::mat& X1,
                               const arma::mat& X2) {
  return X1.t() * X2;
}

arma::mat matrix_product(const arma::mat& A,
                         const arma::mat& B) {
  return A * B;
}
