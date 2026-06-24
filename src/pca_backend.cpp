#include <RcppArmadillo.h>
#ifdef __APPLE__
#include <dlfcn.h>
#endif

static int component_count(int requested, int features, int cells) {
  int limit = std::min(features, cells - 1);
  if (limit < 1) return 0;
  return std::max(1, std::min(requested, limit));
}

#ifdef __APPLE__
namespace {
enum {
  scop_cblas_col_major = 102,
  scop_cblas_no_trans = 111,
  scop_cblas_trans = 112,
  scop_cblas_lower = 122
};

using scop_dsyrk_ptr = void (*)(int, int, int, int, int, double,
                                const double*, int, double, double*, int);
using scop_dgemm_ptr = void (*)(int, int, int, int, int, int, double,
                                const double*, int, const double*, int,
                                double, double*, int);
using scop_dsyevr_ptr = void (*)(const char*, const char*, const char*, const int*,
                                 double*, const int*, const double*, const double*,
                                 const int*, const int*, const double*, int*,
                                 double*, double*, const int*, int*, double*,
                                 const int*, int*, const int*, int*);

struct PcaBlas {
  void* handle = nullptr;
  scop_dsyrk_ptr dsyrk = nullptr;
  scop_dgemm_ptr dgemm = nullptr;
  scop_dsyevr_ptr dsyevr = nullptr;
};

static PcaBlas scop_pca_blas;
static bool scop_pca_blas_checked = false;

static PcaBlas* resolve_pca_blas() {
  if (scop_pca_blas_checked) {
    return (scop_pca_blas.dsyrk && scop_pca_blas.dgemm) ? &scop_pca_blas : nullptr;
  }
  scop_pca_blas_checked = true;
  void* handle = dlopen(
    "/System/Library/Frameworks/Accelerate.framework/Accelerate",
    RTLD_LAZY | RTLD_LOCAL
  );
  if (handle == nullptr) {
    return nullptr;
  }
  scop_pca_blas.handle = handle;
  scop_pca_blas.dsyrk = reinterpret_cast<scop_dsyrk_ptr>(dlsym(handle, "cblas_dsyrk"));
  scop_pca_blas.dgemm = reinterpret_cast<scop_dgemm_ptr>(dlsym(handle, "cblas_dgemm"));
  scop_pca_blas.dsyevr = reinterpret_cast<scop_dsyevr_ptr>(dlsym(handle, "dsyevr_"));
  return (scop_pca_blas.dsyrk && scop_pca_blas.dgemm) ? &scop_pca_blas : nullptr;
}

static bool fast_gram(const arma::mat& X, arma::mat& gram) {
  PcaBlas* blas = resolve_pca_blas();
  if (blas == nullptr) {
    return false;
  }
  const int features = X.n_rows;
  const int cells = X.n_cols;

  gram.zeros(features, features);
  blas->dsyrk(
    scop_cblas_col_major,
    scop_cblas_lower,
    scop_cblas_no_trans,
    features,
    cells,
    1.0,
    X.memptr(),
    features,
    0.0,
    gram.memptr(),
    features
  );
  for (int col = 0; col < features; ++col) {
    for (int row = col + 1; row < features; ++row) {
      gram(col, row) = gram(row, col);
    }
  }
  return true;
}

static bool fast_embeddings(const arma::mat& X,
                            const arma::mat& loadings,
                            arma::mat& embeddings) {
  PcaBlas* blas = resolve_pca_blas();
  if (blas == nullptr) {
    return false;
  }
  const int features = X.n_rows;
  const int cells = X.n_cols;
  const int keep = loadings.n_cols;
  embeddings.set_size(cells, keep);
  blas->dgemm(
    scop_cblas_col_major,
    scop_cblas_trans,
    scop_cblas_no_trans,
    cells,
    keep,
    features,
    1.0,
    X.memptr(),
    features,
    loadings.memptr(),
    features,
    0.0,
    embeddings.memptr(),
    cells
  );
  return true;
}

static bool fast_top_eigen(arma::mat& gram,
                           int keep,
                           arma::vec& eigvals,
                           arma::mat& loadings) {
  PcaBlas* blas = resolve_pca_blas();
  if (blas == nullptr || blas->dsyevr == nullptr) {
    return false;
  }
  const int features = gram.n_rows;
  if (keep < 1 || keep > features) {
    return false;
  }

  const int il = features - keep + 1;
  const int iu = features;
  int m = 0;
  int info = 0;
  const int ldz = features;
  const double vl = 0.0;
  const double vu = 0.0;
  const double abstol = -1.0;
  std::vector<double> values(features);
  std::vector<double> vectors(static_cast<size_t>(features) * keep);
  std::vector<int> support(static_cast<size_t>(2) * keep);

  double work_query = 0.0;
  int iwork_query = 0;
  int lwork = -1;
  int liwork = -1;
  blas->dsyevr(
    "V", "I", "L",
    &features,
    gram.memptr(),
    &features,
    &vl,
    &vu,
    &il,
    &iu,
    &abstol,
    &m,
    values.data(),
    vectors.data(),
    &ldz,
    support.data(),
    &work_query,
    &lwork,
    &iwork_query,
    &liwork,
    &info
  );
  if (info != 0 || !std::isfinite(work_query) || iwork_query < 1) {
    return false;
  }
  lwork = static_cast<int>(work_query);
  liwork = iwork_query;
  std::vector<double> work(static_cast<size_t>(lwork));
  std::vector<int> iwork(static_cast<size_t>(liwork));
  blas->dsyevr(
    "V", "I", "L",
    &features,
    gram.memptr(),
    &features,
    &vl,
    &vu,
    &il,
    &iu,
    &abstol,
    &m,
    values.data(),
    vectors.data(),
    &ldz,
    support.data(),
    work.data(),
    &lwork,
    iwork.data(),
    &liwork,
    &info
  );
  if (info != 0 || m < keep) {
    return false;
  }

  eigvals.set_size(keep);
  loadings.set_size(features, keep);
  for (int j = 0; j < keep; ++j) {
    const int source = keep - 1 - j;
    eigvals[j] = values[static_cast<size_t>(source)];
    std::copy(
      vectors.begin() + static_cast<size_t>(source) * features,
      vectors.begin() + static_cast<size_t>(source + 1) * features,
      loadings.colptr(j)
    );
  }
  return true;
}
}
#endif

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

  arma::mat gram;
#ifdef __APPLE__
  if (!fast_gram(X, gram)) {
    gram = X * X.t();
  }
#else
  gram = X * X.t();
#endif
  arma::vec values;
  arma::mat vectors;
  arma::mat loadings(features, keep);
  arma::vec eigvals(keep);
#ifdef __APPLE__
  arma::mat gram_for_eigen = gram;
  if (!fast_top_eigen(gram_for_eigen, keep, eigvals, loadings)) {
    if (!arma::eig_sym(values, vectors, gram)) {
      Rcpp::stop("pca_backend_run: eigensolver failed");
    }
    for (int j = 0; j < keep; ++j) {
      const arma::uword source = static_cast<arma::uword>(values.n_elem - 1 - j);
      eigvals[j] = values[source];
      loadings.col(j) = vectors.col(source);
    }
  }
#else
  if (!arma::eig_sym(values, vectors, gram)) {
    Rcpp::stop("pca_backend_run: eigensolver failed");
  }
  for (int j = 0; j < keep; ++j) {
    const arma::uword source = static_cast<arma::uword>(values.n_elem - 1 - j);
    eigvals[j] = values[source];
    loadings.col(j) = vectors.col(source);
  }
#endif

  arma::mat embeddings;
#ifdef __APPLE__
  if (!fast_embeddings(X, loadings, embeddings)) {
    embeddings = X.t() * loadings;
  }
#else
  embeddings = X.t() * loadings;
#endif
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
