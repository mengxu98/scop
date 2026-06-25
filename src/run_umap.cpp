#include <Rcpp.h>
#include <cstdint>
#include <cstdlib>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

extern "C" {
#include "knn.h"
#include "umap.h"
}

static int loaded_pid = -1;

static int process_id() {
#ifdef _WIN32
  return static_cast<int>(_getpid());
#else
  return static_cast<int>(getpid());
#endif
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((constructor))
#endif
static void remember_loader_process() {
  loaded_pid = process_id();
}

extern "C" int loader_process_is_current(void) {
  if (loaded_pid < 0) {
    remember_loader_process();
  }
  return process_id() == loaded_pid ? 1 : 0;
}

static std::vector<float> numeric_matrix_to_float(Rcpp::NumericMatrix x) {
  const int rows = x.nrow();
  const int cols = x.ncol();
  std::vector<float> out(static_cast<size_t>(rows) * cols);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      out[static_cast<size_t>(row) * cols + col] =
        static_cast<float>(x(row, col));
    }
  }
  return out;
}

static Rcpp::NumericMatrix float_layout_to_matrix(const std::vector<float>& x,
                                                  int rows,
                                                  int cols) {
  Rcpp::NumericMatrix out(rows, cols);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      out(row, col) = x[static_cast<size_t>(row) * cols + col];
    }
  }
  return out;
}

Rcpp::List runumap_knn(Rcpp::NumericMatrix X,
                       int k,
                       int n_trees,
                       int n_iters,
                       int leaf_size,
                       double seed,
                       int cores) {
  const int rows = X.nrow();
  const int dims = X.ncol();
  if (rows <= 1 || dims <= 0 || k <= 0 || k >= rows) {
    Rcpp::stop("invalid kNN dimensions");
  }

  std::vector<float> input = numeric_matrix_to_float(X);
  std::vector<int32_t> indices(static_cast<size_t>(rows) * k);
  std::vector<float> distances(static_cast<size_t>(rows) * k);
  const int code = knn_descent_f32(
    rows, dims, k, input.data(), n_trees, n_iters, leaf_size,
    static_cast<uint64_t>(seed), indices.data(), distances.data(), cores);
  if (code != 0) {
    Rcpp::stop("kNN failed with code %d", code);
  }

  Rcpp::IntegerMatrix idx(rows, k);
  Rcpp::NumericMatrix dist(rows, k);
  for (int row = 0; row < rows; ++row) {
    for (int rank = 0; rank < k; ++rank) {
      const size_t offset = static_cast<size_t>(row) * k + rank;
      idx(row, rank) = indices[offset] + 1;
      dist(row, rank) = distances[offset];
    }
  }
  return Rcpp::List::create(
    Rcpp::_["idx"] = idx,
    Rcpp::_["dist"] = dist,
    Rcpp::_["rc"] = code
  );
}

Rcpp::NumericMatrix runumap_layout(Rcpp::NumericMatrix init,
                                  Rcpp::IntegerVector head,
                                  Rcpp::IntegerVector tail,
                                  Rcpp::NumericVector epochs_per_sample,
                                  int n_epochs,
                                  double a,
                                  double b,
                                  double gamma,
                                  double initial_alpha,
                                  double negative_sample_rate,
                                  double seed,
                                  int cores) {
  const int rows = init.nrow();
  const int dims = init.ncol();
  const R_xlen_t edges = head.size();
  if (tail.size() != edges || epochs_per_sample.size() != edges) {
    Rcpp::stop("head, tail, and epochs_per_sample must have equal length");
  }
  if (rows <= 0 || dims <= 0 || edges <= 0 || n_epochs <= 0) {
    Rcpp::stop("invalid UMAP layout dimensions");
  }

  std::vector<float> embedding = numeric_matrix_to_float(init);
  std::vector<int32_t> from(static_cast<size_t>(edges));
  std::vector<int32_t> to(static_cast<size_t>(edges));
  std::vector<float> sample_every(static_cast<size_t>(edges));
  for (R_xlen_t edge = 0; edge < edges; ++edge) {
    from[static_cast<size_t>(edge)] = static_cast<int32_t>(head[edge]);
    to[static_cast<size_t>(edge)] = static_cast<int32_t>(tail[edge]);
    sample_every[static_cast<size_t>(edge)] =
      static_cast<float>(epochs_per_sample[edge]);
  }

  const int code = umap_optimize_layout_euclidean_parallel_f32(
    embedding.data(), rows, dims, from.data(), to.data(),
    static_cast<int64_t>(edges), sample_every.data(), n_epochs,
    static_cast<float>(a), static_cast<float>(b), static_cast<float>(gamma),
    static_cast<float>(initial_alpha),
    static_cast<float>(negative_sample_rate),
    static_cast<uint64_t>(seed), cores);

  Rcpp::NumericMatrix out = float_layout_to_matrix(embedding, rows, dims);
  out.attr("layout_return_code") = code;
  return out;
}
