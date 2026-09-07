#include <Rcpp.h>
#include <thisutils/log_message.h>
#include "numeric_utils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>
#include "dynload.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static int core_count(int requested, int jobs) {
  int n = requested > 0 ? requested :
    static_cast<int>(std::thread::hardware_concurrency());
  if (n < 1) n = 1;
  if (jobs > 0 && n > jobs) n = jobs;
  return n;
}

static float squared_distance(const float* a, const float* b, int width) {
  float total = 0.0f;
  for (int pos = 0; pos < width; ++pos) {
    const float delta = a[pos] - b[pos];
    total += delta * delta;
  }
  return total;
}

struct Candidate {
  float distance;
  int index;
};

static bool worse_candidate(const Candidate& a, const Candidate& b) {
  if (a.distance == b.distance) return a.index < b.index;
  return a.distance < b.distance;
}

typedef void (*matrix_multiply_f32)(int, int, int, int, int, int, float,
                                    const float*, int, const float*, int,
                                    float, float*, int);

static matrix_multiply_f32 runtime_sgemm() {
  static matrix_multiply_f32 fn = NULL;
  static bool loaded = false;
  if (!loaded) {
    loaded = true;
    // A cblas already linked into the process (Accelerate on macOS, the
    // OpenBLAS/MKL/BLIS that R itself links on Linux) is the reliable hit.
    void* self = scop_dlopen(NULL);
    if (self != NULL) {
      fn = reinterpret_cast<matrix_multiply_f32>(scop_dlsym(self, "cblas_sgemm"));
    }
    const char* candidates[] = {
      "/System/Library/Frameworks/Accelerate.framework/Accelerate",
      "libopenblas.so.0", "libopenblas.so", "libmkl_rt.so",
      "libblis.so.4", "libblis.so.3", "libblis.so",
      "libblas.so.3", "libblas.so",
      "Rblas.dll", "libopenblas.dll", "openblas.dll", "mkl_rt.dll", NULL
    };
    for (int i = 0; fn == NULL && candidates[i] != NULL; ++i) {
      void* handle = scop_dlopen(candidates[i]);
      if (handle != NULL) {
        fn = reinterpret_cast<matrix_multiply_f32>(scop_dlsym(handle, "cblas_sgemm"));
      }
    }
  }
  return fn;
}

static void blas_set_num_threads(int n) {
  typedef void (*set_fn)(int);
  static set_fn fn = NULL;
  static bool loaded = false;
  if (!loaded) {
    loaded = true;
    void* self = scop_dlopen(NULL);
    if (self != NULL) {
      fn = reinterpret_cast<set_fn>(scop_dlsym(self, "openblas_set_num_threads"));
      if (fn == NULL) {
        fn = reinterpret_cast<set_fn>(scop_dlsym(self, "MKL_Set_Num_Threads"));
      }
    }
  }
  if (fn != NULL && n > 0) {
    fn(n);
  }
}

static int blas_get_num_threads() {
  typedef int (*get_fn)();
  static get_fn fn = NULL;
  static bool loaded = false;
  if (!loaded) {
    loaded = true;
    void* self = scop_dlopen(NULL);
    if (self != NULL) {
      fn = reinterpret_cast<get_fn>(scop_dlsym(self, "openblas_get_num_threads"));
    }
  }
  return fn == NULL ? 0 : fn();
}

static void sorted_from_heap(std::vector<Candidate>& heap,
                             int* idx_out,
                             int query,
                             const std::vector<float>& norms,
                             float* dist_out,
                             int rows) {
  std::sort_heap(heap.begin(), heap.end(), worse_candidate);
  for (int rank = 0; rank < static_cast<int>(heap.size()); ++rank) {
    const int ref = heap[rank].index;
    float value = norms[query] + heap[rank].distance;
    if (value < 0.0f) value = 0.0f;
    idx_out[query + rank * rows] = ref + 1;
    dist_out[query + rank * rows] = value;
  }
}

static bool exact_blas(const std::vector<float>& data,
                       int rows,
                       int cols,
                       int neighbors,
                       int cores,
                       int* idx_out,
                       float* dist_out) {
  matrix_multiply_f32 sgemm = runtime_sgemm();
  if (sgemm == NULL) return false;

  const int workers = core_count(cores, rows);
  std::vector<float> norms(static_cast<size_t>(rows));
  for (int row = 0; row < rows; ++row) {
    const float* current = data.data() + static_cast<size_t>(row) * cols;
    float total = 0.0f;
    for (int col = 0; col < cols; ++col) {
      total += current[col] * current[col];
    }
    norms[row] = total;
  }

  const long tile_mb = 32L;
  long tile = (tile_mb << 20) /
    (static_cast<long>(rows) * static_cast<long>(sizeof(float)));
  int tile_rows = tile < 64L ? 64 : static_cast<int>(tile);
  if (tile_rows > rows) tile_rows = rows;

  std::vector<float> scores(static_cast<size_t>(tile_rows) * rows);
  enum { row_major = 101, no_trans = 111, trans = 112 };

  // Skinny GEMMs (kNN on 20–50 PCs) lose time with 48 BLAS threads. Cap the
  // GEMM pool; the heap search still uses `workers` OpenMP threads after each
  // tile, when BLAS is idle.
  const int gemm_threads = std::max(1, std::min(workers, std::max(1, (cols + 7) / 8)));
  const int previous_blas = blas_get_num_threads();
  blas_set_num_threads(gemm_threads);

  for (int first = 0; first < rows; first += tile_rows) {
    const int count = std::min(tile_rows, rows - first);
    sgemm(
      row_major,
      no_trans,
      trans,
      count,
      rows,
      cols,
      1.0f,
      data.data() + static_cast<size_t>(first) * cols,
      cols,
      data.data(),
      cols,
      0.0f,
      scores.data(),
      rows
    );

#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
    {
      std::vector<Candidate> heap(static_cast<size_t>(neighbors));
#pragma omp for schedule(static)
      for (int local = 0; local < count; ++local) {
        const int query = first + local;
        const float* score_row = scores.data() + static_cast<size_t>(local) * rows;
        for (int ref = 0; ref < neighbors; ++ref) {
          heap[ref] = Candidate{
            norms[ref] - 2.0f * score_row[ref],
            ref
          };
        }
        std::make_heap(heap.begin(), heap.end(), worse_candidate);
        for (int ref = neighbors; ref < rows; ++ref) {
          float value = norms[ref] - 2.0f * score_row[ref];
          if (value < heap.front().distance) {
            std::pop_heap(heap.begin(), heap.end(), worse_candidate);
            heap.back() = Candidate{value, ref};
            std::push_heap(heap.begin(), heap.end(), worse_candidate);
          }
        }
        sorted_from_heap(heap, idx_out, query, norms, dist_out, rows);
      }
    }
#else
    std::vector<std::thread> workers_pool;
    workers_pool.reserve(static_cast<size_t>(workers));
    for (int worker = 0; worker < workers; ++worker) {
      workers_pool.emplace_back([&, worker]() {
        std::vector<Candidate> heap(static_cast<size_t>(neighbors));
        for (int local = worker; local < count; local += workers) {
          const int query = first + local;
          const float* score_row = scores.data() + static_cast<size_t>(local) * rows;
          for (int ref = 0; ref < neighbors; ++ref) {
            heap[ref] = Candidate{
              norms[ref] - 2.0f * score_row[ref],
              ref
            };
          }
          std::make_heap(heap.begin(), heap.end(), worse_candidate);
          for (int ref = neighbors; ref < rows; ++ref) {
            float value = norms[ref] - 2.0f * score_row[ref];
            if (value < heap.front().distance) {
              std::pop_heap(heap.begin(), heap.end(), worse_candidate);
              heap.back() = Candidate{value, ref};
              std::push_heap(heap.begin(), heap.end(), worse_candidate);
            }
          }
          sorted_from_heap(heap, idx_out, query, norms, dist_out, rows);
        }
      });
    }
    for (std::thread& worker : workers_pool) worker.join();
#endif
  }
  if (previous_blas > 0) {
    blas_set_num_threads(previous_blas);
  }
  return true;
}

static void insert_candidate(std::vector<Candidate>& heap, Candidate candidate) {
  if (candidate.distance < heap.front().distance ||
      (candidate.distance == heap.front().distance &&
       candidate.index < heap.front().index)) {
    std::pop_heap(heap.begin(), heap.end(), worse_candidate);
    heap.back() = candidate;
    std::push_heap(heap.begin(), heap.end(), worse_candidate);
  }
}

static void exact_worker(const std::vector<float>& data,
                         int rows,
                         int cols,
                         int neighbors,
                         int begin,
                         int end,
                         int* output) {
  std::vector<Candidate> heap(static_cast<size_t>(neighbors));
  for (int query = begin; query < end; ++query) {
    const float* query_row = data.data() + static_cast<size_t>(query) * cols;
    for (int ref = 0; ref < neighbors; ++ref) {
      heap[ref] = Candidate{
        squared_distance(query_row, data.data() + static_cast<size_t>(ref) * cols, cols),
        ref
      };
    }
    std::make_heap(heap.begin(), heap.end(), worse_candidate);
    for (int ref = neighbors; ref < rows; ++ref) {
      insert_candidate(heap, Candidate{
        squared_distance(query_row, data.data() + static_cast<size_t>(ref) * cols, cols),
        ref
      });
    }
    std::sort_heap(heap.begin(), heap.end(), worse_candidate);
    for (int rank = 0; rank < neighbors; ++rank) {
      output[query + rank * rows] = heap[rank].index + 1;
    }
  }
}

static void cross_knn_worker(const std::vector<float>& reference,
                             const std::vector<float>& query,
                             const std::vector<float>& reference_norms,
                             const std::vector<float>& query_norms,
                             int reference_rows,
                             int query_rows,
                             int cols,
                             int neighbors,
                             bool cosine,
                             int begin,
                             int end,
                             int* idx_output,
                             double* distance_output) {
  std::vector<Candidate> heap(static_cast<size_t>(neighbors));
  for (int query_row = begin; query_row < end; ++query_row) {
    const float* current = query.data() +
      static_cast<size_t>(query_row) * cols;
    for (int reference_row = 0; reference_row < neighbors; ++reference_row) {
      const float* candidate = reference.data() +
        static_cast<size_t>(reference_row) * cols;
      float distance = 0.0f;
      if (cosine) {
        for (int col = 0; col < cols; ++col) {
          distance += current[col] * candidate[col];
        }
        distance = 1.0f - distance /
          (query_norms[query_row] * reference_norms[reference_row]);
      } else {
        distance = squared_distance(current, candidate, cols);
      }
      heap[reference_row] = Candidate{distance, reference_row};
    }
    std::make_heap(heap.begin(), heap.end(), worse_candidate);
    for (int reference_row = neighbors; reference_row < reference_rows;
         ++reference_row) {
      const float* candidate = reference.data() +
        static_cast<size_t>(reference_row) * cols;
      float distance = 0.0f;
      if (cosine) {
        for (int col = 0; col < cols; ++col) {
          distance += current[col] * candidate[col];
        }
        distance = 1.0f - distance /
          (query_norms[query_row] * reference_norms[reference_row]);
      } else {
        distance = squared_distance(current, candidate, cols);
      }
      insert_candidate(heap, Candidate{distance, reference_row});
    }
    std::sort_heap(heap.begin(), heap.end(), worse_candidate);
    for (int rank = 0; rank < neighbors; ++rank) {
      float distance = heap[rank].distance;
      if (!cosine) distance = std::sqrt(std::max(0.0f, distance));
      idx_output[query_row + rank * query_rows] = heap[rank].index + 1;
      distance_output[query_row + rank * query_rows] = distance;
    }
  }
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix exact_knn_f32(Rcpp::NumericMatrix data,
                                 int k,
                                 int cores) {
  const int rows = data.nrow();
  const int dims = data.ncol();
  if (rows <= 0 || dims <= 0 || k <= 0 || k > rows) {
    thisutils::log_message("invalid exact kNN dimensions", "error");
  }
  std::vector<float> packed = matrix_as_row_float(data);
  Rcpp::IntegerMatrix neighbors(rows, k);
  int* out = INTEGER(neighbors);
  std::vector<float> distances(static_cast<size_t>(rows) * k);
  if (exact_blas(packed, rows, dims, k, cores, out, distances.data())) {
    return neighbors;
  }
  const int workers = core_count(cores, rows);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<size_t>(workers));
  for (int worker = 0; worker < workers; ++worker) {
    const int begin = static_cast<int>(
      static_cast<int64_t>(worker) * rows / workers);
    const int end = static_cast<int>(
      static_cast<int64_t>(worker + 1) * rows / workers);
    pool.emplace_back(exact_worker, std::cref(packed), rows, dims, k,
                      begin, end, out);
  }
  for (std::thread& worker : pool) worker.join();
  return neighbors;
}

// [[Rcpp::export]]
Rcpp::List cross_knn_f32(Rcpp::NumericMatrix reference,
                         Rcpp::NumericMatrix query,
                         int k,
                         std::string metric,
                         int cores) {
  const int reference_rows = reference.nrow();
  const int query_rows = query.nrow();
  const int dims = reference.ncol();
  const bool cosine = metric == "cosine";
  if (reference_rows <= 0 || query_rows <= 0 || dims <= 0 ||
      query.ncol() != dims || k <= 0 || k > reference_rows ||
      (!cosine && metric != "euclidean")) {
    thisutils::log_message("invalid cross kNN dimensions or metric", "error");
  }

  std::vector<float> reference_packed = matrix_as_row_float(reference);
  std::vector<float> query_packed = matrix_as_row_float(query);
  for (size_t i = 0; i < reference_packed.size(); ++i) {
    if (!std::isfinite(reference_packed[i])) {
      thisutils::log_message("cross kNN requires finite reference values", "error");
    }
  }
  for (size_t i = 0; i < query_packed.size(); ++i) {
    if (!std::isfinite(query_packed[i])) {
      thisutils::log_message("cross kNN requires finite query values", "error");
    }
  }

  std::vector<float> reference_norms(static_cast<size_t>(reference_rows), 1.0f);
  std::vector<float> query_norms(static_cast<size_t>(query_rows), 1.0f);
  if (cosine) {
    for (int row = 0; row < reference_rows; ++row) {
      const float* current = reference_packed.data() +
        static_cast<size_t>(row) * dims;
      float norm_sq = 0.0f;
      for (int col = 0; col < dims; ++col) norm_sq += current[col] * current[col];
      if (!(norm_sq > 0.0f) || !std::isfinite(norm_sq)) {
        thisutils::log_message("cross cosine kNN requires non-zero finite reference rows", "error");
      }
      reference_norms[row] = std::sqrt(norm_sq);
    }
    for (int row = 0; row < query_rows; ++row) {
      const float* current = query_packed.data() +
        static_cast<size_t>(row) * dims;
      float norm_sq = 0.0f;
      for (int col = 0; col < dims; ++col) norm_sq += current[col] * current[col];
      if (!(norm_sq > 0.0f) || !std::isfinite(norm_sq)) {
        thisutils::log_message("cross cosine kNN requires non-zero finite query rows", "error");
      }
      query_norms[row] = std::sqrt(norm_sq);
    }
  }

  Rcpp::IntegerMatrix idx(query_rows, k);
  Rcpp::NumericMatrix distance(query_rows, k);
  const int workers = core_count(cores, query_rows);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<size_t>(workers));
  for (int worker = 0; worker < workers; ++worker) {
    const int begin = static_cast<int>(
      static_cast<int64_t>(worker) * query_rows / workers);
    const int end = static_cast<int>(
      static_cast<int64_t>(worker + 1) * query_rows / workers);
    pool.emplace_back(
      cross_knn_worker,
      std::cref(reference_packed),
      std::cref(query_packed),
      std::cref(reference_norms),
      std::cref(query_norms),
      reference_rows,
      query_rows,
      dims,
      k,
      cosine,
      begin,
      end,
      INTEGER(idx),
      REAL(distance)
    );
  }
  for (std::thread& worker : pool) worker.join();
  return Rcpp::List::create(
    Rcpp::Named("idx") = idx,
    Rcpp::Named("distance") = distance
  );
}
