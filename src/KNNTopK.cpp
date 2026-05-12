#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

using namespace Rcpp;

struct KnnTopKEntry {
  double score;
  int index;
};

static bool knn_topk_less(const KnnTopKEntry& a, const KnnTopKEntry& b) {
  if (a.score < b.score) {
    return true;
  }
  if (a.score > b.score) {
    return false;
  }
  return a.index < b.index;
}

struct KnnTopKWorker {
  const double* reference;
  const double* query;
  const double* ref_norm;
  int* idx;
  double* dist;
  const int n_ref;
  const int n_query;
  const int dims;
  const int k;
  const int metric_id;
  const bool exclude_self;

  KnnTopKWorker(
    const double* reference,
    const double* query,
    const double* ref_norm,
    int* idx,
    double* dist,
    int n_ref,
    int n_query,
    int dims,
    int k,
    int metric_id,
    bool exclude_self
  ) : reference(reference),
      query(query),
      ref_norm(ref_norm),
      idx(idx),
      dist(dist),
      n_ref(n_ref),
      n_query(n_query),
      dims(dims),
      k(k),
      metric_id(metric_id),
      exclude_self(exclude_self) {}

  inline double reference_at(int row, int col) const {
    return reference[row + n_ref * col];
  }

  inline double query_at(int row, int col) const {
    return query[row + n_query * col];
  }

  void run(std::size_t begin, std::size_t end) const {
    std::vector<KnnTopKEntry> entries;
    entries.reserve(n_ref);

    for (std::size_t qi = begin; qi < end; ++qi) {
      entries.clear();
      double q_norm = 0.0;
      if (metric_id == 2) {
        for (int d = 0; d < dims; ++d) {
          const double qv = query_at(qi, d);
          q_norm += qv * qv;
        }
        q_norm = std::sqrt(q_norm);
      }

      for (int ri = 0; ri < n_ref; ++ri) {
        if (exclude_self && static_cast<int>(qi) == ri) {
          continue;
        }

        double score = 0.0;
        if (metric_id == 1) {
          for (int d = 0; d < dims; ++d) {
            const double diff = query_at(qi, d) - reference_at(ri, d);
            score += diff * diff;
          }
        } else {
          double dot = 0.0;
          for (int d = 0; d < dims; ++d) {
            dot += query_at(qi, d) * reference_at(ri, d);
          }
          const double denom = q_norm * ref_norm[ri];
          const double sim = denom > 0.0 ? dot / denom : 0.0;
          score = 1.0 - sim;
        }

        KnnTopKEntry entry;
        entry.score = score;
        entry.index = ri + 1;
        entries.push_back(entry);
      }

      const int take = std::min(k, static_cast<int>(entries.size()));
      if (take > 0) {
        std::partial_sort(
          entries.begin(),
          entries.begin() + take,
          entries.end(),
          knn_topk_less
        );
      }

      for (int j = 0; j < take; ++j) {
        idx[qi + n_query * j] = entries[j].index;
        dist[qi + n_query * j] = metric_id == 1 ? std::sqrt(entries[j].score) : entries[j].score;
      }
    }
  }
};

// [[Rcpp::export]]
List knn_topk_cpp(
  NumericMatrix reference,
  NumericMatrix query,
  int k,
  std::string metric = "euclidean",
  bool exclude_self = false,
  int n_threads = 0
) {
  if (reference.ncol() != query.ncol()) {
    stop("reference and query must have the same number of columns");
  }
  if (k < 1) {
    stop("k must be positive");
  }

  const int n_ref = reference.nrow();
  const int n_query = query.nrow();
  const int dims = reference.ncol();
  const int metric_id = metric == "cosine" ? 2 : 1;
  const int take = std::min(k, n_ref - static_cast<int>(exclude_self));
  if (take < 1 || n_ref < 1 || n_query < 1) {
    return List::create(
      _["idx"] = IntegerMatrix(n_query, 0),
      _["dist"] = NumericMatrix(n_query, 0)
    );
  }

  NumericVector ref_norm(n_ref);
  if (metric_id == 2) {
    for (int ri = 0; ri < n_ref; ++ri) {
      double norm = 0.0;
      for (int d = 0; d < dims; ++d) {
        const double value = reference(ri, d);
        norm += value * value;
      }
      ref_norm[ri] = std::sqrt(norm);
    }
  }

  IntegerMatrix idx(n_query, take);
  NumericMatrix dist(n_query, take);
  std::fill(idx.begin(), idx.end(), NA_INTEGER);
  std::fill(dist.begin(), dist.end(), NA_REAL);

  KnnTopKWorker worker(
    REAL(reference),
    REAL(query),
    REAL(ref_norm),
    INTEGER(idx),
    REAL(dist),
    n_ref,
    n_query,
    dims,
    take,
    metric_id,
    exclude_self
  );

  unsigned int thread_count = n_threads > 0 ?
    static_cast<unsigned int>(n_threads) :
    std::thread::hardware_concurrency();
  if (thread_count == 0) {
    thread_count = 1;
  }
  thread_count = std::min<unsigned int>(thread_count, static_cast<unsigned int>(n_query));

  if (thread_count <= 1) {
    worker.run(0, n_query);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    const std::size_t block = (n_query + thread_count - 1) / thread_count;
    for (unsigned int t = 0; t < thread_count; ++t) {
      const std::size_t begin = t * block;
      const std::size_t end = std::min<std::size_t>(begin + block, n_query);
      if (begin >= end) {
        break;
      }
      threads.emplace_back([&worker, begin, end]() {
        worker.run(begin, end);
      });
    }
    for (std::thread& thread : threads) {
      thread.join();
    }
  }

  return List::create(
    _["idx"] = idx,
    _["dist"] = dist
  );
}
