#include <Rcpp.h>
#include <thisutils/log_message.h>

using namespace Rcpp;

// [[Rcpp::export]]
List knn_vote_labels_cpp(const IntegerMatrix& labels, int n_levels) {
  const int n = labels.nrow();
  const int k = labels.ncol();
  if (n_levels < 1) {
    thisutils::log_message("n_levels must be positive", "error");
  }

  NumericMatrix probability(n, n_levels);
  IntegerVector best(n, NA_INTEGER);
  for (int i = 0; i < n; ++i) {
    int valid = 0;
    for (int j = 0; j < k; ++j) {
      const int label = labels(i, j);
      if (label == NA_INTEGER) {
        continue;
      }
      if (label < 1 || label > n_levels) {
        thisutils::log_message("labels contains an out-of-range level code", "error");
      }
      probability(i, label - 1) += 1.0;
      ++valid;
    }
    if (valid == 0) {
      for (int level = 0; level < n_levels; ++level) {
        probability(i, level) = NA_REAL;
      }
      continue;
    }

    int best_level = 0;
    double best_count = probability(i, 0);
    for (int level = 0; level < n_levels; ++level) {
      probability(i, level) /= static_cast<double>(valid);
      if (level > 0 && probability(i, level) > best_count / valid) {
        best_count = probability(i, level) * valid;
        best_level = level;
      }
    }
    best[i] = best_level + 1;
  }

  return List::create(
    _["probability"] = probability,
    _["best"] = best
  );
}
