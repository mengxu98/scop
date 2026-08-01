#include <Rcpp.h>
#include <R_ext/Random.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

inline void score_row(
    const double* values,
    int n_values,
    const int* edge_from,
    const int* edge_to,
    int n_edges,
    int method,
    double& statistic,
    double& score) {
  int n_finite = 0;
  long double sum = 0.0L;
  for (int j = 0; j < n_values; ++j) {
    if (std::isfinite(values[j])) {
      ++n_finite;
      sum += static_cast<long double>(values[j]);
    }
  }
  if (n_finite < 3) {
    return;
  }

  long double mean_extended = sum / static_cast<long double>(n_finite);
  long double mean_correction = 0.0L;
  for (int j = 0; j < n_values; ++j) {
    if (std::isfinite(values[j])) {
      mean_correction += static_cast<long double>(values[j]) - mean_extended;
    }
  }
  mean_extended += mean_correction / static_cast<long double>(n_finite);
  const double mean = static_cast<double>(mean_extended);
  long double denominator_sum = 0.0L;
  for (int j = 0; j < n_values; ++j) {
    if (std::isfinite(values[j])) {
      const double centered = values[j] - mean;
      denominator_sum += static_cast<long double>(centered * centered);
    }
  }
  const double denominator = static_cast<double>(denominator_sum);
  if (!std::isfinite(denominator) || denominator <= 0.0) {
    return;
  }

  int finite_edges = 0;
  long double numerator_sum = 0.0L;
  for (int edge = 0; edge < n_edges; ++edge) {
    const double from_value = values[edge_from[edge]];
    const double to_value = values[edge_to[edge]];
    if (!std::isfinite(from_value) || !std::isfinite(to_value)) {
      continue;
    }
    ++finite_edges;
    if (method == 1) {
      numerator_sum += static_cast<long double>(
        (from_value - mean) * (to_value - mean)
      );
    } else {
      const double difference = from_value - to_value;
      numerator_sum += static_cast<long double>(difference * difference);
    }
  }
  if (finite_edges == 0) {
    return;
  }

  const double numerator = static_cast<double>(numerator_sum);
  if (method == 1) {
    statistic = static_cast<double>(n_finite) /
      static_cast<double>(finite_edges) * numerator / denominator;
    score = statistic;
  } else {
    statistic = static_cast<double>(n_finite - 1) /
      (2.0 * static_cast<double>(finite_edges)) * numerator / denominator;
    score = 1.0 - statistic;
  }
}

inline void score_row_with_permutations(
    const double* values,
    int n_values,
    const int* edge_from,
    const int* edge_to,
    int n_edges,
    int method,
    int n_permutations,
    double& statistic,
    double& score,
    double& p_value,
    std::vector<double>& pool,
    std::vector<double>& permuted) {
  score_row(
    values,
    n_values,
    edge_from,
    edge_to,
    n_edges,
    method,
    statistic,
    score
  );
  if (n_permutations <= 0 || !std::isfinite(statistic)) {
    return;
  }

  int valid = 0;
  int extreme = 0;
  pool.assign(values, values + n_values);
  permuted.resize(n_values);
  for (int permutation = 0; permutation < n_permutations; ++permutation) {
    std::copy(values, values + n_values, pool.begin());
    int remaining = n_values;
    for (int position = 0; position < n_values; ++position) {
      const int selected = static_cast<int>(R_unif_index(
        static_cast<double>(remaining)
      ));
      permuted[position] = pool[selected];
      pool[selected] = pool[remaining - 1];
      --remaining;
    }
    double permutation_statistic = NA_REAL;
    double permutation_score = NA_REAL;
    score_row(
      permuted.data(),
      n_values,
      edge_from,
      edge_to,
      n_edges,
      method,
      permutation_statistic,
      permutation_score
    );
    if (!std::isfinite(permutation_statistic)) {
      continue;
    }
    ++valid;
    if (
      (method == 1 && permutation_statistic >= statistic) ||
      (method == 2 && permutation_statistic <= statistic)
    ) {
      ++extreme;
    }
  }
  if (valid > 0) {
    p_value = static_cast<double>(extreme + 1) /
      static_cast<double>(valid + 1);
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List spatial_variable_score_cpp(
    SEXP expr,
    Rcpp::IntegerVector edge_from,
    Rcpp::IntegerVector edge_to,
    int method,
    int n_permutations = 0) {
  if (edge_from.size() != edge_to.size()) {
    Rcpp::stop("`edge_from` and `edge_to` must have the same length.");
  }
  if (method != 1 && method != 2) {
    Rcpp::stop("`method` must be 1 (Moran) or 2 (Geary).");
  }
  if (n_permutations < 0) {
    Rcpp::stop("`n_permutations` must be non-negative.");
  }

  const bool is_sparse = Rf_inherits(expr, "dgCMatrix");
  const bool is_dense = Rf_isMatrix(expr) && TYPEOF(expr) == REALSXP;
  if (!is_sparse && !is_dense) {
    Rcpp::stop("`expr` must be a numeric matrix or a `dgCMatrix`.");
  }

  int n_features = 0;
  int n_spots = 0;
  if (is_sparse) {
    Rcpp::S4 sparse(expr);
    Rcpp::IntegerVector dimensions = sparse.slot("Dim");
    n_features = dimensions[0];
    n_spots = dimensions[1];
  } else {
    Rcpp::NumericMatrix dense(expr);
    n_features = dense.nrow();
    n_spots = dense.ncol();
  }

  const int n_edges = edge_from.size();
  for (int edge = 0; edge < n_edges; ++edge) {
    if (edge_from[edge] < 0 || edge_from[edge] >= n_spots ||
        edge_to[edge] < 0 || edge_to[edge] >= n_spots) {
      Rcpp::stop("Edge indices are outside the expression matrix columns.");
    }
  }

  Rcpp::NumericVector statistic(n_features, NA_REAL);
  Rcpp::NumericVector score(n_features, NA_REAL);
  Rcpp::NumericVector p_value(n_features, NA_REAL);
  double* statistic_ptr = statistic.begin();
  double* score_ptr = score.begin();
  double* p_value_ptr = p_value.begin();
  const int* edge_from_ptr = edge_from.begin();
  const int* edge_to_ptr = edge_to.begin();

  if (!is_sparse) {
    Rcpp::NumericMatrix dense(expr);
    const double* dense_ptr = dense.begin();

#ifdef _OPENMP
#pragma omp parallel if(n_permutations == 0)
#endif
    {
      std::vector<double> row(n_spots);
      std::vector<double> pool;
      std::vector<double> permuted;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int feature = 0; feature < n_features; ++feature) {
        for (int spot = 0; spot < n_spots; ++spot) {
          row[spot] = dense_ptr[feature + n_features * spot];
        }
        score_row_with_permutations(
          row.data(),
          n_spots,
          edge_from_ptr,
          edge_to_ptr,
          n_edges,
          method,
          n_permutations,
          statistic_ptr[feature],
          score_ptr[feature],
          p_value_ptr[feature],
          pool,
          permuted
        );
      }
    }
  } else {
    Rcpp::S4 sparse(expr);
    Rcpp::IntegerVector column_ptr = sparse.slot("p");
    Rcpp::IntegerVector row_index = sparse.slot("i");
    Rcpp::NumericVector sparse_value = sparse.slot("x");
    const int nnz = sparse_value.size();

    std::vector<int> row_ptr(n_features + 1, 0);
    for (int index = 0; index < nnz; ++index) {
      ++row_ptr[row_index[index] + 1];
    }
    for (int feature = 0; feature < n_features; ++feature) {
      row_ptr[feature + 1] += row_ptr[feature];
    }

    std::vector<int> cursor(row_ptr);
    std::vector<int> spot_index(nnz);
    std::vector<double> values(nnz);
    for (int spot = 0; spot < n_spots; ++spot) {
      for (int index = column_ptr[spot]; index < column_ptr[spot + 1]; ++index) {
        const int feature = row_index[index];
        const int destination = cursor[feature]++;
        spot_index[destination] = spot;
        values[destination] = sparse_value[index];
      }
    }

#ifdef _OPENMP
#pragma omp parallel if(n_permutations == 0)
#endif
    {
      std::vector<double> row(n_spots, 0.0);
      std::vector<double> pool;
      std::vector<double> permuted;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int feature = 0; feature < n_features; ++feature) {
        std::fill(row.begin(), row.end(), 0.0);
        for (int index = row_ptr[feature]; index < row_ptr[feature + 1]; ++index) {
          row[spot_index[index]] += values[index];
        }
        score_row_with_permutations(
          row.data(),
          n_spots,
          edge_from_ptr,
          edge_to_ptr,
          n_edges,
          method,
          n_permutations,
          statistic_ptr[feature],
          score_ptr[feature],
          p_value_ptr[feature],
          pool,
          permuted
        );
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["statistic"] = statistic,
    Rcpp::_["score"] = score,
    Rcpp::_["p_value"] = p_value
  );
}
