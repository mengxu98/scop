// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <R_ext/RS.h>
#include <thisutils/log_message.h>
#include "velocity_utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace arma;

extern "C" {
void F77_CALL(dtrexc)(
    const char* compq,
    const arma::blas_int* n,
    double* t,
    const arma::blas_int* ldt,
    double* q,
    const arma::blas_int* ldq,
    arma::blas_int* ifst,
    arma::blas_int* ilst,
    double* work,
    arma::blas_int* info
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
    , const arma::blas_len compq_len
#endif
);
}

namespace {

// CellRank 2.0.7 ranks each CSR row with NumPy's default argsort before
// reversing the order.  Equal connectivity weights therefore follow NumPy's
// intentionally unstable scalar aquicksort order, not a stable column order.
// Reproducing that order matters whenever the hard-threshold cutoff crosses a
// tie (a common case for shared-nearest-neighbour graphs).
inline bool numpy126_double_less(double lhs, double rhs) {
  return lhs < rhs || (std::isnan(rhs) && !std::isnan(lhs));
}

std::vector<int> numpy126_aquicksort_indices(
    const std::vector<double>& values)
{
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(values.size());
  std::vector<int> order(static_cast<size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  if (n <= 1) return order;

  const auto less_index = [&](int lhs, int rhs) {
    return numpy126_double_less(values[lhs], values[rhs]);
  };
  const auto heap_sort = [&](std::ptrdiff_t offset, std::ptrdiff_t count) {
    const auto at = [&](std::ptrdiff_t one_based) -> int& {
      return order[static_cast<size_t>(offset + one_based - 1)];
    };
    std::ptrdiff_t i = 0;
    std::ptrdiff_t j = 0;
    std::ptrdiff_t l = 0;
    int tmp = 0;
    for (l = count >> 1; l > 0; --l) {
      tmp = at(l);
      for (i = l, j = l << 1; j <= count;) {
        if (j < count && less_index(at(j), at(j + 1))) ++j;
        if (less_index(tmp, at(j))) {
          at(i) = at(j);
          i = j;
          j += j;
        } else {
          break;
        }
      }
      at(i) = tmp;
    }
    for (std::ptrdiff_t remaining = count; remaining > 1;) {
      tmp = at(remaining);
      at(remaining) = at(1);
      --remaining;
      for (i = 1, j = 2; j <= remaining;) {
        if (j < remaining && less_index(at(j), at(j + 1))) ++j;
        if (less_index(tmp, at(j))) {
          at(i) = at(j);
          i = j;
          j += j;
        } else {
          break;
        }
      }
      at(i) = tmp;
    }
  };

  int most_significant_bit = 0;
  for (std::ptrdiff_t value = n; (value >>= 1);) ++most_significant_bit;
  int depth = most_significant_bit * 2;
  std::ptrdiff_t left = 0;
  std::ptrdiff_t right = n - 1;
  std::vector<std::ptrdiff_t> left_stack;
  std::vector<std::ptrdiff_t> right_stack;
  std::vector<int> depth_stack;
  left_stack.reserve(128);
  right_stack.reserve(128);
  depth_stack.reserve(128);

  for (;;) {
    if (depth < 0) {
      heap_sort(left, right - left + 1);
      goto stack_pop;
    }
    while ((right - left) > 15) {
      const std::ptrdiff_t middle = left + ((right - left) >> 1);
      if (less_index(order[middle], order[left])) {
        std::swap(order[middle], order[left]);
      }
      if (less_index(order[right], order[middle])) {
        std::swap(order[right], order[middle]);
      }
      if (less_index(order[middle], order[left])) {
        std::swap(order[middle], order[left]);
      }

      const int pivot = order[middle];
      std::ptrdiff_t lower = left;
      std::ptrdiff_t upper = right - 1;
      std::swap(order[middle], order[upper]);
      for (;;) {
        do { ++lower; } while (less_index(order[lower], pivot));
        do { --upper; } while (less_index(pivot, order[upper]));
        if (lower >= upper) break;
        std::swap(order[lower], order[upper]);
      }
      std::swap(order[lower], order[right - 1]);

      --depth;
      if (lower - left < right - lower) {
        left_stack.push_back(lower + 1);
        right_stack.push_back(right);
        depth_stack.push_back(depth);
        right = lower - 1;
      } else {
        left_stack.push_back(left);
        right_stack.push_back(lower - 1);
        depth_stack.push_back(depth);
        left = lower + 1;
      }
    }

    for (std::ptrdiff_t current = left + 1; current <= right; ++current) {
      const int value = order[current];
      std::ptrdiff_t destination = current;
      std::ptrdiff_t previous = current - 1;
      while (destination > left && less_index(value, order[previous])) {
        order[destination--] = order[previous--];
      }
      order[destination] = value;
    }

  stack_pop:
    if (left_stack.empty()) break;
    left = left_stack.back();
    left_stack.pop_back();
    right = right_stack.back();
    right_stack.pop_back();
    depth = depth_stack.back();
    depth_stack.pop_back();
  }
  return order;
}

arma::blas_int gpcca_select_no_eigenvalue(const double*, const double*) {
  return 0;
}

struct GpccaSchurBlock {
  int start;
  int size;
  double magnitude;
};

std::vector<GpccaSchurBlock> gpcca_schur_blocks(
    const std::vector<double>& schur_form,
    int n,
    int start)
{
  std::vector<GpccaSchurBlock> blocks;
  const double tolerance = 100.0 * std::numeric_limits<double>::epsilon();
  int index = start;
  while (index < n) {
    const bool is_pair = index + 1 < n &&
      std::abs(schur_form[(index + 1) + index * n]) > tolerance;
    if (!is_pair) {
      blocks.push_back({
        index,
        1,
        std::abs(schur_form[index + index * n])
      });
      ++index;
      continue;
    }
    const double a = schur_form[index + index * n];
    const double b = schur_form[index + (index + 1) * n];
    const double c = schur_form[(index + 1) + index * n];
    const double d = schur_form[(index + 1) + (index + 1) * n];
    blocks.push_back({
      index,
      2,
      std::sqrt(std::abs(a * d - b * c))
    });
    index += 2;
  }
  return blocks;
}

bool gpcca_dominant_real_schur(
    const NumericMatrix& matrix,
    int n_components,
    NumericMatrix& schur_vectors,
    NumericVector& eigenvalues_real,
    NumericVector& eigenvalues_imaginary)
{
  const int n = matrix.nrow();
  if (matrix.ncol() != n || n < 1) return false;
  const int requested = std::max(1, std::min(n_components, n));
  arma::blas_int lapack_n = static_cast<arma::blas_int>(n);
  arma::blas_int leading_dimension = lapack_n;
  char job_vectors = 'V';
  char no_sort = 'N';
  arma::blas_int selected_dimension = 0;
  arma::blas_int info = 0;

  std::vector<double> schur_form(matrix.begin(), matrix.end());
  std::vector<double> schur_basis(static_cast<size_t>(n) * n, 0.0);
  std::vector<double> real_values(n, 0.0);
  std::vector<double> imaginary_values(n, 0.0);
  std::vector<arma::blas_int> selected(n, 0);

  arma::blas_int workspace_size = -1;
  double workspace_query = 0.0;
  arma::lapack::gees(
    &job_vectors,
    &no_sort,
    reinterpret_cast<void*>(gpcca_select_no_eigenvalue),
    &lapack_n,
    schur_form.data(),
    &leading_dimension,
    &selected_dimension,
    real_values.data(),
    imaginary_values.data(),
    schur_basis.data(),
    &leading_dimension,
    &workspace_query,
    &workspace_size,
    selected.data(),
    &info
  );
  if (info != 0 || !std::isfinite(workspace_query)) return false;

  workspace_size = std::max<arma::blas_int>(
    3 * lapack_n,
    static_cast<arma::blas_int>(std::ceil(workspace_query))
  );
  std::vector<double> workspace(workspace_size, 0.0);
  schur_form.assign(matrix.begin(), matrix.end());
  info = 0;
  arma::lapack::gees(
    &job_vectors,
    &no_sort,
    reinterpret_cast<void*>(gpcca_select_no_eigenvalue),
    &lapack_n,
    schur_form.data(),
    &leading_dimension,
    &selected_dimension,
    real_values.data(),
    imaginary_values.data(),
    schur_basis.data(),
    &leading_dimension,
    workspace.data(),
    &workspace_size,
    selected.data(),
    &info
  );
  if (info != 0) return false;

  const char update_vectors = 'V';
  std::vector<double> reorder_workspace(std::max(1, n), 0.0);
  int destination = 0;
  while (destination < requested) {
    const int remaining = requested - destination;
    std::vector<GpccaSchurBlock> blocks = gpcca_schur_blocks(
      schur_form, n, destination
    );
    int best = -1;
    for (int block = 0; block < static_cast<int>(blocks.size()); ++block) {
      if (blocks[block].size > remaining) continue;
      if (best < 0 || blocks[block].magnitude > blocks[best].magnitude) {
        best = block;
      }
    }
    if (best < 0) return false;

    arma::blas_int from = static_cast<arma::blas_int>(blocks[best].start + 1);
    arma::blas_int to = static_cast<arma::blas_int>(destination + 1);
    info = 0;
    if (from != to) {
      F77_CALL(dtrexc)(
        &update_vectors,
        &lapack_n,
        schur_form.data(),
        &leading_dimension,
        schur_basis.data(),
        &leading_dimension,
        &from,
        &to,
        reorder_workspace.data(),
        &info
#if defined(ARMA_USE_FORTRAN_HIDDEN_ARGS)
        , static_cast<arma::blas_len>(1)
#endif
      );
      if (info != 0) return false;
    }
    const std::vector<GpccaSchurBlock> moved = gpcca_schur_blocks(
      schur_form, n, destination
    );
    if (moved.empty()) return false;
    destination += moved.front().size;
  }

  schur_vectors = NumericMatrix(n, requested);
  eigenvalues_real = NumericVector(requested);
  eigenvalues_imaginary = NumericVector(requested);
  int column = 0;
  while (column < requested) {
    for (int row = 0; row < n; ++row) {
      schur_vectors(row, column) = schur_basis[row + column * n];
    }
    const bool is_pair = column + 1 < requested &&
      std::abs(schur_form[(column + 1) + column * n]) >
        100.0 * std::numeric_limits<double>::epsilon();
    if (!is_pair) {
      eigenvalues_real[column] = schur_form[column + column * n];
      ++column;
      continue;
    }
    for (int row = 0; row < n; ++row) {
      schur_vectors(row, column + 1) = schur_basis[row + (column + 1) * n];
    }
    const double a = schur_form[column + column * n];
    const double b = schur_form[column + (column + 1) * n];
    const double c = schur_form[(column + 1) + column * n];
    const double d = schur_form[(column + 1) + (column + 1) * n];
    const double real_part = 0.5 * (a + d);
    const double discriminant =
      0.25 * (a - d) * (a - d) + b * c;
    const double imaginary_part = std::sqrt(std::max(0.0, -discriminant));
    eigenvalues_real[column] = real_part;
    eigenvalues_real[column + 1] = real_part;
    eigenvalues_imaginary[column] = imaginary_part;
    eigenvalues_imaginary[column + 1] = -imaginary_part;
    column += 2;
  }
  return true;
}

std::vector<int> select_gpcca_terminal_representatives(
    const arma::mat& chi,
    int n_cells_terminal)
{
  const int n = chi.n_rows;
  const int M = chi.n_cols;
  const int requested = std::max(1, n_cells_terminal);
  std::vector<int> representative_lineage(n, -1);
  std::vector<char> used(n, 0);

  std::vector<std::pair<double, int>> lineage_priority;
  lineage_priority.reserve(M);
  for (int m = 0; m < M; ++m) {
    double best = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (std::isfinite(chi(i, m))) best = std::max(best, chi(i, m));
    }
    lineage_priority.push_back({best, m});
  }
  std::sort(
    lineage_priority.begin(), lineage_priority.end(),
    [](const std::pair<double, int>& left, const std::pair<double, int>& right) {
      if (left.first == right.first) return left.second < right.second;
      return left.first > right.first;
    }
  );

  int remaining_cells = n;
  for (int position = 0; position < M; ++position) {
    const int m = lineage_priority[position].second;
    const int remaining_lineages = M - position;
    const int target = std::max(
      1,
      std::min(requested, remaining_cells - (remaining_lineages - 1))
    );

    std::vector<std::pair<double, int>> candidates;
    candidates.reserve(n);
    for (int i = 0; i < n; ++i) {
      if (!used[i]) candidates.push_back({chi(i, m), i});
    }
    std::sort(
      candidates.begin(), candidates.end(),
      [](const std::pair<double, int>& left, const std::pair<double, int>& right) {
        if (left.first == right.first) return left.second < right.second;
        return left.first > right.first;
      }
    );

    int selected = 0;
    for (const auto& candidate : candidates) {
      if (selected >= target) break;
      const int cell = candidate.second;
      if (used[cell]) continue;
      used[cell] = 1;
      representative_lineage[cell] = m;
      ++selected;
      --remaining_cells;
    }
  }
  return representative_lineage;
}

bool fill_gpcca_rotation(
    const arma::vec& alpha,
    const arma::mat& features,
    arma::mat& rotation,
    double& objective)
{
  const int M = features.n_cols;
  const int k = M - 1;
  if (k < 1 || alpha.n_elem != static_cast<arma::uword>(k * k)) return false;

  rotation.zeros(M, M);
  for (int row = 0; row < k; ++row) {
    for (int column = 0; column < k; ++column) {
      rotation(row + 1, column + 1) = alpha[row * k + column];
    }
  }
  for (int row = 1; row < M; ++row) {
    rotation(row, 0) = -arma::accu(rotation.submat(row, 1, row, M - 1));
  }

  const arma::mat dummy = -features.cols(1, M - 1) * rotation.rows(1, M - 1);
  for (int column = 0; column < M; ++column) {
    rotation(0, column) = dummy.col(column).max();
    if (!std::isfinite(rotation(0, column))) return false;
  }
  const double scale = arma::accu(rotation.row(0));
  if (!std::isfinite(scale) || std::abs(scale) <= 1e-15) return false;
  rotation /= scale;
  if (rotation.row(0).min() <= 0.0) return false;

  objective = static_cast<double>(M);
  for (int column = 0; column < M; ++column) {
    const double numerator = arma::dot(
      rotation.col(column), rotation.col(column)
    );
    objective -= numerator / rotation(0, column);
  }
  return std::isfinite(objective);
}

bool optimize_gpcca_rotation(
    const arma::mat& features,
    arma::mat& rotation,
    int& evaluations,
    double& initial_objective,
    double& final_objective)
{
  const int M = features.n_cols;
  const int k = M - 1;
  const int dimension = k * k;
  if (dimension < 1) return false;

  arma::vec initial(dimension);
  for (int row = 0; row < k; ++row) {
    for (int column = 0; column < k; ++column) {
      initial[row * k + column] = rotation(row + 1, column + 1);
    }
  }

  arma::mat feasible;
  if (!fill_gpcca_rotation(
      initial, features, feasible, initial_objective)) return false;

  arma::mat simplex(dimension + 1, dimension);
  simplex.row(0) = initial.t();
  for (int coordinate = 0; coordinate < dimension; ++coordinate) {
    simplex.row(coordinate + 1) = initial.t();
    double& value = simplex(coordinate + 1, coordinate);
    value = value != 0.0 ? 1.05 * value : 0.00025;
  }

  arma::vec values(dimension + 1);
  arma::mat candidate_rotation;
  evaluations = 0;
  for (int point = 0; point <= dimension; ++point) {
    double value = std::numeric_limits<double>::infinity();
    fill_gpcca_rotation(
      simplex.row(point).t(), features, candidate_rotation, value
    );
    values[point] = value;
    ++evaluations;
  }

  const int max_evaluations = std::max(200, 200 * dimension);
  const double x_tolerance = 1e-4;
  const double f_tolerance = 1e-4;
  bool converged = false;

  while (evaluations < max_evaluations) {
    const arma::uvec ordering = arma::sort_index(values);
    simplex = simplex.rows(ordering);
    values = values(ordering);

    double coordinate_spread = 0.0;
    for (int point = 1; point <= dimension; ++point) {
      coordinate_spread = std::max(
        coordinate_spread,
        arma::abs(simplex.row(point) - simplex.row(0)).max()
      );
    }
    const double objective_spread =
      arma::abs(values.tail(dimension) - values[0]).max();
    if (coordinate_spread <= x_tolerance && objective_spread <= f_tolerance) {
      converged = true;
      break;
    }

    const arma::rowvec centroid =
      arma::mean(simplex.rows(0, dimension - 1), 0);
    const arma::rowvec reflected =
      2.0 * centroid - simplex.row(dimension);
    double reflected_value = std::numeric_limits<double>::infinity();
    fill_gpcca_rotation(
      reflected.t(), features, candidate_rotation, reflected_value
    );
    ++evaluations;

    bool shrink = false;
    if (reflected_value < values[0]) {
      const arma::rowvec expanded =
        3.0 * centroid - 2.0 * simplex.row(dimension);
      double expanded_value = std::numeric_limits<double>::infinity();
      fill_gpcca_rotation(
        expanded.t(), features, candidate_rotation, expanded_value
      );
      ++evaluations;
      if (expanded_value < reflected_value) {
        simplex.row(dimension) = expanded;
        values[dimension] = expanded_value;
      } else {
        simplex.row(dimension) = reflected;
        values[dimension] = reflected_value;
      }
    } else if (reflected_value < values[dimension - 1]) {
      simplex.row(dimension) = reflected;
      values[dimension] = reflected_value;
    } else if (reflected_value < values[dimension]) {
      const arma::rowvec contracted =
        1.5 * centroid - 0.5 * simplex.row(dimension);
      double contracted_value = std::numeric_limits<double>::infinity();
      fill_gpcca_rotation(
        contracted.t(), features, candidate_rotation, contracted_value
      );
      ++evaluations;
      if (contracted_value <= reflected_value) {
        simplex.row(dimension) = contracted;
        values[dimension] = contracted_value;
      } else {
        shrink = true;
      }
    } else {
      const arma::rowvec contracted =
        0.5 * centroid + 0.5 * simplex.row(dimension);
      double contracted_value = std::numeric_limits<double>::infinity();
      fill_gpcca_rotation(
        contracted.t(), features, candidate_rotation, contracted_value
      );
      ++evaluations;
      if (contracted_value < values[dimension]) {
        simplex.row(dimension) = contracted;
        values[dimension] = contracted_value;
      } else {
        shrink = true;
      }
    }

    if (shrink) {
      for (int point = 1;
           point <= dimension && evaluations < max_evaluations;
           ++point) {
        simplex.row(point) = simplex.row(0) +
          0.5 * (simplex.row(point) - simplex.row(0));
        double value = std::numeric_limits<double>::infinity();
        fill_gpcca_rotation(
          simplex.row(point).t(), features, candidate_rotation, value
        );
        values[point] = value;
        ++evaluations;
      }
    }
  }

  const arma::uword best = values.index_min();
  if (!fill_gpcca_rotation(
      simplex.row(best).t(), features, rotation, final_objective)) return false;
  return converged;
}

bool build_gpcca_isa_membership(
    const NumericMatrix& schur_vectors,
    const NumericVector& input_distribution,
    arma::mat& membership,
    IntegerVector& macro,
    IntegerVector& simplex_indices,
    int& optimization_evaluations,
    double& objective_initial,
    double& objective_final,
    bool& optimization_converged)
{
  const int n = schur_vectors.nrow();
  const int M = schur_vectors.ncol();
  if (M < 2 || input_distribution.size() != n) return false;

  arma::vec weights(n);
  double weight_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    double value = input_distribution[i];
    if (!std::isfinite(value) || value < 0.0) value = 0.0;
    weights(i) = value;
    weight_sum += value;
  }
  if (!std::isfinite(weight_sum) || weight_sum <= 1e-15) {
    weights.fill(1.0 / static_cast<double>(n));
  } else {
    weights /= weight_sum;
  }

  // Transform the balanced eigenvectors back to right eigenvectors and use an
  // input-distribution-weighted modified Gram-Schmidt basis with a constant Perron
  // vector.  The subsequent inner-simplex rotation is invariant to sign
  // choices and substantially more faithful than abs(Schur) memberships.
  arma::mat features(n, M, arma::fill::zeros);
  features.col(0).ones();
  for (int column = 1; column < M; ++column) {
    arma::vec value(n);
    for (int i = 0; i < n; ++i) {
      const double scale = std::sqrt(std::max(weights(i), 1e-15));
      double entry = schur_vectors(i, column) / scale;
      value(i) = std::isfinite(entry) ? entry : 0.0;
    }
    for (int previous = 0; previous < column; ++previous) {
      const double projection = arma::dot(weights % value, features.col(previous));
      value -= projection * features.col(previous);
    }
    const double norm = std::sqrt(arma::dot(weights, arma::square(value)));
    if (!std::isfinite(norm) || norm <= 1e-12) return false;
    features.col(column) = value / norm;
  }

  arma::mat working = features;
  std::vector<arma::uword> indices(M, 0);
  arma::vec distances = arma::sqrt(arma::sum(arma::square(working), 1));
  indices[0] = distances.index_max();
  working.each_row() -= working.row(indices[0]);
  for (int column = 1; column < M; ++column) {
    arma::rowvec direction = working.row(indices[column - 1]);
    for (int i = 0; i < n; ++i) {
      const double projection = arma::dot(working.row(i), direction);
      working.row(i) -= projection * direction;
    }
    distances = arma::sqrt(arma::sum(arma::square(working), 1));
    indices[column] = distances.index_max();
    const double scale = distances(indices[column]);
    if (!std::isfinite(scale) || scale <= 1e-12) return false;
    working /= scale;
  }

  arma::mat vertices(M, M);
  for (int row = 0; row < M; ++row) vertices.row(row) = features.row(indices[row]);
  arma::mat rotation;
  if (!arma::pinv(rotation, vertices) || !rotation.is_finite()) return false;
  optimization_converged = optimize_gpcca_rotation(
    features,
    rotation,
    optimization_evaluations,
    objective_initial,
    objective_final
  );
  membership = features * rotation;
  if (!membership.is_finite()) return false;

  macro = IntegerVector(n);
  for (int i = 0; i < n; ++i) {
    double total = 0.0;
    for (int lineage = 0; lineage < M; ++lineage) {
      double value = membership(i, lineage);
      if (!std::isfinite(value) || value < 0.0) value = 0.0;
      membership(i, lineage) = value;
      total += value;
    }
    if (total <= 1e-12) {
      for (int lineage = 0; lineage < M; ++lineage) {
        membership(i, lineage) = 1.0 / static_cast<double>(M);
      }
    } else {
      membership.row(i) /= total;
    }
    macro[i] = membership.row(i).index_max() + 1;
  }

  simplex_indices = IntegerVector(M);
  for (int lineage = 0; lineage < M; ++lineage) {
    simplex_indices[lineage] = static_cast<int>(indices[lineage]) + 1;
  }
  return true;
}

bool solve_gpcca_cell_absorption(
    const NumericMatrix& transition,
    const arma::mat& chi,
    const std::vector<int>& representative_lineage,
    NumericMatrix& absorption,
    int& iterations,
    double& residual,
    std::string& method)
{
  const int n = transition.nrow();
  const int M = chi.n_cols;
  std::vector<int> transient_index;
  transient_index.reserve(n);
  std::vector<int> transient_position(n, -1);
  for (int i = 0; i < n; ++i) {
    if (representative_lineage[i] < 0) {
      transient_position[i] = transient_index.size();
      transient_index.push_back(i);
    }
  }

  const int nQ = transient_index.size();
  absorption = NumericMatrix(n, M);
  if (nQ == 0) {
    for (int i = 0; i < n; ++i) {
      const int lineage = representative_lineage[i];
      if (lineage >= 0 && lineage < M) absorption(i, lineage) = 1.0;
    }
    iterations = 0;
    residual = 0.0;
    method = "terminal_identity";
    return true;
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(nQ * 24);
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nQ, M);
  for (int qi = 0; qi < nQ; ++qi) {
    const int i = transient_index[qi];
    triplets.emplace_back(qi, qi, 1.0);
    for (int j = 0; j < n; ++j) {
      const double weight = transition(i, j);
      if (!std::isfinite(weight) || weight <= 1e-15) continue;
      const int lineage = representative_lineage[j];
      if (lineage >= 0 && lineage < M) {
        rhs(qi, lineage) += weight;
      } else {
        const int qj = transient_position[j];
        if (qj >= 0) triplets.emplace_back(qi, qj, -weight);
      }
    }
  }

  Eigen::SparseMatrix<double> system(nQ, nQ);
  system.setFromTriplets(triplets.begin(), triplets.end());
  system.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(system);
  solver.factorize(system);
  bool solved = solver.info() == Eigen::Success;
  Eigen::MatrixXd solution;
  if (solved) {
    solution = solver.solve(rhs);
    solved = solver.info() == Eigen::Success && solution.allFinite();
  }

  if (solved) {
    for (int qi = 0; qi < nQ; ++qi) {
      const int i = transient_index[qi];
      for (int m = 0; m < M; ++m) absorption(i, m) = solution(qi, m);
    }
    iterations = 1;
    residual = (system * solution - rhs).cwiseAbs().maxCoeff();
    method = "sparse_lu";
  } else {
    arma::mat probabilities = chi;
    for (int i = 0; i < n; ++i) {
      const int lineage = representative_lineage[i];
      if (lineage < 0) continue;
      probabilities.row(i).zeros();
      probabilities(i, lineage) = 1.0;
    }
    arma::mat next = probabilities;
    residual = std::numeric_limits<double>::infinity();
    const int max_iterations = 5000;
    for (iterations = 1; iterations <= max_iterations; ++iterations) {
      residual = 0.0;
      for (int qi = 0; qi < nQ; ++qi) {
        const int i = transient_index[qi];
        next.row(i).zeros();
        for (int j = 0; j < n; ++j) {
          const double weight = transition(i, j);
          if (std::isfinite(weight) && weight > 1e-15) {
            next.row(i) += weight * probabilities.row(j);
          }
        }
        residual = std::max(
          residual,
          arma::abs(next.row(i) - probabilities.row(i)).max()
        );
      }
      probabilities.swap(next);
      if (residual < 1e-8) break;
    }
    for (int i = 0; i < n; ++i)
      for (int m = 0; m < M; ++m)
        absorption(i, m) = probabilities(i, m);
    method = "fixed_point";
    solved = residual < 1e-6;
  }

  for (int i = 0; i < n; ++i) {
    const int lineage = representative_lineage[i];
    if (lineage >= 0) {
      for (int m = 0; m < M; ++m) absorption(i, m) = (m == lineage) ? 1.0 : 0.0;
      continue;
    }
    double total = 0.0;
    for (int m = 0; m < M; ++m) {
      double value = absorption(i, m);
      if (!std::isfinite(value) || value < 0.0) value = 0.0;
      absorption(i, m) = value;
      total += value;
    }
    if (total <= 1e-12) {
      for (int m = 0; m < M; ++m) absorption(i, m) = chi(i, m);
      total = 1.0;
    }
    for (int m = 0; m < M; ++m) absorption(i, m) /= total;
  }
  return solved;
}

} // namespace

// ── 0. Sparse hard-threshold pseudotime kernel ─────────────────────────────

// [[Rcpp::export]]
Eigen::SparseMatrix<double> cellrank_hard_threshold_kernel_cpp(
    const Eigen::MappedSparseMatrix<double> connectivities,
    NumericVector pseudotime,
    double frac_to_keep = 0.3,
    bool backward = false)
{
  const int n_cells = connectivities.rows();
  if (connectivities.cols() != n_cells || pseudotime.size() != n_cells) {
    stop("CellRank connectivities and pseudotime dimensions do not agree");
  }
  if (!std::isfinite(frac_to_keep) || frac_to_keep < 0.0 || frac_to_keep > 1.0) {
    stop("frac_to_keep must be between 0 and 1");
  }

  std::vector<double> time(static_cast<size_t>(n_cells));
  double maximum_time = -std::numeric_limits<double>::infinity();
  for (int cell = 0; cell < n_cells; ++cell) {
    const double value = pseudotime[cell];
    if (!std::isfinite(value)) stop("CellRank pseudotime contains non-finite values");
    time[cell] = value;
    maximum_time = std::max(maximum_time, value);
  }
  if (backward) {
    for (double& value : time) value = maximum_time - value;
  }

  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> RowSparseMatrix;
  RowSparseMatrix row_connectivities = connectivities;
  row_connectivities.makeCompressed();
  std::vector<Eigen::Triplet<double> > triplets;
  triplets.reserve(static_cast<size_t>(connectivities.nonZeros()) + n_cells);

  std::vector<int> columns;
  std::vector<double> weights;
  std::vector<unsigned char> retained;
  for (int cell = 0; cell < n_cells; ++cell) {
    columns.clear();
    weights.clear();
    for (RowSparseMatrix::InnerIterator entry(row_connectivities, cell);
         entry; ++entry) {
      const double weight = entry.value();
      if (!std::isfinite(weight) || weight <= 0.0) continue;
      columns.push_back(entry.col());
      weights.push_back(weight);
    }

    const int degree = static_cast<int>(weights.size());
    if (degree == 0) {
      triplets.emplace_back(cell, cell, 1.0);
      continue;
    }
    const int n_close = std::max(
      0,
      std::min(30, static_cast<int>(std::floor(degree * frac_to_keep)))
    );
    const std::vector<int> ascending = numpy126_aquicksort_indices(weights);
    retained.assign(static_cast<size_t>(degree), 0);
    double row_sum = 0.0;
    for (int rank = 0; rank < degree; ++rank) {
      const int index = ascending[static_cast<size_t>(degree - rank - 1)];
      const bool keep = rank < n_close ||
        time[cell] <= time[columns[static_cast<size_t>(index)]];
      if (!keep) continue;
      retained[static_cast<size_t>(index)] = 1;
      row_sum += weights[static_cast<size_t>(index)];
    }

    if (!std::isfinite(row_sum) || row_sum <= 1e-12) {
      triplets.emplace_back(cell, cell, 1.0);
      continue;
    }
    const double inverse_sum = 1.0 / row_sum;
    for (int index = 0; index < degree; ++index) {
      if (retained[static_cast<size_t>(index)] == 0) continue;
      triplets.emplace_back(
        cell,
        columns[static_cast<size_t>(index)],
        weights[static_cast<size_t>(index)] * inverse_sum
      );
    }
  }

  Eigen::SparseMatrix<double> transition(n_cells, n_cells);
  transition.setFromTriplets(triplets.begin(), triplets.end());
  transition.makeCompressed();
  return transition;
}

// ── 1. Transition matrix validation & normalization ──────────────────────────

// [[Rcpp::export]]
List cellrank_validate_transition_matrix_cpp(
    NumericMatrix T_,
    double eps = 1e-10,
    double min_self_loop = 0.01)
{
  int n = T_.nrow();
  if (T_.ncol() != n) thisutils::log_message("Transition matrix must be square", "error");

  NumericMatrix T = clone(T_);
  int nans = 0, negs = 0, zeros = 0, loops = 0;
  scop_util::validate_transition_matrix(T, eps, min_self_loop, &nans, &negs, &zeros, &loops);

  return List::create(
    _["transition_matrix"] = T,
    _["nans_fixed"] = nans, _["negs_clipped"] = negs,
    _["zero_rows_fixed"] = zeros, _["self_loops_added"] = loops
  );
}

// ── 2. Stationary distribution ───────────────────────────────────────────────

// [[Rcpp::export]]
NumericVector cellrank_stationary_distribution_cpp(
    NumericMatrix T_, int max_iter = 1000, double tol = 1e-10)
{
  return scop_util::stationary_distribution(T_, max_iter, tol);
}

// ── 3. GPCCA real Schur decomposition ──────────────────────────────────────

// [[Rcpp::export]]
List cellrank_schur_cpp(NumericMatrix T_, int n_components = 10)
{
  int n = T_.nrow();
  if (n_components < 2) n_components = 2;
  if (n_components > n) n_components = n;

  // CellRank/pyGPCCA use a uniform input distribution by default for the
  // weighted Schur decomposition. This is distinct from the stationary
  // distribution, which is still returned as an estimator diagnostic.
  NumericVector pi = cellrank_stationary_distribution_cpp(T_, 200, 1e-8);
  NumericVector eta(n, 1.0 / static_cast<double>(n));

  // Build T_bar = diag(sqrt(eta)) * T * diag(1/sqrt(eta)).
  NumericMatrix T_bar(n, n);
  for (int i = 0; i < n; ++i) {
    double si = std::sqrt(std::max(eta[i], 1e-15));
    for (int j = 0; j < n; ++j) {
      double sj = std::sqrt(std::max(eta[j], 1e-15));
      T_bar(i, j) = si * T_(i, j) / sj;
    }
  }

  NumericMatrix schur_vecs;
  NumericVector eigenvalues;
  NumericVector eigenvalues_imaginary;
  if (!gpcca_dominant_real_schur(
      T_bar,
      n_components,
      schur_vecs,
      eigenvalues,
      eigenvalues_imaginary)) {
    stop("Unable to compute the dominant real Schur basis for CellRank GPCCA");
  }
  const int nc = schur_vecs.ncol();

  // Macrostate assignment: largest non-trivial component per row. The leading
  // Perron vector is close to constant and can otherwise collapse all cells into
  // a single macrostate on connectivity kernels.
  IntegerVector macro(n);
  int start_component = (nc > 1 && std::abs(eigenvalues[0] - 1.0) < 1e-6) ? 1 : 0;
  for (int i = 0; i < n; ++i) {
    int best = start_component; double best_val = -1;
    for (int c = start_component; c < nc; ++c) {
      double v = std::abs(schur_vecs(i, c));
      if (v > best_val) { best_val = v; best = c; }
    }
    macro[i] = best + 1;
  }

  return List::create(
    _["eigenvalues"] = eigenvalues, _["schur_vectors"] = schur_vecs,
    _["eigenvalues_imaginary"] = eigenvalues_imaginary,
    _["macrostate_assignment"] = macro, _["n_macrostates"] = nc,
    _["stationary_distribution"] = pi,
    _["method"] = "lapack_real_schur"
  );
}

// ── 4. Auto-detect number of macrostates via eigengap ────────────────────────

// [[Rcpp::export]]
int cellrank_auto_n_states_cpp(NumericVector eigenvalues, int min_states = 2, int max_states = 20)
{
  int n = eigenvalues.size();
  if (n <= min_states) return min_states;

  double max_gap = 0.0;
  int best_k = min_states;
  for (int k = min_states; k < std::min(max_states, n - 1); ++k) {
    double gap = std::abs(eigenvalues[k] - eigenvalues[k + 1]);
    if (gap > max_gap) {
      max_gap = gap;
      best_k = k + 1;
    }
  }

  for (int k = min_states; k < std::min(max_states, n - 1); ++k) {
    double denom = std::abs(eigenvalues[k]);
    if (denom > 1e-10) {
      double rel_gap = std::abs(eigenvalues[k] - eigenvalues[k + 1]) / denom;
      if (rel_gap > 0.3) return k + 1;
    }
  }

  return best_k;
}

// ── 5. Velocity kernel ───────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_velocity_kernel_cpp(
    NumericMatrix velocity_embedding,   // cells × dims
    NumericMatrix embedding,            // cells × dims
    IntegerMatrix knn_idx,              // cells × k (1-based)
    bool backward = false,
    double softmax_scale = 4.0,
    int n_neighbors_velo = -1)
{
  const int n_cells = velocity_embedding.nrow();
  const int n_dims = velocity_embedding.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    double vn = 0.0;
    for (int d = 0; d < n_dims; ++d)
      vn += velocity_embedding(cell, d) * velocity_embedding(cell, d);
    vn = std::sqrt(vn);
    if (vn < 1e-10) { T(cell, cell) = 1.0; continue; }

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      double dot = 0.0, dn = 0.0;
      for (int d = 0; d < n_dims; ++d) {
        double delta = embedding(nb, d) - embedding(cell, d);
        if (backward) {
          dot += (-velocity_embedding(cell, d)) * delta;
        } else {
          dot += velocity_embedding(cell, d) * delta;
        }
        dn += delta * delta;
      }
      dn = std::sqrt(dn);
      if (dn < 1e-10) continue;

      double cosine = dot / (vn * dn);
      if (cosine > 0 && std::isfinite(cosine)) {
        double weight = std::exp(cosine * softmax_scale);
        T(cell, nb) = weight;
        row_sum += weight;
      }
    }
    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// [[Rcpp::export]]
NumericMatrix cellrank_connectivity_kernel_cpp(
    IntegerMatrix knn_idx,
    NumericMatrix knn_dist)
{
  const int n_cells = knn_idx.nrow();
  const int n_neighbors = knn_idx.ncol();
  if (knn_dist.nrow() != n_cells || knn_dist.ncol() != n_neighbors) {
    thisutils::log_message("knn_idx and knn_dist must have the same dimensions", "error");
  }

  NumericMatrix T(n_cells, n_cells);
  std::vector<double> distances;
  distances.reserve(n_neighbors);
  for (int cell = 0; cell < n_cells; ++cell) {
    distances.clear();
    for (int k = 0; k < n_neighbors; ++k) {
      const double d = knn_dist(cell, k);
      if (!NumericVector::is_na(d)) distances.push_back(d);
    }
    double sigma = 1.0;
    if (!distances.empty()) {
      std::sort(distances.begin(), distances.end());
      const size_t middle = distances.size() / 2;
      sigma = distances.size() % 2 == 0 ?
        (distances[middle - 1] + distances[middle]) / 2.0 :
        distances[middle];
      sigma += 1e-10;
    }

    double row_sum = 0.0;
    for (int k = 0; k < n_neighbors; ++k) {
      int neighbor = knn_idx(cell, k);
      if (neighbor == NA_INTEGER) continue;
      neighbor -= 1;
      if (neighbor < 0 || neighbor >= n_cells || neighbor == cell) continue;
      const double d = knn_dist(cell, k);
      if (NumericVector::is_na(d)) continue;
      const double weight = std::exp(-(d * d) / (2.0 * sigma * sigma));
      T(cell, neighbor) = weight;
      row_sum += weight;
    }
    if (row_sum > 0.0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }
  return T;
}

// ── 5b. Gene-space velocity kernel (matching Python CellRank) ────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_velocity_kernel_gene_cpp(
    NumericMatrix gene_velocity,         // genes × cells
    NumericMatrix expression,            // genes × cells (Ms or counts)
    IntegerMatrix knn_idx,               // cells × k (1-based)
    bool backward = false,
    double softmax_scale = 4.0,
    int n_neighbors_velo = -1)
{
  const int n_genes = gene_velocity.nrow();
  const int n_cells = gene_velocity.ncol();
  const int n_neighbors = knn_idx.ncol();
  if (n_neighbors_velo <= 0) n_neighbors_velo = n_neighbors;

  // Precompute velocity norms per cell
  std::vector<double> vn_cells(n_cells, 0.0);
  for (int c = 0; c < n_cells; ++c) {
    double sq = 0.0;
    for (int g = 0; g < n_genes; ++g)
      sq += gene_velocity(g, c) * gene_velocity(g, c);
    vn_cells[c] = std::sqrt(sq);
  }

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    if (vn_cells[cell] < 1e-10) { T(cell, cell) = 1.0; continue; }

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors_velo && col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      // Cosine correlation in gene space:
      // cos(v_cell, expr_nb - expr_cell)
      double dot = 0.0, ndsq = 0.0;
      for (int g = 0; g < n_genes; ++g) {
        double delta = expression(g, nb) - expression(g, cell);
        double v = backward ? -gene_velocity(g, cell) : gene_velocity(g, cell);
        dot += v * delta;
        ndsq += delta * delta;
      }
      double nd = std::sqrt(ndsq);
      if (nd < 1e-10) continue;

      double cosine = dot / (vn_cells[cell] * nd);
      if (!std::isfinite(cosine)) continue;
      // CellRank uses exp(cosine / scale) for ALL neighbors
      double weight = std::exp(cosine / softmax_scale);
      T(cell, nb) = weight;
      row_sum += weight;
    }
    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// ── 6. Pseudotime kernel ─────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_pseudotime_kernel_cpp(
    NumericVector pseudotime,
    IntegerMatrix knn_idx,
    NumericVector cell_weights,
    double bandwidth = 1.0,
    bool backward = false)
{
  const int n_cells = pseudotime.size();
  const int n_neighbors = knn_idx.ncol();

  NumericMatrix T(n_cells, n_cells);

  for (int cell = 0; cell < n_cells; ++cell) {
    double pt_cell = pseudotime[cell];
    double w_cell = (cell_weights.size() >= n_cells) ? cell_weights[cell] : 1.0;

    double row_sum = 0.0;
    for (int col = 0; col < n_neighbors; ++col) {
      int nb = knn_idx(cell, col);
      if (nb == NA_INTEGER) continue;
      nb -= 1;
      if (nb < 0 || nb >= n_cells || nb == cell) continue;

      double pt_nb = pseudotime[nb];
      double w_nb = (cell_weights.size() >= n_cells) ? cell_weights[nb] : 1.0;

      double dt;
      if (backward) {
        dt = pt_cell - pt_nb;
      } else {
        dt = pt_nb - pt_cell;
      }

      if (dt <= 0) continue;

      double weight = w_nb * std::exp(-dt * dt / (2.0 * bandwidth * bandwidth));
      T(cell, nb) = weight;
      row_sum += weight;
    }

    if (row_sum > 0) {
      for (int j = 0; j < n_cells; ++j) T(cell, j) /= row_sum;
    } else {
      T(cell, cell) = 1.0;
    }
  }

  return T;
}

// ── 7. CytoTRACE kernel ───────────────────────────────────────────────────────

// [[Rcpp::export]]
NumericMatrix cellrank_cytotrace_kernel_cpp(
    NumericVector gene_counts,
    IntegerMatrix knn_idx,
    double bandwidth = 1.0,
    bool backward = false)
{
  const int n_cells = gene_counts.size();
  const int n_neighbors = knn_idx.ncol();

  double gc_min = gene_counts[0], gc_max = gene_counts[0];
  for (int i = 1; i < n_cells; ++i) {
    if (gene_counts[i] < gc_min) gc_min = gene_counts[i];
    if (gene_counts[i] > gc_max) gc_max = gene_counts[i];
  }
  NumericVector cytotrace_score(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    cytotrace_score[i] = (gc_max > gc_min)
      ? (gene_counts[i] - gc_min) / (gc_max - gc_min) : 0.5;
  }

  NumericVector weights(n_cells, 1.0);
  NumericMatrix T = cellrank_pseudotime_kernel_cpp(
    cytotrace_score, knn_idx, weights, bandwidth, backward);

  return T;
}

// ── 8. CFLARE estimator ──────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_cflare_cpp(
    NumericMatrix T_,
    int n_states = 5,
    int max_iter = 1000,
    double tol = 1e-6)
{
  int n = T_.nrow();
  if (n_states > n) n_states = n;

  List validated = cellrank_validate_transition_matrix_cpp(T_);
  NumericMatrix T = validated["transition_matrix"];

  NumericVector pi = cellrank_stationary_distribution_cpp(T, max_iter, tol);

  List schur = cellrank_schur_cpp(T, n_states);
  NumericVector eigenvalues = schur["eigenvalues"];
  NumericMatrix schur_vecs = schur["schur_vectors"];
  IntegerVector macro = schur["macrostate_assignment"];
  int M = schur["n_macrostates"];

  // Build coarse transition matrix
  NumericMatrix P_coarse(M, M);
  NumericVector sizes(M);
  for (int i = 0; i < n; ++i) sizes[macro[i] - 1] += 1.0;
  for (int i = 0; i < n; ++i) {
    int ai = macro[i] - 1;
    for (int j = 0; j < n; ++j) {
      int aj = macro[j] - 1;
      P_coarse(ai, aj) += T(i, j);
    }
  }
  for (int a = 0; a < M; ++a) {
    double rs = 0.0;
    for (int b = 0; b < M; ++b) rs += P_coarse(a, b);
    if (rs > 0) for (int b = 0; b < M; ++b) P_coarse(a, b) /= rs;
  }

  std::vector<int> is_terminal, term_idx, trans_idx;
  int n_terminal;
  scop_util::detect_terminal_states(P_coarse, M, is_terminal, term_idx, trans_idx, n_terminal);
  int nT = term_idx.size();
  int nQ = trans_idx.size();

  NumericMatrix abs_prob(n, std::max(nT, 1));
  IntegerVector term_state(n, NA_INTEGER);
  NumericVector fate_conf(n, 0.0);
  scop_util::compute_absorption_probabilities(P_coarse, macro, n, M, is_terminal, term_idx, trans_idx, nT, nQ, abs_prob, term_state, fate_conf);
  IntegerVector terminal_obs(n);
  for (int i = 0; i < n; ++i) {
    int m = macro[i] - 1;
    terminal_obs[i] = (m >= 0 && m < M && is_terminal[m]) ? (m + 1) : 0;
  }

  return List::create(
    _["transition_matrix"] = T,
    _["stationary_distribution"] = pi,
    _["eigenvalues"] = eigenvalues,
    _["schur_vectors"] = schur_vecs,
    _["macrostate_assignment"] = macro,
    _["n_macrostates"] = M,
    _["terminal_states"] = terminal_obs,
    _["lineage_assignment"] = term_state,
    _["fate_confidence"] = fate_conf,
    _["absorption_probabilities"] = abs_prob,
    _["n_terminal_states"] = nT,
    _["method"] = "CFLARE"
  );
}

// ── 9. GPCCA Estimator ────────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_gpcca_cpp(
    NumericMatrix T_,
    int n_states = 5,
    int n_cells_terminal = 10,
    bool skip_perron = false)
{
  int n = T_.nrow();

  List validated = cellrank_validate_transition_matrix_cpp(T_);
  NumericMatrix T = validated["transition_matrix"];

  NumericVector pi = cellrank_stationary_distribution_cpp(T, 200, 1e-8);
  NumericVector eta(n, 1.0 / static_cast<double>(n));

  int schur_components = skip_perron ? (n_states + 2) : n_states;
  List schur_result = cellrank_schur_cpp(T, schur_components);
  NumericVector eigenvalues_all = schur_result["eigenvalues"];
  NumericMatrix schur_vecs_all = schur_result["schur_vectors"];
  IntegerVector macro_all = schur_result["macrostate_assignment"];

  // If skip_perron, drop the first component (Perron vector, eigenvalue ≈ 1)
  int start_comp = skip_perron ? 1 : 0;
  int M = std::min(n_states, schur_vecs_all.ncol() - start_comp);
  NumericVector eigenvalues(M);
  NumericMatrix schur_vecs(n, M);
  for (int j = 0; j < M; ++j) {
    eigenvalues[j] = eigenvalues_all[j + start_comp];
    for (int i = 0; i < n; ++i)
      schur_vecs(i, j) = schur_vecs_all(i, j + start_comp);
  }
  // Recompute macro from the selected Schur vectors
  IntegerVector macro(n);
  for (int i = 0; i < n; ++i) {
    int best = 0; double best_val = -1;
    for (int c = 0; c < M; ++c) {
      double v = std::abs(schur_vecs(i, c));
      if (v > best_val) { best_val = v; best = c; }
    }
    macro[i] = best + 1;
  }

  // Approximate the PCCA rotation with the Inner Simplex Algorithm.  Fall back
  // to absolute Schur magnitudes only when the selected invariant basis is
  // numerically rank deficient.
  arma::mat chi(n, M);
  IntegerVector simplex_indices;
  int optimization_evaluations = 0;
  double objective_initial = NA_REAL;
  double objective_final = NA_REAL;
  bool optimization_converged = false;
  bool isa_converged = build_gpcca_isa_membership(
    schur_vecs,
    eta,
    chi,
    macro,
    simplex_indices,
    optimization_evaluations,
    objective_initial,
    objective_final,
    optimization_converged
  );
  if (!isa_converged) {
    chi.zeros();
    simplex_indices = IntegerVector();
    for (int i = 0; i < n; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < M; ++j) {
        double v = std::abs(schur_vecs(i, j));
        if (!std::isfinite(v)) v = 0.0;
        chi(i, j) = v;
        row_sum += v;
      }
      if (row_sum > 1e-10) {
        chi.row(i) /= row_sum;
      } else {
        chi.row(i).fill(1.0 / static_cast<double>(M));
      }
      macro[i] = chi.row(i).index_max() + 1;
    }
  }

  // Match pyGPCCA's coarse graining:
  // pinv(chi^T D_eta chi) * chi^T D_eta T chi.
  arma::vec eta_arma(n);
  for (int i = 0; i < n; ++i) eta_arma(i) = eta[i];
  arma::mat D_eta = arma::diagmat(eta_arma);
  arma::mat T_arma(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      T_arma(i, j) = T(i, j);

  arma::mat weighted_membership = chi.t() * D_eta;
  arma::mat P_coarse_arma = arma::pinv(weighted_membership * chi) *
    weighted_membership * T_arma * chi;
  for (int i = 0; i < M; ++i) {
    double rs = arma::accu(P_coarse_arma.row(i));
    if (rs > 1e-12) {
      P_coarse_arma.row(i) /= rs;
    } else {
      P_coarse_arma(i, i) = 1.0;
    }
  }
  NumericMatrix P_coarse = Rcpp::wrap(P_coarse_arma);

  // CellRank's GPCCA wrapper promotes representative cells from every computed
  // macrostate to terminal states. Solve absorption on the original cell-level
  // transition matrix instead of assigning one coarse probability vector to all
  // cells in the same hard macrostate.
  std::vector<int> representative_lineage =
    select_gpcca_terminal_representatives(chi, n_cells_terminal);
  NumericMatrix abs_prob;
  int absorption_iterations = 0;
  double absorption_residual = NA_REAL;
  std::string absorption_method;
  bool absorption_converged = solve_gpcca_cell_absorption(
    T, chi, representative_lineage, abs_prob,
    absorption_iterations, absorption_residual, absorption_method
  );

  IntegerVector terminal_obs(n, 0);
  IntegerVector term_state(n, NA_INTEGER);
  NumericVector fate_conf(n, 0.0);
  int terminal_cell_count = 0;
  for (int i = 0; i < n; ++i) {
    const int representative = representative_lineage[i];
    if (representative >= 0 && representative < M) {
      terminal_obs[i] = representative + 1;
      ++terminal_cell_count;
    }
    int best_lineage = 0;
    double best_probability = abs_prob(i, 0);
    for (int m = 1; m < M; ++m) {
      if (abs_prob(i, m) > best_probability) {
        best_probability = abs_prob(i, m);
        best_lineage = m;
      }
    }
    term_state[i] = best_lineage + 1;
    fate_conf[i] = best_probability;
  }

  return List::create(
    _["transition_matrix"] = T,
    _["stationary_distribution"] = pi,
    _["eigenvalues"] = eigenvalues,
    _["schur_vectors"] = schur_vecs,
    _["chi"] = wrap(chi),
    _["simplex_indices"] = simplex_indices,
    _["membership_method"] = !isa_converged ?
      "absolute_schur_fallback" :
      (optimization_converged ?
        "optimized_inner_simplex" : "inner_simplex_optimization_limit"),
    _["membership_optimization_evaluations"] = optimization_evaluations,
    _["membership_objective_initial"] = objective_initial,
    _["membership_objective_final"] = objective_final,
    _["membership_optimization_converged"] = optimization_converged,
    _["coarse_transition"] = wrap(P_coarse),
    _["macrostate_assignment"] = macro,
    _["n_macrostates"] = M,
    _["terminal_states"] = terminal_obs,
    _["lineage_assignment"] = term_state,
    _["fate_confidence"] = fate_conf,
    _["absorption_probabilities"] = abs_prob,
    _["n_terminal_states"] = M,
    _["n_terminal_cells"] = terminal_cell_count,
    _["absorption_method"] = absorption_method,
    _["absorption_iterations"] = absorption_iterations,
    _["absorption_residual"] = absorption_residual,
    _["absorption_converged"] = absorption_converged,
    _["method"] = "GPCCA"
  );
}

// ── 10. Lineage drivers ──────────────────────────────────────────────────────────

// [[Rcpp::export]]
List cellrank_lineage_drivers_cpp(
    NumericMatrix expression,
    NumericMatrix abs_probs,
    IntegerVector lineage_idx = IntegerVector())
{
  const int n_genes = expression.nrow();
  const int n_cells = expression.ncol();
  const int n_lineages = abs_probs.ncol();

  if (lineage_idx.size() == 0) {
    lineage_idx = IntegerVector(n_lineages);
    for (int i = 0; i < n_lineages; ++i) lineage_idx[i] = i + 1;
  }

  int n_out = lineage_idx.size();
  NumericMatrix corr(n_genes, n_out);
  NumericMatrix pval(n_genes, n_out);
  NumericVector means(n_genes);
  NumericVector vars(n_genes);

  for (int g = 0; g < n_genes; ++g) {
    double m = 0.0;
    for (int c = 0; c < n_cells; ++c) m += expression(g, c);
    m /= n_cells;
    means[g] = m;
    double v = 0.0;
    for (int c = 0; c < n_cells; ++c) v += (expression(g, c) - m) * (expression(g, c) - m);
    v /= (n_cells - 1);
    vars[g] = v;
  }

  for (int li = 0; li < n_out; ++li) {
    int l = lineage_idx[li] - 1;
    if (l < 0 || l >= n_lineages) continue;

    double mean_prob = 0.0;
    for (int c = 0; c < n_cells; ++c) mean_prob += abs_probs(c, l);
    mean_prob /= n_cells;
    double var_prob = 0.0;
    for (int c = 0; c < n_cells; ++c) {
      double d = abs_probs(c, l) - mean_prob;
      var_prob += d * d;
    }
    var_prob /= (n_cells - 1);

    for (int g = 0; g < n_genes; ++g) {
      double cov = 0.0;
      for (int c = 0; c < n_cells; ++c)
        cov += (expression(g, c) - means[g]) * (abs_probs(c, l) - mean_prob);
      cov /= (n_cells - 1);

      double denom = std::sqrt(vars[g] * var_prob);
      double r = denom > 1e-12 ? cov / denom : 0.0;
      if (r > 1.0) r = 1.0;
      if (r < -1.0) r = -1.0;
      corr(g, li) = r;

      double t_stat = n_cells > 3 ? r * std::sqrt((n_cells - 2.0) / (1.0 - r * r + 1e-15)) : 0.0;
      pval(g, li) = t_stat > 0 ? 1.0 / (1.0 + t_stat * t_stat) : 1.0;
    }
  }

  return List::create(
    _["correlation"] = corr,
    _["pval"] = pval,
    _["lineage_idx"] = lineage_idx
  );
}
