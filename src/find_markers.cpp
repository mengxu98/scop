// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

struct RowSparseBlock {
  std::vector<int> row_ptr;
  std::vector<int> col_index;
  std::vector<double> value;
  int cells;
  int features;
};

static RowSparseBlock make_row_blocks(S4 matrix) {
  NumericVector values = matrix.slot("x");
  IntegerVector colptr = matrix.slot("p");
  IntegerVector rows = matrix.slot("i");
  IntegerVector dims = matrix.slot("Dim");

  RowSparseBlock block;
  block.features = dims[0];
  block.cells = dims[1];
  block.row_ptr.assign(block.features + 1, 0);
  block.col_index.resize(values.size());
  block.value.resize(values.size());

  for (int pos = 0; pos < rows.size(); ++pos) {
    ++block.row_ptr[rows[pos] + 1];
  }
  for (int feature = 0; feature < block.features; ++feature) {
    block.row_ptr[feature + 1] += block.row_ptr[feature];
  }

  std::vector<int> cursor = block.row_ptr;
  for (int cell = 0; cell < block.cells; ++cell) {
    for (int pos = colptr[cell]; pos < colptr[cell + 1]; ++pos) {
      int feature = rows[pos];
      int target = cursor[feature]++;
      block.col_index[target] = cell;
      block.value[target] = values[pos];
    }
  }
  return block;
}

struct MarkerSummaryWorker : public Worker {
  const std::vector<int>& row_ptr;
  const std::vector<int>& col_index;
  const std::vector<double>& value;
  const RVector<int> groups;
  const RVector<int> group_size;
  RMatrix<double> p_value;
  RMatrix<double> sum_value;
  RMatrix<double> detected;
  int cells;
  int group_count;
  double tie_scale;
  double rank_center;

  MarkerSummaryWorker(const RowSparseBlock& block,
                      IntegerVector groups,
                      IntegerVector group_size,
                      NumericMatrix p_value,
                      NumericMatrix sum_value,
                      NumericMatrix detected)
      : row_ptr(block.row_ptr), col_index(block.col_index), value(block.value),
        groups(groups),
        group_size(group_size), p_value(p_value), sum_value(sum_value),
        detected(detected), cells(block.cells), group_count(group_size.size()) {
    const double n = static_cast<double>(cells);
    tie_scale = 1.0 / (12.0 * (n * n - n));
    rank_center = n * n * n - n;
  }

  void operator()(std::size_t first_feature, std::size_t last_feature) {
    std::vector<int> counts(group_count);
    std::vector<double> sums(group_count);
    std::vector<double> ranks(group_count);
    std::vector<int> order;

    for (std::size_t feature = first_feature; feature < last_feature; ++feature) {
      std::fill(counts.begin(), counts.end(), 0);
      std::fill(sums.begin(), sums.end(), 0.0);
      std::fill(ranks.begin(), ranks.end(), 0.0);

      const int begin = row_ptr[feature];
      const int end = row_ptr[feature + 1];
      const int nonzero = end - begin;
      const int zero = cells - nonzero;
      order.resize(nonzero);

      for (int k = 0; k < nonzero; ++k) {
        const int pos = begin + k;
        order[k] = pos;
        int group = groups[col_index[pos]] - 1;
        if (group >= 0 && group < group_count) {
          ++counts[group];
          sums[group] += std::expm1(value[pos]);
        }
      }

      const double zero_rank = (zero + 1.0) * 0.5;
      for (int group = 0; group < group_count; ++group) {
        detected(feature, group) = counts[group];
        sum_value(feature, group) = sums[group];
        ranks[group] = (static_cast<double>(group_size[group]) - counts[group]) *
          zero_rank;
      }

      std::sort(order.begin(), order.end(),
                [&](int a, int b) {
                  if (value[a] == value[b]) {
                    return col_index[a] < col_index[b];
                  }
                  return value[a] < value[b];
                });

      double tie_penalty = zero > 1 ?
        static_cast<double>(zero) * zero * zero - zero : 0.0;
      int run_begin = 0;
      int rank_base = zero + 1;
      while (run_begin < nonzero) {
        int run_end = run_begin + 1;
        while (run_end < nonzero &&
               value[order[run_end]] == value[order[run_begin]]) {
          ++run_end;
        }
        const int width = run_end - run_begin;
        const double rank = rank_base + (width - 1.0) * 0.5;
        for (int pos = run_begin; pos < run_end; ++pos) {
          int group = groups[col_index[order[pos]]] - 1;
          if (group >= 0 && group < group_count) {
            ranks[group] += rank;
          }
        }
        if (width > 1) {
          tie_penalty += static_cast<double>(width) * width * width - width;
        }
        rank_base += width;
        run_begin = run_end;
      }

      const double variance_component = (rank_center - tie_penalty) * tie_scale;
      for (int group = 0; group < group_count; ++group) {
        const double n1 = group_size[group];
        const double n2 = cells - n1;
        const double pairs = n1 * n2;
        double z = ranks[group] - n1 * (n1 + 1.0) * 0.5 - pairs * 0.5;
        if (z > 0.0) {
          z -= 0.5;
        } else if (z < 0.0) {
          z += 0.5;
        }
        const double sigma = std::sqrt(pairs * variance_component);
        p_value(feature, group) = 2.0 *
          R::pnorm5(-std::abs(z / sigma), 0.0, 1.0, 1, 0);
      }
    }
  }
};

// [[Rcpp::export]]
List parallel_all_in_one_dgc(SEXP x_sexp,
                             IntegerVector groups,
                             IntegerVector group_sizes) {
  RowSparseBlock rows = make_row_blocks(S4(x_sexp));
  const int groups_n = group_sizes.size();
  NumericMatrix p(rows.features, groups_n);
  NumericMatrix sums(rows.features, groups_n);
  NumericMatrix counts(rows.features, groups_n);

  MarkerSummaryWorker worker(rows, groups, group_sizes, p, sums, counts);
  parallelFor(0, rows.features, worker);

  return List::create(
    Named("pval_by_group") = p,
    Named("sum_by_group") = sums,
    Named("detected_by_group") = counts
  );
}
