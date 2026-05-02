// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

struct SparseWilcoxEntry {
  double value;
  bool group1;
};

struct SparseWilcoxGroup {
  double value;
  int n1;
  int n2;
};

static bool sparse_wilcox_entry_less(
  const SparseWilcoxEntry& a,
  const SparseWilcoxEntry& b
) {
  if (a.value < b.value) {
    return true;
  }
  if (a.value > b.value) {
    return false;
  }
  return a.group1 && !b.group1;
}

static bool sparse_wilcox_group_less(
  const SparseWilcoxGroup& a,
  const SparseWilcoxGroup& b
) {
  return a.value < b.value;
}

static double sparse_wilcox_two_sided_p_from_groups(
  std::vector<SparseWilcoxGroup>& groups,
  int n1,
  int n2
) {
  const int n = n1 + n2;
  if (n <= 1 || n1 == 0 || n2 == 0) {
    return 1.0;
  }

  std::sort(groups.begin(), groups.end(), sparse_wilcox_group_less);

  double rank_sum1 = 0.0;
  double tie_adjust_num = 0.0;
  int seen = 0;
  int i = 0;
  const int n_groups = static_cast<int>(groups.size());
  while (i < n_groups) {
    int j = i + 1;
    int tie_n1 = groups[i].n1;
    int tie_n2 = groups[i].n2;
    while (j < n_groups && groups[j].value == groups[i].value) {
      tie_n1 += groups[j].n1;
      tie_n2 += groups[j].n2;
      ++j;
    }

    const int tie_n = tie_n1 + tie_n2;
    if (tie_n > 0) {
      const double avg_rank =
        (static_cast<double>(seen + 1) + static_cast<double>(seen + tie_n)) / 2.0;
      rank_sum1 += avg_rank * static_cast<double>(tie_n1);
      if (tie_n > 1) {
        tie_adjust_num += static_cast<double>(tie_n) *
          static_cast<double>(tie_n + 1) *
          static_cast<double>(tie_n - 1);
      }
      seen += tie_n;
    }
    i = j;
  }

  if (seen != n) {
    return NA_REAL;
  }

  const double n1d = static_cast<double>(n1);
  const double n2d = static_cast<double>(n2);
  const double nd = static_cast<double>(n);
  const double u = n1d * n2d + n1d * (n1d + 1.0) / 2.0 - rank_sum1;
  const double mu = n1d * n2d / 2.0;
  double sigma2 = n1d * n2d * (nd + 1.0) / 12.0;
  if (tie_adjust_num > 0.0 && n > 1) {
    const double adjustment = tie_adjust_num / (nd * (nd + 1.0) * (nd - 1.0));
    sigma2 *= (1.0 - adjustment);
  }
  if (!R_finite(sigma2) || sigma2 <= 0.0) {
    return 1.0;
  }

  const double sigma = std::sqrt(sigma2);
  const double z_lower_tail = (u + 0.5 - mu) / sigma;
  const double z_upper_tail = (u - 0.5 - mu) / sigma;
  const double less = R::pnorm(z_upper_tail, 0.0, 1.0, false, false);
  const double greater = R::pnorm(z_lower_tail, 0.0, 1.0, true, false);
  double p = 2.0 * std::min(less, greater);
  if (!R_finite(p)) {
    return NA_REAL;
  }
  if (p > 1.0) {
    p = 1.0;
  }
  if (p < 0.0) {
    p = 0.0;
  }
  return p;
}

static double sparse_wilcox_two_sided_p(std::vector<SparseWilcoxEntry>& entries) {
  const int n = static_cast<int>(entries.size());
  if (n <= 1) {
    return 1.0;
  }

  int n1 = 0;
  for (int i = 0; i < n; ++i) {
    if (entries[i].group1) {
      ++n1;
    }
  }
  const int n2 = n - n1;
  if (n1 == 0 || n2 == 0) {
    return 1.0;
  }

  std::sort(entries.begin(), entries.end(), sparse_wilcox_entry_less);

  double rank_sum1 = 0.0;
  double tie_adjust_num = 0.0;
  int i = 0;
  while (i < n) {
    int j = i + 1;
    while (j < n && entries[j].value == entries[i].value) {
      ++j;
    }

    const int tie_n = j - i;
    const double avg_rank = (static_cast<double>(i + 1) + static_cast<double>(j)) / 2.0;
    int tie_n1 = 0;
    for (int k = i; k < j; ++k) {
      if (entries[k].group1) {
        ++tie_n1;
      }
    }
    rank_sum1 += avg_rank * static_cast<double>(tie_n1);
    if (tie_n > 1) {
      tie_adjust_num += static_cast<double>(tie_n) *
        static_cast<double>(tie_n + 1) *
        static_cast<double>(tie_n - 1);
    }
    i = j;
  }

  const double n1d = static_cast<double>(n1);
  const double n2d = static_cast<double>(n2);
  const double nd = static_cast<double>(n);
  const double u = n1d * n2d + n1d * (n1d + 1.0) / 2.0 - rank_sum1;
  const double mu = n1d * n2d / 2.0;
  double sigma2 = n1d * n2d * (nd + 1.0) / 12.0;
  if (tie_adjust_num > 0.0 && n > 1) {
    const double adjustment = tie_adjust_num / (nd * (nd + 1.0) * (nd - 1.0));
    sigma2 *= (1.0 - adjustment);
  }
  if (!R_finite(sigma2) || sigma2 <= 0.0) {
    return 1.0;
  }

  const double sigma = std::sqrt(sigma2);
  const double z_lower_tail = (u + 0.5 - mu) / sigma;
  const double z_upper_tail = (u - 0.5 - mu) / sigma;
  const double less = R::pnorm(z_upper_tail, 0.0, 1.0, false, false);
  const double greater = R::pnorm(z_lower_tail, 0.0, 1.0, true, false);
  double p = 2.0 * std::min(less, greater);
  if (!R_finite(p)) {
    return NA_REAL;
  }
  if (p > 1.0) {
    p = 1.0;
  }
  if (p < 0.0) {
    p = 0.0;
  }
  return p;
}

static double sparse_wilcox_two_sided_p_all_cells(
  std::vector<SparseWilcoxEntry>& entries,
  int n1,
  int n2
) {
  int nonzero_n1 = 0;
  int nonzero_n2 = 0;
  for (int i = 0; i < static_cast<int>(entries.size()); ++i) {
    if (entries[i].group1) {
      ++nonzero_n1;
    } else {
      ++nonzero_n2;
    }
  }

  std::sort(entries.begin(), entries.end(), sparse_wilcox_entry_less);

  std::vector<SparseWilcoxGroup> groups;
  const int zero_n1 = n1 - nonzero_n1;
  const int zero_n2 = n2 - nonzero_n2;
  if (zero_n1 + zero_n2 > 0) {
    groups.push_back(SparseWilcoxGroup{0.0, zero_n1, zero_n2});
  }

  int i = 0;
  const int n_entries = static_cast<int>(entries.size());
  while (i < n_entries) {
    int j = i + 1;
    int tie_n1 = entries[i].group1 ? 1 : 0;
    int tie_n2 = entries[i].group1 ? 0 : 1;
    while (j < n_entries && entries[j].value == entries[i].value) {
      if (entries[j].group1) {
        ++tie_n1;
      } else {
        ++tie_n2;
      }
      ++j;
    }
    groups.push_back(SparseWilcoxGroup{entries[i].value, tie_n1, tie_n2});
    i = j;
  }

  return sparse_wilcox_two_sided_p_from_groups(groups, n1, n2);
}

// [[Rcpp::export]]
NumericVector wilcox_rank_sum_sparse_cpp(
  S4 mat,
  int n_group1,
  double min_expression = 0.0
) {
  IntegerVector dims = mat.slot("Dim");
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  if (n_group1 < 0 || n_group1 > n_cols) {
    stop("n_group1 must be between 0 and ncol(mat)");
  }

  IntegerVector row_idx = mat.slot("i");
  IntegerVector col_ptr = mat.slot("p");
  NumericVector values = mat.slot("x");

  std::vector<std::vector<SparseWilcoxEntry> > rows(n_rows);
  for (int col = 0; col < n_cols; ++col) {
    const bool is_group1 = col < n_group1;
    for (int ptr = col_ptr[col]; ptr < col_ptr[col + 1]; ++ptr) {
      const double value = values[ptr];
      if (!R_finite(value) || value <= min_expression) {
        continue;
      }
      rows[row_idx[ptr]].push_back(SparseWilcoxEntry{value, is_group1});
    }
  }

  NumericVector p_values(n_rows);
  for (int row = 0; row < n_rows; ++row) {
    p_values[row] = sparse_wilcox_two_sided_p(rows[row]);
  }
  return p_values;
}

// [[Rcpp::export]]
NumericVector wilcox_rank_sum_sparse_all_cells_cpp(
  S4 mat,
  int n_group1
) {
  IntegerVector dims = mat.slot("Dim");
  const int n_rows = dims[0];
  const int n_cols = dims[1];
  if (n_group1 < 0 || n_group1 > n_cols) {
    stop("n_group1 must be between 0 and ncol(mat)");
  }

  IntegerVector row_idx = mat.slot("i");
  IntegerVector col_ptr = mat.slot("p");
  NumericVector values = mat.slot("x");

  std::vector<std::vector<SparseWilcoxEntry> > rows(n_rows);
  for (int col = 0; col < n_cols; ++col) {
    const bool is_group1 = col < n_group1;
    for (int ptr = col_ptr[col]; ptr < col_ptr[col + 1]; ++ptr) {
      const double value = values[ptr];
      if (!R_finite(value) || value == 0.0) {
        continue;
      }
      rows[row_idx[ptr]].push_back(SparseWilcoxEntry{value, is_group1});
    }
  }

  NumericVector p_values(n_rows);
  const int n_group2 = n_cols - n_group1;
  for (int row = 0; row < n_rows; ++row) {
    p_values[row] = sparse_wilcox_two_sided_p_all_cells(
      rows[row],
      n_group1,
      n_group2
    );
  }
  return p_values;
}
