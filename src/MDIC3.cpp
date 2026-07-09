#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

namespace {

arma::mat mdic3_log10_positive(const arma::mat& x) {
  arma::mat out(x.n_rows, x.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      const double value = x(i, j);
      if (std::isfinite(value) && value > 0.0) {
        const double logged = std::log10(value);
        if (logged > 0.0) out(i, j) = logged;
      }
    }
  }
  return out;
}

arma::mat mdic3_cellular_raw(const arma::mat& m) {
  arma::mat out(m.n_rows, m.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < m.n_rows; ++i) {
    for (arma::uword j = 0; j < m.n_cols; ++j) {
      const double value = m(i, j);
      if (!std::isfinite(value) || value == 0.0) continue;
      if (value > 0.0) {
        out(i, j) += value;
      } else {
        out(j, i) += std::abs(value);
      }
    }
  }
  return out;
}

arma::mat mdic3_type_raw(const arma::mat& cellular, const IntegerVector& group) {
  const int n_cells = static_cast<int>(cellular.n_rows);
  int n_groups = 0;
  for (int i = 0; i < group.size(); ++i) {
    if (IntegerVector::is_na(group[i]) || group[i] < 1) {
      stop("mdic3_score_cpp: group ids must be positive integers without NA");
    }
    if (group[i] > n_groups) n_groups = group[i];
  }
  arma::mat out(n_groups, n_groups, arma::fill::zeros);
  for (int i = 0; i < n_cells; ++i) {
    const int gi = group[i] - 1;
    for (int j = 0; j < n_cells; ++j) {
      const int gj = group[j] - 1;
      const double value = cellular(i, j);
      if (!std::isfinite(value) || value == 0.0) continue;
      if (gi == gj) {
        out(gi, gj) += std::abs(value);
      } else if (value > 0.0) {
        out(gi, gj) += value;
      } else {
        out(gj, gi) += std::abs(value);
      }
    }
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
List mdic3_score_cpp(NumericMatrix expression, NumericMatrix grn, IntegerVector group) {
  const int n_genes = expression.nrow();
  const int n_cells = expression.ncol();
  if (grn.nrow() != n_genes || grn.ncol() != n_genes) {
    stop("mdic3_score_cpp: grn must be a square gene-by-gene matrix aligned to expression rows");
  }
  if (group.size() != n_cells) {
    stop("mdic3_score_cpp: group length must match expression columns");
  }
  if (n_genes < 1 || n_cells < 2) {
    stop("mdic3_score_cpp: expression must contain at least one gene and two cells");
  }

  arma::mat aa = as<arma::mat>(expression);
  arma::mat grn_mat = as<arma::mat>(grn);
  aa.replace(arma::datum::nan, 0.0);
  grn_mat.replace(arma::datum::nan, 0.0);

  arma::mat u;
  arma::vec s;
  arma::mat v;
  if (!arma::svd_econ(u, s, v, aa, "both", "std")) {
    stop("mdic3_score_cpp: SVD failed");
  }
  arma::mat s_rect(n_genes, n_cells, arma::fill::zeros);
  const arma::uword diag_n = std::min<arma::uword>(
    static_cast<arma::uword>(s.n_elem),
    std::min<arma::uword>(s_rect.n_rows, s_rect.n_cols)
  );
  for (arma::uword i = 0; i < diag_n; ++i) {
    s_rect(i, i) = s(i);
  }

  const arma::mat t_mat = grn_mat * s_rect;
  const arma::mat m = arma::pinv(t_mat) * aa;
  const arma::mat cellular = mdic3_cellular_raw(m);
  const arma::mat cellular_log = mdic3_log10_positive(cellular);
  const arma::mat type_raw = mdic3_type_raw(cellular, group);
  const arma::mat type_log = mdic3_log10_positive(type_raw);

  return List::create(
    _["mdic3_matrix"] = wrap(m),
    _["cellular_communication"] = wrap(cellular),
    _["cellular_communication_log"] = wrap(cellular_log),
    _["celltype_communication_raw"] = wrap(type_raw),
    _["celltype_communication"] = wrap(type_log)
  );
}
