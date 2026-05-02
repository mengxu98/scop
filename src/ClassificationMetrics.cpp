// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

static double metric_choose2_cpp(double x) {
  return x * (x - 1.0) / 2.0;
}

static double metric_mean_na_rm(const NumericVector& x) {
  double total = 0.0;
  int n = 0;
  for (int i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (!NumericVector::is_na(value) && !R_IsNaN(value)) {
      total += value;
      ++n;
    }
  }
  if (n == 0) {
    return R_NaN;
  }
  return total / static_cast<double>(n);
}

static double metric_mean_na_rm_std(const std::vector<double>& x) {
  double total = 0.0;
  int n = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double value = x[i];
    if (!R_IsNA(value) && !R_IsNaN(value)) {
      total += value;
      ++n;
    }
  }
  if (n == 0) {
    return R_NaN;
  }
  return total / static_cast<double>(n);
}

static double metric_entropy_cpp(const std::vector<double>& counts, double n) {
  if (n <= 0.0) {
    return 0.0;
  }
  double out = 0.0;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] > 0.0) {
      const double p = counts[i] / n;
      out -= p * std::log(p);
    }
  }
  return out;
}

// [[Rcpp::export]]
List classification_metrics_cpp(
  CharacterVector predicted,
  CharacterVector truth,
  CharacterVector classes,
  double rare_threshold = 0.05
) {
  const int n_classes = classes.size();
  std::unordered_map<std::string, int> class_index;
  class_index.reserve(n_classes);
  for (int i = 0; i < n_classes; ++i) {
    class_index[as<std::string>(classes[i])] = i;
  }

  std::vector<double> tab(static_cast<std::size_t>(n_classes) * n_classes, 0.0);
  double n = 0.0;
  double correct = 0.0;

  const int n_obs = std::min(predicted.size(), truth.size());
  for (int i = 0; i < n_obs; ++i) {
    if (predicted[i] == NA_STRING || truth[i] == NA_STRING) {
      continue;
    }
    const std::string pred_i = as<std::string>(predicted[i]);
    const std::string truth_i = as<std::string>(truth[i]);
    std::unordered_map<std::string, int>::const_iterator pred_it = class_index.find(pred_i);
    std::unordered_map<std::string, int>::const_iterator truth_it = class_index.find(truth_i);
    if (pred_it == class_index.end() || truth_it == class_index.end()) {
      continue;
    }
    const int pred_idx = pred_it->second;
    const int truth_idx = truth_it->second;
    tab[static_cast<std::size_t>(pred_idx) * n_classes + truth_idx] += 1.0;
    n += 1.0;
    if (pred_idx == truth_idx) {
      correct += 1.0;
    }
  }

  if (n <= 0.0 || n_classes == 0) {
    NumericVector empty_num(0);
    IntegerVector empty_int(0);
    CharacterVector empty_char(0);
    DataFrame class_table = DataFrame::create(
      _["class"] = empty_char,
      _["precision"] = empty_num,
      _["recall"] = empty_num,
      _["f1"] = empty_num,
      _["support"] = empty_int,
      _["stringsAsFactors"] = false
    );
    return List::create(
      _["accuracy"] = NA_REAL,
      _["macro_f1"] = NA_REAL,
      _["purity"] = NA_REAL,
      _["nmi"] = NA_REAL,
      _["ari"] = NA_REAL,
      _["rare_recall"] = NA_REAL,
      _["class_table"] = class_table
    );
  }

  std::vector<double> row_sum(n_classes, 0.0);
  std::vector<double> col_sum(n_classes, 0.0);
  for (int row = 0; row < n_classes; ++row) {
    for (int col = 0; col < n_classes; ++col) {
      const double value = tab[static_cast<std::size_t>(row) * n_classes + col];
      row_sum[row] += value;
      col_sum[col] += value;
    }
  }

  NumericVector precision(n_classes);
  NumericVector recall(n_classes);
  NumericVector f1(n_classes);
  IntegerVector support(n_classes);

  double purity_num = 0.0;
  for (int cls = 0; cls < n_classes; ++cls) {
    const double tp = tab[static_cast<std::size_t>(cls) * n_classes + cls];
    const double fp = row_sum[cls] - tp;
    const double fn = col_sum[cls] - tp;

    precision[cls] = (tp + fp) == 0.0 ? NA_REAL : tp / (tp + fp);
    recall[cls] = (tp + fn) == 0.0 ? NA_REAL : tp / (tp + fn);
    if (
      NumericVector::is_na(precision[cls]) ||
      NumericVector::is_na(recall[cls]) ||
      (precision[cls] + recall[cls]) == 0.0
    ) {
      f1[cls] = NA_REAL;
    } else {
      f1[cls] = 2.0 * precision[cls] * recall[cls] / (precision[cls] + recall[cls]);
    }
    support[cls] = static_cast<int>(col_sum[cls]);

    double row_max = 0.0;
    for (int col = 0; col < n_classes; ++col) {
      row_max = std::max(row_max, tab[static_cast<std::size_t>(cls) * n_classes + col]);
    }
    purity_num += row_max;
  }

  double sum_nij = 0.0;
  for (std::size_t i = 0; i < tab.size(); ++i) {
    sum_nij += metric_choose2_cpp(tab[i]);
  }
  double sum_a = 0.0;
  double sum_b = 0.0;
  for (int cls = 0; cls < n_classes; ++cls) {
    sum_a += metric_choose2_cpp(row_sum[cls]);
    sum_b += metric_choose2_cpp(col_sum[cls]);
  }

  const double choose_n = metric_choose2_cpp(n);
  double ari = NA_REAL;
  if (choose_n > 0.0) {
    const double expected = sum_a * sum_b / choose_n;
    const double max_index = (sum_a + sum_b) / 2.0;
    if (max_index != expected) {
      ari = (sum_nij - expected) / (max_index - expected);
    }
  }

  double mi = 0.0;
  for (int row = 0; row < n_classes; ++row) {
    for (int col = 0; col < n_classes; ++col) {
      const double nij = tab[static_cast<std::size_t>(row) * n_classes + col];
      if (nij > 0.0 && row_sum[row] > 0.0 && col_sum[col] > 0.0) {
        const double pij = nij / n;
        const double pi = row_sum[row] / n;
        const double pj = col_sum[col] / n;
        mi += pij * std::log(pij / (pi * pj));
      }
    }
  }
  const double hx = metric_entropy_cpp(row_sum, n);
  const double hy = metric_entropy_cpp(col_sum, n);
  const double nmi = (hx + hy) == 0.0 ? NA_REAL : 2.0 * mi / (hx + hy);

  std::vector<double> rare_recalls;
  for (int cls = 0; cls < n_classes; ++cls) {
    const double support_prop = static_cast<double>(support[cls]) / n;
    if (support_prop <= rare_threshold) {
      rare_recalls.push_back(recall[cls]);
    }
  }
  const double rare_recall = rare_recalls.size() == 0 ? NA_REAL : metric_mean_na_rm_std(rare_recalls);

  DataFrame class_table = DataFrame::create(
    _["class"] = classes,
    _["precision"] = precision,
    _["recall"] = recall,
    _["f1"] = f1,
    _["support"] = support,
    _["stringsAsFactors"] = false
  );

  return List::create(
    _["accuracy"] = correct / n,
    _["macro_f1"] = metric_mean_na_rm(f1),
    _["purity"] = purity_num / n,
    _["nmi"] = nmi,
    _["ari"] = ari,
    _["rare_recall"] = rare_recall,
    _["class_table"] = class_table
  );
}
