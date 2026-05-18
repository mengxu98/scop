#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

struct ScoreIndex {
  double score;
  int index;
};

bool score_desc(const ScoreIndex& a, const ScoreIndex& b) {
  if (a.score == b.score) {
    return a.index < b.index;
  }
  return a.score > b.score;
}

void keep_top_scores(std::vector<ScoreIndex>& scores, int n_keep) {
  n_keep = std::max(0, std::min(n_keep, static_cast<int>(scores.size())));
  if (n_keep == 0) {
    scores.clear();
    return;
  }
  if (n_keep < static_cast<int>(scores.size())) {
    std::partial_sort(scores.begin(), scores.begin() + n_keep, scores.end(), score_desc);
    scores.resize(n_keep);
  } else {
    std::sort(scores.begin(), scores.end(), score_desc);
  }
}

std::vector<int> unique_sorted_indices(std::vector<int> x) {
  std::sort(x.begin(), x.end());
  x.erase(std::unique(x.begin(), x.end()), x.end());
  return x;
}

struct DgcMatrixView {
  int n_gene;
  int n_cell;
  IntegerVector row_idx;
  IntegerVector col_ptr;
  NumericVector values;
};

DgcMatrixView dgc_view(S4 x) {
  IntegerVector dims = x.slot("Dim");
  if (dims.size() != 2) {
    stop("sparse matrix must have two dimensions");
  }
  DgcMatrixView out;
  out.n_gene = dims[0];
  out.n_cell = dims[1];
  out.row_idx = x.slot("i");
  out.col_ptr = x.slot("p");
  out.values = x.slot("x");
  return out;
}

std::vector<int> selected_gene_indices(
    const NumericMatrix& ref,
    const IntegerVector& labels,
    int n_labels,
    int n_top,
    int additional_per_label) {
  const int n_gene = ref.nrow();
  const int n_cell = ref.ncol();
  std::vector<int> label_n(n_labels, 0);
  for (int cell = 0; cell < n_cell; ++cell) {
    const int lab = labels[cell] - 1;
    if (lab >= 0 && lab < n_labels) {
      label_n[lab]++;
    }
  }

  std::vector<double> type_mean(n_gene * n_labels, 0.0);
  std::vector<ScoreIndex> scores;
  scores.reserve(n_gene);
  for (int gene = 0; gene < n_gene; ++gene) {
    double mean_sum = 0.0;
    double log_sum = 0.0;
    for (int cell = 0; cell < n_cell; ++cell) {
      const int lab = labels[cell] - 1;
      if (lab < 0 || lab >= n_labels) {
        continue;
      }
      double value = ref(gene, cell);
      if (R_finite(value) && value > 0.0) {
        type_mean[gene * n_labels + lab] += value;
      }
    }
    for (int lab = 0; lab < n_labels; ++lab) {
      double mean = 0.0;
      if (label_n[lab] > 0) {
        mean = type_mean[gene * n_labels + lab] / label_n[lab];
      }
      type_mean[gene * n_labels + lab] = mean;
      mean_sum += mean;
      log_sum += std::log(mean + 1.0);
    }
    const double score = std::log(mean_sum / n_labels + 1.0) * n_labels - log_sum;
    scores.push_back({score, gene});
  }

  std::vector<int> selected;
  const int n_keep = std::min(std::max(n_top, 0), n_gene);
  selected.reserve(n_keep + std::max(0, additional_per_label) * n_labels);
  keep_top_scores(scores, n_keep);
  for (int i = 0; i < static_cast<int>(scores.size()); ++i) {
    selected.push_back(scores[i].index);
  }

  if (additional_per_label > 0) {
    for (int lab = 0; lab < n_labels; ++lab) {
      std::vector<ScoreIndex> means;
      means.reserve(n_gene);
      for (int gene = 0; gene < n_gene; ++gene) {
        means.push_back({type_mean[gene * n_labels + lab], gene});
      }
      const int add_keep = std::min(additional_per_label, n_gene);
      keep_top_scores(means, add_keep);
      for (int i = 0; i < static_cast<int>(means.size()); ++i) {
        selected.push_back(means[i].index);
      }
    }
  }

  return unique_sorted_indices(selected);
}

std::vector<int> selected_gene_indices_sparse(
    const DgcMatrixView& ref,
    const IntegerVector& labels,
    int n_labels,
    int n_top,
    int additional_per_label) {
  const int n_gene = ref.n_gene;
  const int n_cell = ref.n_cell;
  std::vector<int> label_n(n_labels, 0);
  for (int cell = 0; cell < n_cell; ++cell) {
    const int lab = labels[cell] - 1;
    if (lab >= 0 && lab < n_labels) {
      label_n[lab]++;
    }
  }

  std::vector<double> type_mean(n_gene * n_labels, 0.0);
  for (int cell = 0; cell < n_cell; ++cell) {
    const int lab = labels[cell] - 1;
    if (lab < 0 || lab >= n_labels) {
      continue;
    }
    for (int ptr = ref.col_ptr[cell]; ptr < ref.col_ptr[cell + 1]; ++ptr) {
      const int gene = ref.row_idx[ptr];
      const double value = ref.values[ptr];
      if (gene >= 0 && gene < n_gene && R_finite(value) && value > 0.0) {
        type_mean[gene * n_labels + lab] += value;
      }
    }
  }

  std::vector<ScoreIndex> scores;
  scores.reserve(n_gene);
  for (int gene = 0; gene < n_gene; ++gene) {
    double mean_sum = 0.0;
    double log_sum = 0.0;
    for (int lab = 0; lab < n_labels; ++lab) {
      double mean = 0.0;
      if (label_n[lab] > 0) {
        mean = type_mean[gene * n_labels + lab] / label_n[lab];
      }
      type_mean[gene * n_labels + lab] = mean;
      mean_sum += mean;
      log_sum += std::log(mean + 1.0);
    }
    const double score = std::log(mean_sum / n_labels + 1.0) * n_labels - log_sum;
    scores.push_back({score, gene});
  }

  std::vector<int> selected;
  const int n_keep = std::min(std::max(n_top, 0), n_gene);
  selected.reserve(n_keep + std::max(0, additional_per_label) * n_labels);
  keep_top_scores(scores, n_keep);
  for (int i = 0; i < static_cast<int>(scores.size()); ++i) {
    selected.push_back(scores[i].index);
  }

  if (additional_per_label > 0) {
    for (int lab = 0; lab < n_labels; ++lab) {
      std::vector<ScoreIndex> means;
      means.reserve(n_gene);
      for (int gene = 0; gene < n_gene; ++gene) {
        means.push_back({type_mean[gene * n_labels + lab], gene});
      }
      keep_top_scores(means, std::min(additional_per_label, n_gene));
      for (int i = 0; i < static_cast<int>(means.size()); ++i) {
        selected.push_back(means[i].index);
      }
    }
  }

  return unique_sorted_indices(selected);
}

NumericMatrix build_log_probability_core(
    const NumericMatrix& ref,
    const IntegerVector& labels,
    const std::vector<int>& genes,
    int n_labels) {
  const int n_cell = ref.ncol();
  const int n_gene = genes.size();
  NumericMatrix core(n_gene, n_labels);
  std::vector<double> sums(n_labels, 0.0);
  const double log2 = std::log(2.0);

  for (int g = 0; g < n_gene; ++g) {
    const int gene = genes[g];
    for (int cell = 0; cell < n_cell; ++cell) {
      const int lab = labels[cell] - 1;
      if (lab < 0 || lab >= n_labels) {
        continue;
      }
      double value = ref(gene, cell);
      if (!R_finite(value) || value < 0.0) {
        value = 0.0;
      }
      const double log_value = std::log(value + 1.0) / log2;
      core(g, lab) += log_value;
      sums[lab] += log_value;
    }
  }

  for (int lab = 0; lab < n_labels; ++lab) {
    const double denom = sums[lab] + n_gene;
    for (int g = 0; g < n_gene; ++g) {
      core(g, lab) = std::log(core(g, lab) + 1.0) - std::log(denom);
    }
  }
  return core;
}

NumericMatrix build_log_probability_core_sparse(
    const DgcMatrixView& ref,
    const IntegerVector& labels,
    const std::vector<int>& genes,
    int n_labels) {
  const int n_gene = genes.size();
  NumericMatrix core(n_gene, n_labels);
  std::vector<double> sums(n_labels, 0.0);
  std::vector<int> gene_to_selected(ref.n_gene, -1);
  const double log2 = std::log(2.0);

  for (int g = 0; g < n_gene; ++g) {
    gene_to_selected[genes[g]] = g;
  }

  for (int cell = 0; cell < ref.n_cell; ++cell) {
    const int lab = labels[cell] - 1;
    if (lab < 0 || lab >= n_labels) {
      continue;
    }
    for (int ptr = ref.col_ptr[cell]; ptr < ref.col_ptr[cell + 1]; ++ptr) {
      const int gene = ref.row_idx[ptr];
      if (gene < 0 || gene >= ref.n_gene) {
        continue;
      }
      const int g = gene_to_selected[gene];
      if (g < 0) {
        continue;
      }
      double value = ref.values[ptr];
      if (!R_finite(value) || value < 0.0) {
        value = 0.0;
      }
      const double log_value = std::log(value + 1.0) / log2;
      core(g, lab) += log_value;
      sums[lab] += log_value;
    }
  }

  for (int lab = 0; lab < n_labels; ++lab) {
    const double denom = sums[lab] + n_gene;
    for (int g = 0; g < n_gene; ++g) {
      core(g, lab) = std::log(core(g, lab) + 1.0) - std::log(denom);
    }
  }
  return core;
}

NumericMatrix predict_probabilities(
    const NumericMatrix& query,
    const NumericMatrix& core,
    const std::vector<int>& genes) {
  const int n_query = query.ncol();
  const int n_labels = core.ncol();
  const int n_gene = genes.size();
  NumericMatrix prob(n_query, n_labels);
  const double log2 = std::log(2.0);

  for (int cell = 0; cell < n_query; ++cell) {
    double max_score = -std::numeric_limits<double>::infinity();
    for (int lab = 0; lab < n_labels; ++lab) {
      double score = 0.0;
      for (int g = 0; g < n_gene; ++g) {
        double value = query(genes[g], cell);
        if (!R_finite(value) || value < 0.0) {
          value = 0.0;
        }
        score += (std::log(value + 1.0) / log2) * core(g, lab);
      }
      prob(cell, lab) = score;
      if (score > max_score) {
        max_score = score;
      }
    }

    double sum_exp = 0.0;
    for (int lab = 0; lab < n_labels; ++lab) {
      const double val = std::exp(prob(cell, lab) - max_score);
      prob(cell, lab) = val;
      sum_exp += val;
    }
    if (sum_exp > 0.0) {
      for (int lab = 0; lab < n_labels; ++lab) {
        prob(cell, lab) /= sum_exp;
      }
    }
  }
  return prob;
}

NumericMatrix predict_probabilities_sparse(
    const DgcMatrixView& query,
    const NumericMatrix& core,
    const std::vector<int>& genes) {
  const int n_query = query.n_cell;
  const int n_labels = core.ncol();
  NumericMatrix prob(n_query, n_labels);
  std::vector<int> gene_to_selected(query.n_gene, -1);
  std::vector<double> scores(n_labels, 0.0);
  const double log2 = std::log(2.0);

  for (int g = 0; g < static_cast<int>(genes.size()); ++g) {
    gene_to_selected[genes[g]] = g;
  }

  for (int cell = 0; cell < n_query; ++cell) {
    std::fill(scores.begin(), scores.end(), 0.0);
    for (int ptr = query.col_ptr[cell]; ptr < query.col_ptr[cell + 1]; ++ptr) {
      const int gene = query.row_idx[ptr];
      if (gene < 0 || gene >= query.n_gene) {
        continue;
      }
      const int g = gene_to_selected[gene];
      if (g < 0) {
        continue;
      }
      double value = query.values[ptr];
      if (!R_finite(value) || value < 0.0) {
        value = 0.0;
      }
      const double log_value = std::log(value + 1.0) / log2;
      for (int lab = 0; lab < n_labels; ++lab) {
        scores[lab] += log_value * core(g, lab);
      }
    }

    double max_score = -std::numeric_limits<double>::infinity();
    for (int lab = 0; lab < n_labels; ++lab) {
      if (scores[lab] > max_score) {
        max_score = scores[lab];
      }
    }
    double sum_exp = 0.0;
    for (int lab = 0; lab < n_labels; ++lab) {
      const double val = std::exp(scores[lab] - max_score);
      prob(cell, lab) = val;
      sum_exp += val;
    }
    if (sum_exp > 0.0) {
      for (int lab = 0; lab < n_labels; ++lab) {
        prob(cell, lab) /= sum_exp;
      }
    }
  }
  return prob;
}

IntegerVector max_label_indices(const NumericMatrix& prob) {
  const int n_query = prob.nrow();
  const int n_labels = prob.ncol();
  IntegerVector out(n_query);
  for (int cell = 0; cell < n_query; ++cell) {
    int best = 0;
    double best_value = prob(cell, 0);
    for (int lab = 1; lab < n_labels; ++lab) {
      if (prob(cell, lab) > best_value) {
        best_value = prob(cell, lab);
        best = lab;
      }
    }
    out[cell] = best + 1;
  }
  return out;
}

}  // namespace

// [[Rcpp::export]]
List scibet_fit_predict(
    NumericMatrix ref,
    NumericMatrix query,
    IntegerVector labels,
    int n_labels,
    int n_top,
    int additional_per_label = 0) {
  if (ref.nrow() != query.nrow()) {
    stop("ref and query must have the same number of rows");
  }
  if (ref.ncol() != labels.size()) {
    stop("labels must have one value per reference cell");
  }
  if (n_labels < 2) {
    stop("at least two labels are required");
  }
  if (n_top < 1) {
    stop("n_top must be positive");
  }

  std::vector<int> genes = selected_gene_indices(
    ref,
    labels,
    n_labels,
    n_top,
    additional_per_label
  );
  if (genes.empty()) {
    stop("no SciBet features were selected");
  }

  NumericMatrix core = build_log_probability_core(ref, labels, genes, n_labels);
  NumericMatrix prob = predict_probabilities(query, core, genes);
  IntegerVector predicted = max_label_indices(prob);

  IntegerVector genes_out(genes.size());
  for (int i = 0; i < static_cast<int>(genes.size()); ++i) {
    genes_out[i] = genes[i] + 1;
  }

  return List::create(
    _["feature_index"] = genes_out,
    _["core"] = core,
    _["probabilities"] = prob,
    _["predicted_index"] = predicted
  );
}

// [[Rcpp::export]]
List scibet_fit_predict_sparse(
    S4 ref,
    S4 query,
    IntegerVector labels,
    int n_labels,
    int n_top,
    int additional_per_label = 0) {
  DgcMatrixView ref_view = dgc_view(ref);
  DgcMatrixView query_view = dgc_view(query);
  if (ref_view.n_gene != query_view.n_gene) {
    stop("ref and query must have the same number of rows");
  }
  if (ref_view.n_cell != labels.size()) {
    stop("labels must have one value per reference cell");
  }
  if (n_labels < 2) {
    stop("at least two labels are required");
  }
  if (n_top < 1) {
    stop("n_top must be positive");
  }

  std::vector<int> genes = selected_gene_indices_sparse(
    ref_view,
    labels,
    n_labels,
    n_top,
    additional_per_label
  );
  if (genes.empty()) {
    stop("no SciBet features were selected");
  }

  NumericMatrix core = build_log_probability_core_sparse(ref_view, labels, genes, n_labels);
  NumericMatrix prob = predict_probabilities_sparse(query_view, core, genes);
  IntegerVector predicted = max_label_indices(prob);

  IntegerVector genes_out(genes.size());
  for (int i = 0; i < static_cast<int>(genes.size()); ++i) {
    genes_out[i] = genes[i] + 1;
  }

  return List::create(
    _["feature_index"] = genes_out,
    _["core"] = core,
    _["probabilities"] = prob,
    _["predicted_index"] = predicted
  );
}

// [[Rcpp::export]]
NumericMatrix scibet_predict(NumericMatrix query, NumericMatrix core, IntegerVector feature_index) {
  std::vector<int> genes(feature_index.size());
  for (int i = 0; i < feature_index.size(); ++i) {
    genes[i] = feature_index[i] - 1;
  }
  return predict_probabilities(query, core, genes);
}
