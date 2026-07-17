// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix estimate_ssgsea_scores_cpp(
    NumericMatrix ranked,
    IntegerMatrix sample_order,
    List gene_sets) {
  const int n_genes = ranked.nrow();
  const int n_samples = ranked.ncol();
  if (sample_order.nrow() != n_genes || sample_order.ncol() != n_samples) {
    stop("sample_order must have the same dimensions as ranked");
  }

  NumericMatrix scores(n_samples, gene_sets.size());
  std::fill(scores.begin(), scores.end(), NA_REAL);
  std::vector<unsigned char> in_set(n_genes, 0);

  for (int set_i = 0; set_i < gene_sets.size(); ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    std::fill(in_set.begin(), in_set.end(), 0);
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        in_set[gene] = 1;
      }
    }

    for (int sample_i = 0; sample_i < n_samples; ++sample_i) {
      double hit_weight_total = 0.0;
      int no_hit_total = 0;
      for (int position = 0; position < n_genes; ++position) {
        const int gene = sample_order(position, sample_i) - 1;
        if (gene < 0 || gene >= n_genes) {
          stop("sample_order contains an out-of-bounds gene index");
        }
        if (in_set[gene]) {
          hit_weight_total += R_pow(ranked(gene, sample_i), 0.25);
        } else {
          ++no_hit_total;
        }
      }
      if (hit_weight_total <= 0.0 || no_hit_total == 0) {
        continue;
      }

      double hit_reward = 0.0;
      double no_hit_penalty = 0.0;
      double score = 0.0;
      for (int position = 0; position < n_genes; ++position) {
        const int gene = sample_order(position, sample_i) - 1;
        if (in_set[gene]) {
          hit_reward += R_pow(ranked(gene, sample_i), 0.25) / hit_weight_total;
        } else {
          no_hit_penalty += 1.0 / static_cast<double>(no_hit_total);
        }
        score += hit_reward - no_hit_penalty;
      }
      scores(sample_i, set_i) = score;
    }
  }
  return scores;
}
