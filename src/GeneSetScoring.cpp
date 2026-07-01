// [[Rcpp::depends(RcppArmadillo, cli)]]
#include <RcppArmadillo.h>
#include "log_message.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace Rcpp;

struct AucEntry {
  int gene;
  double value;
  double tiebreak;
};

// Deterministic gene-index hash for tie-breaking.
// Replaces unif_rand() so that the same gene always gets the same tiebreaker
// regardless of strategy (sparse / topk / full), making them produce identical
// rankings for the same input.  Also eliminates RNG-call-order sensitivity.
static inline double gene_hash_tiebreak(int gene) {
  // SplitMix64-style hash, returns a deterministic value in [0, 1)
  unsigned long long x = static_cast<unsigned long long>(gene) + 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return static_cast<double>(x & 0x7fffffffffffffffULL) / static_cast<double>(0x8000000000000000ULL);
}

static bool aucell_entry_before(const AucEntry& a, const AucEntry& b) {
  if (a.value > b.value) {
    return true;
  }
  if (a.value < b.value) {
    return false;
  }
  return a.tiebreak < b.tiebreak;
}

static double aucell_max_auc(int n_genes, int auc_threshold, bool norm_auc) {
  if (!norm_auc) {
    return static_cast<double>(auc_threshold) * static_cast<double>(n_genes);
  }

  const int m = std::min(n_genes, auc_threshold - 1);
  if (m <= 0) {
    return 0.0;
  }

  double out = 0.0;
  for (int i = 1; i < m; ++i) {
    out += static_cast<double>(i);
  }
  out += static_cast<double>(auc_threshold - m) * static_cast<double>(m);
  return out;
}

static double aucell_auc_from_ranks(
  const std::vector<int>& ranks,
  int auc_threshold,
  double max_auc
) {
  if (max_auc <= 0.0) {
    return R_NaN;
  }

  std::vector<int> x;
  x.reserve(ranks.size());
  for (std::vector<int>::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
    if (*it > 0 && *it < auc_threshold) {
      x.push_back(*it);
    }
  }
  if (x.empty()) {
    return 0.0;
  }
  std::sort(x.begin(), x.end());

  double auc = 0.0;
  int prev = x[0];
  for (std::size_t i = 1; i < x.size(); ++i) {
    auc += static_cast<double>(x[i] - prev) * static_cast<double>(i);
    prev = x[i];
  }
  auc += static_cast<double>(auc_threshold - prev) * static_cast<double>(x.size());
  return auc / max_auc;
}

static double ctxcore_auc_from_ranks(
  const std::vector<int>& ranks,
  int auc_threshold,
  int n_genes
) {
  const int rank_cutoff = std::max(0, auc_threshold - 1);
  if (rank_cutoff <= 0 || ranks.empty()) {
    return 0.0;
  }

  std::vector<int> x;
  x.reserve(ranks.size());
  for (std::vector<int>::const_iterator it = ranks.begin(); it != ranks.end(); ++it) {
    if (*it > 0 && (*it - 1) < rank_cutoff) {
      x.push_back(*it - 1);
    }
  }
  if (x.empty()) {
    return 0.0;
  }
  std::sort(x.begin(), x.end());

  double auc = 0.0;
  int prev = x[0];
  double cumulative = 1.0;
  for (std::size_t i = 1; i < x.size(); ++i) {
    auc += static_cast<double>(x[i] - prev) * cumulative;
    prev = x[i];
    cumulative += 1.0;
  }
  auc += static_cast<double>(rank_cutoff - prev) * cumulative;
  const double max_auc = static_cast<double>(rank_cutoff + 1) * static_cast<double>(n_genes);
  return max_auc > 0.0 ? auc / max_auc : R_NaN;
}

// [[Rcpp::export]]
NumericMatrix aucell_auc_sparse(
  S4 expr,
  List gene_sets,
  int auc_max_rank,
  bool norm_auc = true,
  int strategy = 1,
  int algorithm = 1
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();
  const int auc_threshold = std::max(1, auc_max_rank);

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<std::vector<int> > sets(n_sets);
  std::vector<double> max_auc(n_sets);
  std::vector<int> set_gene_union;
  set_gene_union.reserve(n_sets * 64);
  std::vector<int> zero_order;

  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
        set_gene_union.push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
    max_auc[set_i] = aucell_max_auc(static_cast<int>(sets[set_i].size()), auc_threshold, norm_auc);
  }
  std::sort(set_gene_union.begin(), set_gene_union.end());
  set_gene_union.erase(std::unique(set_gene_union.begin(), set_gene_union.end()), set_gene_union.end());

  NumericMatrix scores(n_cells, n_sets);
  std::vector<int> rank_by_gene(n_genes, 0);
  std::vector<double> value_by_gene(n_genes, 0.0);
  std::vector<int> touched_values;
  std::vector<int> touched_ranks;
  std::vector<AucEntry> entries;
  entries.reserve(n_genes);
  const int top_n = std::max(0, std::min(n_genes, auc_threshold - 1));
  if (strategy == 1 && top_n > 0) {
    zero_order.resize(n_genes);
    std::iota(zero_order.begin(), zero_order.end(), 0);
    std::sort(zero_order.begin(), zero_order.end(), [](int a, int b) {
      return gene_hash_tiebreak(a) < gene_hash_tiebreak(b);
    });
  }

  for (int cell = 0; cell < n_cells; ++cell) {
    touched_values.clear();
    touched_ranks.clear();
    entries.clear();

    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (!R_finite(value) || value == 0.0) {
        continue;
      }
      if (strategy == 2) {
        entries.push_back(AucEntry{gene, value, gene_hash_tiebreak(gene)});
      } else {
        value_by_gene[gene] = value;
        touched_values.push_back(gene);
      }
    }

    if (strategy == 2) {
      std::sort(entries.begin(), entries.end(), aucell_entry_before);

      for (std::size_t rank_i = 0; rank_i < entries.size(); ++rank_i) {
        const int gene = entries[rank_i].gene;
        rank_by_gene[gene] = static_cast<int>(rank_i) + 1;
        touched_ranks.push_back(gene);
      }

      if (static_cast<int>(entries.size()) < auc_threshold) {
        const int n_zero = n_genes - static_cast<int>(entries.size());
        for (std::vector<int>::const_iterator it = set_gene_union.begin(); it != set_gene_union.end(); ++it) {
          const int gene = *it;
          if (rank_by_gene[gene] == 0 && n_zero > 0) {
            // Deterministic zero-gene rank: hash-based offset into [K+1, n_genes]
            const int offset = static_cast<int>(gene_hash_tiebreak(gene + 7777777) * static_cast<double>(n_zero));
            rank_by_gene[gene] = static_cast<int>(entries.size()) + 1 + offset;
            touched_ranks.push_back(gene);
          }
        }
      }
    } else {
      if (strategy == 1) {
        if (top_n > 0) {
          entries.reserve(touched_values.size() + top_n);
          for (std::vector<int>::const_iterator it = touched_values.begin(); it != touched_values.end(); ++it) {
            const int gene = *it;
            entries.push_back(AucEntry{gene, value_by_gene[gene], gene_hash_tiebreak(gene)});
          }
          if (top_n < static_cast<int>(entries.size())) {
            std::nth_element(entries.begin(), entries.begin() + top_n, entries.end(), aucell_entry_before);
            entries.resize(top_n);
          } else if (static_cast<int>(entries.size()) < top_n) {
            for (std::vector<int>::const_iterator it = zero_order.begin(); it != zero_order.end() && static_cast<int>(entries.size()) < top_n; ++it) {
              const int gene = *it;
              if (value_by_gene[gene] == 0.0) {
                entries.push_back(AucEntry{gene, 0.0, gene_hash_tiebreak(gene)});
              }
            }
          }
          std::sort(entries.begin(), entries.end(), aucell_entry_before);
        }
      } else {
        for (int gene = 0; gene < n_genes; ++gene) {
          entries.push_back(AucEntry{gene, value_by_gene[gene], gene_hash_tiebreak(gene)});
        }
        std::sort(entries.begin(), entries.end(), aucell_entry_before);
      }

      for (std::size_t rank_i = 0; rank_i < entries.size(); ++rank_i) {
        const int gene = entries[rank_i].gene;
        rank_by_gene[gene] = static_cast<int>(rank_i) + 1;
        touched_ranks.push_back(gene);
      }
    }

    for (int set_i = 0; set_i < n_sets; ++set_i) {
      std::vector<int> ranks;
      ranks.reserve(sets[set_i].size());
      for (std::vector<int>::const_iterator it = sets[set_i].begin(); it != sets[set_i].end(); ++it) {
        const int rank = rank_by_gene[*it];
        if (rank > 0) {
          ranks.push_back(rank);
        }
      }
      if (algorithm == 2) {
        scores(cell, set_i) = ctxcore_auc_from_ranks(
          ranks,
          auc_threshold,
          static_cast<int>(sets[set_i].size())
        );
      } else {
        scores(cell, set_i) = aucell_auc_from_ranks(ranks, auc_threshold, max_auc[set_i]);
      }
    }

    for (std::vector<int>::const_iterator it = touched_ranks.begin(); it != touched_ranks.end(); ++it) {
      rank_by_gene[*it] = 0;
    }
    for (std::vector<int>::const_iterator it = touched_values.begin(); it != touched_values.end(); ++it) {
      value_by_gene[*it] = 0.0;
    }
  }

  return scores;
}

// [[Rcpp::export]]
NumericMatrix aucell_auc_ranked(
  NumericMatrix rankings,
  List gene_sets,
  int auc_max_rank
) {
  const int n_cells = rankings.nrow();
  const int n_genes = rankings.ncol();
  const int n_sets = gene_sets.size();
  const int auc_threshold = std::max(1, auc_max_rank);

  std::vector<std::vector<int> > sets(n_sets);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
  }

  NumericMatrix scores(n_cells, n_sets);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int set_i = 0; set_i < n_sets; ++set_i) {
      std::vector<int> ranks;
      ranks.reserve(sets[set_i].size());
      for (std::vector<int>::const_iterator it = sets[set_i].begin(); it != sets[set_i].end(); ++it) {
        const double rank0 = rankings(cell, *it);
        if (R_finite(rank0)) {
          ranks.push_back(static_cast<int>(rank0) + 1);
        }
      }
      scores(cell, set_i) = ctxcore_auc_from_ranks(
        ranks,
        auc_threshold,
        static_cast<int>(sets[set_i].size())
      );
    }
  }
  return scores;
}

// [[Rcpp::export]]
DataFrame ora_hypergeom(
  CharacterVector genes,
  CharacterVector term_ids,
  CharacterVector term_genes,
  CharacterVector term_name_ids,
  CharacterVector term_names,
  int min_size = 10,
  int max_size = 2147483647
) {
  if (term_ids.size() != term_genes.size()) {
    stop("term_ids and term_genes must have the same length");
  }

  std::unordered_map<std::string, std::unordered_set<std::string> > term_to_genes;
  std::unordered_map<std::string, std::string> term_to_name;
  std::vector<std::string> term_order;
  std::unordered_set<std::string> term_seen;
  std::unordered_set<std::string> universe;

  for (int i = 0; i < term_name_ids.size() && i < term_names.size(); ++i) {
    if (CharacterVector::is_na(term_name_ids[i]) || CharacterVector::is_na(term_names[i])) {
      continue;
    }
    term_to_name[as<std::string>(term_name_ids[i])] = as<std::string>(term_names[i]);
  }

  for (int i = 0; i < term_ids.size(); ++i) {
    if (CharacterVector::is_na(term_ids[i]) || CharacterVector::is_na(term_genes[i])) {
      continue;
    }
    const std::string term = as<std::string>(term_ids[i]);
    const std::string gene = as<std::string>(term_genes[i]);
    if (term.empty() || gene.empty()) {
      continue;
    }
    if (term_seen.insert(term).second) {
      term_order.push_back(term);
    }
    term_to_genes[term].insert(gene);
    universe.insert(gene);
  }

  std::vector<std::string> query_order;
  std::unordered_set<std::string> query_seen;
  std::unordered_map<std::string, int> query_rank;
  query_order.reserve(genes.size());
  for (int i = 0; i < genes.size(); ++i) {
    if (CharacterVector::is_na(genes[i])) {
      continue;
    }
    const std::string gene = as<std::string>(genes[i]);
    if (gene.empty() || universe.find(gene) == universe.end()) {
      continue;
    }
    if (query_seen.insert(gene).second) {
      query_rank[gene] = static_cast<int>(query_order.size());
      query_order.push_back(gene);
    }
  }

  const int universe_size = static_cast<int>(universe.size());
  const int query_size = static_cast<int>(query_order.size());
  std::vector<std::string> out_id;
  std::vector<std::string> out_description;
  std::vector<std::string> out_gene_id;
  std::vector<int> out_count;
  std::vector<int> out_query_size;
  std::vector<int> out_term_size;
  std::vector<int> out_universe_size;
  std::vector<double> out_pvalue;

  if (universe_size <= 0 || query_size <= 0) {
    return DataFrame::create(
      _["ID"] = CharacterVector(0),
      _["Description"] = CharacterVector(0),
      _["geneID"] = CharacterVector(0),
      _["Count"] = IntegerVector(0),
      _["query_size"] = IntegerVector(0),
      _["term_size"] = IntegerVector(0),
      _["universe_size"] = IntegerVector(0),
      _["pvalue"] = NumericVector(0),
      _["stringsAsFactors"] = false
    );
  }

  for (std::vector<std::string>::const_iterator term_it = term_order.begin(); term_it != term_order.end(); ++term_it) {
    const std::string& term = *term_it;
    const std::unordered_set<std::string>& term_set = term_to_genes[term];
    const int term_size = static_cast<int>(term_set.size());
    if (term_size < min_size || term_size > max_size || term_size <= 0 || term_size >= universe_size) {
      continue;
    }

    std::vector<std::pair<int, std::string> > hits;
    hits.reserve(std::min(term_size, query_size));
    for (std::unordered_set<std::string>::const_iterator gene_it = term_set.begin(); gene_it != term_set.end(); ++gene_it) {
      std::unordered_map<std::string, int>::const_iterator rank_it = query_rank.find(*gene_it);
      if (rank_it != query_rank.end()) {
        hits.push_back(std::make_pair(rank_it->second, *gene_it));
      }
    }
    const int count = static_cast<int>(hits.size());
    if (count <= 0) {
      continue;
    }
    std::sort(hits.begin(), hits.end());

    std::string hit_text;
    for (std::size_t hit_i = 0; hit_i < hits.size(); ++hit_i) {
      if (hit_i > 0) {
        hit_text += "/";
      }
      hit_text += hits[hit_i].second;
    }

    const double pvalue = R::phyper(
      static_cast<double>(count - 1),
      static_cast<double>(term_size),
      static_cast<double>(universe_size - term_size),
      static_cast<double>(query_size),
      false,
      false
    );

    out_id.push_back(term);
    std::unordered_map<std::string, std::string>::const_iterator name_it = term_to_name.find(term);
    out_description.push_back(name_it == term_to_name.end() ? term : name_it->second);
    out_gene_id.push_back(hit_text);
    out_count.push_back(count);
    out_query_size.push_back(query_size);
    out_term_size.push_back(term_size);
    out_universe_size.push_back(universe_size);
    out_pvalue.push_back(pvalue);
  }

  return DataFrame::create(
    _["ID"] = out_id,
    _["Description"] = out_description,
    _["geneID"] = out_gene_id,
    _["Count"] = out_count,
    _["query_size"] = out_query_size,
    _["term_size"] = out_term_size,
    _["universe_size"] = out_universe_size,
    _["pvalue"] = out_pvalue,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
NumericMatrix module_score_sparse(
  S4 expr,
  List feature_sets,
  List control_sets
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = feature_sets.size();
  if (control_sets.size() != n_sets) {
    stop("feature_sets and control_sets must have the same length");
  }

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  NumericMatrix scores(n_cells, n_sets);
  std::vector<std::vector<std::pair<int, double> > > gene_to_scores(n_genes);
  std::vector<bool> valid_sets(n_sets, true);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector feature_idx = feature_sets[set_i];
    IntegerVector control_idx = control_sets[set_i];
    std::vector<int> features;
    std::vector<int> controls;
    features.reserve(feature_idx.size());
    controls.reserve(control_idx.size());
    for (int i = 0; i < feature_idx.size(); ++i) {
      const int gene = feature_idx[i] - 1;
      if (gene >= 0 && gene < n_genes) {
        features.push_back(gene);
      }
    }
    for (int i = 0; i < control_idx.size(); ++i) {
      const int gene = control_idx[i] - 1;
      if (gene >= 0 && gene < n_genes) {
        controls.push_back(gene);
      }
    }
    if (features.empty() || controls.empty()) {
      valid_sets[set_i] = false;
      continue;
    }
    const double feature_weight = 1.0 / static_cast<double>(features.size());
    const double control_weight = -1.0 / static_cast<double>(controls.size());
    for (std::vector<int>::const_iterator it = features.begin(); it != features.end(); ++it) {
      gene_to_scores[*it].push_back(std::make_pair(set_i, feature_weight));
    }
    for (std::vector<int>::const_iterator it = controls.begin(); it != controls.end(); ++it) {
      gene_to_scores[*it].push_back(std::make_pair(set_i, control_weight));
    }
  }

  for (int set_i = 0; set_i < n_sets; ++set_i) {
    if (!valid_sets[set_i]) {
      for (int cell = 0; cell < n_cells; ++cell) {
        scores(cell, set_i) = R_NaN;
      }
    }
  }

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (!R_finite(value) || value == 0.0) {
        continue;
      }
      const std::vector<std::pair<int, double> >& memberships = gene_to_scores[gene];
      for (std::vector<std::pair<int, double> >::const_iterator it = memberships.begin(); it != memberships.end(); ++it) {
        scores(cell, it->first) += value * it->second;
      }
    }
  }

  return scores;
}

static double quantile_type7(std::vector<double>& x, double prob) {
  if (x.empty()) {
    return R_NaN;
  }
  std::sort(x.begin(), x.end());
  if (x.size() == 1) {
    return x[0];
  }
  const double h = 1.0 + (static_cast<double>(x.size()) - 1.0) * prob;
  const int h_floor = static_cast<int>(std::floor(h));
  const double frac = h - static_cast<double>(h_floor);
  const double lower = x[static_cast<std::size_t>(h_floor - 1)];
  if (h_floor >= static_cast<int>(x.size())) {
    return lower;
  }
  const double upper = x[static_cast<std::size_t>(h_floor)];
  return lower + frac * (upper - lower);
}

static double log2_fraction_diff(
  int count_1,
  int total_1,
  int count_2,
  int total_2,
  double pseudocount
) {
  const double frac_1 = static_cast<double>(count_1) / static_cast<double>(total_1);
  const double frac_2 = static_cast<double>(count_2) / static_cast<double>(total_2);
  return std::log((frac_2 + pseudocount) / (frac_1 + pseudocount)) / std::log(2.0);
}

// [[Rcpp::export]]
DataFrame proportion_permutation(
  IntegerVector sample_ids,
  IntegerVector cluster_ids,
  CharacterVector cluster_levels,
  int n_permutations,
  double pseudocount = 1e-8,
  bool verbose = false
) {
  if (sample_ids.size() != cluster_ids.size()) {
    stop("sample_ids and cluster_ids must have the same length");
  }
  const int n_cells = sample_ids.size();
  const int n_clusters = cluster_levels.size();
  if (n_cells == 0 || n_clusters == 0) {
    return DataFrame::create(
      _["clusters"] = CharacterVector(0),
      _["fraction_1"] = NumericVector(0),
      _["fraction_2"] = NumericVector(0),
      _["obs_log2FD"] = NumericVector(0),
      _["pval"] = NumericVector(0),
      _["boot_mean_log2FD"] = NumericVector(0),
      _["boot_CI_2.5"] = NumericVector(0),
      _["boot_CI_97.5"] = NumericVector(0),
      _["stringsAsFactors"] = false
    );
  }

  std::vector<int> labels(n_cells);
  std::vector<int> clusters(n_cells);
  std::vector<int> group1_clusters;
  std::vector<int> group2_clusters;
  int total_1 = 0;
  int total_2 = 0;
  for (int i = 0; i < n_cells; ++i) {
    const int label = sample_ids[i];
    const int cluster = cluster_ids[i] - 1;
    if ((label != 1 && label != 2) || cluster < 0 || cluster >= n_clusters) {
      stop("sample_ids must be 1/2 and cluster_ids must be within cluster_levels");
    }
    labels[i] = label;
    clusters[i] = cluster;
    if (label == 1) {
      ++total_1;
      group1_clusters.push_back(cluster);
    } else {
      ++total_2;
      group2_clusters.push_back(cluster);
    }
  }
  if (total_1 == 0 || total_2 == 0) {
    stop("Both comparison groups must contain cells");
  }

  std::vector<int> obs_1(n_clusters, 0);
  std::vector<int> obs_2(n_clusters, 0);
  for (int i = 0; i < n_cells; ++i) {
    if (labels[i] == 1) {
      ++obs_1[clusters[i]];
    } else {
      ++obs_2[clusters[i]];
    }
  }

  std::vector<double> obs_log2(n_clusters);
  std::vector<int> increased(n_clusters, 0);
  std::vector<int> decreased(n_clusters, 0);
  for (int k = 0; k < n_clusters; ++k) {
    obs_log2[k] = log2_fraction_diff(obs_1[k], total_1, obs_2[k], total_2, pseudocount);
  }

  std::vector<int> perm_labels(labels);
  std::vector<int> perm_1(n_clusters, 0);
  std::vector<int> perm_2(n_clusters, 0);
  scop::CliProgress progress(
    2 * n_permutations,
    verbose,
    "Run proportion permutation and bootstrap"
  );
  for (int perm = 0; perm < n_permutations; ++perm) {
    if (scop::should_check_interrupt(perm, 2 * n_permutations, verbose)) {
      Rcpp::checkUserInterrupt();
    }
    progress.set(perm);
    perm_labels = labels;
    for (int i = n_cells - 1; i > 0; --i) {
      const int j = static_cast<int>(std::floor(unif_rand() * static_cast<double>(i + 1)));
      std::swap(perm_labels[i], perm_labels[j]);
    }
    std::fill(perm_1.begin(), perm_1.end(), 0);
    std::fill(perm_2.begin(), perm_2.end(), 0);
    for (int i = 0; i < n_cells; ++i) {
      if (perm_labels[i] == 1) {
        ++perm_1[clusters[i]];
      } else {
        ++perm_2[clusters[i]];
      }
    }
    for (int k = 0; k < n_clusters; ++k) {
      const double perm_log2 = log2_fraction_diff(perm_1[k], total_1, perm_2[k], total_2, pseudocount);
      if (obs_log2[k] <= perm_log2) {
        ++increased[k];
      }
      if (obs_log2[k] >= perm_log2) {
        ++decreased[k];
      }
    }
  }
  progress.set(n_permutations);

  std::vector<double> boot_values(static_cast<std::size_t>(n_clusters) * std::max(1, n_permutations), R_NaN);
  std::vector<int> boot_1(n_clusters, 0);
  std::vector<int> boot_2(n_clusters, 0);
  for (int perm = 0; perm < n_permutations; ++perm) {
    const int progress_value = n_permutations + perm;
    if (scop::should_check_interrupt(progress_value, 2 * n_permutations, verbose)) {
      Rcpp::checkUserInterrupt();
    }
    progress.set(progress_value);
    std::fill(boot_1.begin(), boot_1.end(), 0);
    std::fill(boot_2.begin(), boot_2.end(), 0);
    for (int i = 0; i < total_1; ++i) {
      const int draw = static_cast<int>(std::floor(unif_rand() * static_cast<double>(total_1)));
      ++boot_1[group1_clusters[draw]];
    }
    for (int i = 0; i < total_2; ++i) {
      const int draw = static_cast<int>(std::floor(unif_rand() * static_cast<double>(total_2)));
      ++boot_2[group2_clusters[draw]];
    }
    for (int k = 0; k < n_clusters; ++k) {
      boot_values[static_cast<std::size_t>(k) * n_permutations + perm] =
        log2_fraction_diff(boot_1[k], total_1, boot_2[k], total_2, pseudocount);
    }
  }
  progress.set(2 * n_permutations, true);

  NumericVector fraction_1(n_clusters);
  NumericVector fraction_2(n_clusters);
  NumericVector pval(n_clusters);
  NumericVector boot_mean(n_clusters);
  NumericVector boot_low(n_clusters);
  NumericVector boot_high(n_clusters);

  for (int k = 0; k < n_clusters; ++k) {
    fraction_1[k] = static_cast<double>(obs_1[k]) / static_cast<double>(total_1);
    fraction_2[k] = static_cast<double>(obs_2[k]) / static_cast<double>(total_2);
    pval[k] = obs_log2[k] > 0.0 ?
      (static_cast<double>(increased[k] + 1) / static_cast<double>(n_permutations + 1)) :
      (static_cast<double>(decreased[k] + 1) / static_cast<double>(n_permutations + 1));

    std::vector<double> one_cluster;
    one_cluster.reserve(n_permutations);
    double sum = 0.0;
    for (int perm = 0; perm < n_permutations; ++perm) {
      const double value = boot_values[static_cast<std::size_t>(k) * n_permutations + perm];
      if (R_finite(value)) {
        one_cluster.push_back(value);
        sum += value;
      }
    }
    if (one_cluster.empty()) {
      boot_mean[k] = R_NaN;
      boot_low[k] = R_NaN;
      boot_high[k] = R_NaN;
    } else {
      boot_mean[k] = sum / static_cast<double>(one_cluster.size());
      std::vector<double> q_values = one_cluster;
      boot_low[k] = quantile_type7(q_values, 0.025);
      q_values = one_cluster;
      boot_high[k] = quantile_type7(q_values, 0.975);
    }
  }

  return DataFrame::create(
    _["clusters"] = cluster_levels,
    _["fraction_1"] = fraction_1,
    _["fraction_2"] = fraction_2,
    _["obs_log2FD"] = obs_log2,
    _["pval"] = pval,
    _["boot_mean_log2FD"] = boot_mean,
    _["boot_CI_2.5"] = boot_low,
    _["boot_CI_97.5"] = boot_high,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
NumericMatrix ssgsea_rank_dense(
  S4 expr,
  List gene_sets,
  double alpha = 0.25,
  bool normalize = true
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<std::vector<int> > sets(n_sets);
  std::vector<int> set_sizes(n_sets, 0);
  std::vector<unsigned char> valid_sets(n_sets, 0);
  std::vector<int> set_gene_union;
  set_gene_union.reserve(n_sets * 64);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
        set_gene_union.push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
    set_sizes[set_i] = static_cast<int>(sets[set_i].size());
    valid_sets[set_i] = (set_sizes[set_i] > 0 && set_sizes[set_i] < n_genes) ? 1 : 0;
  }
  std::sort(set_gene_union.begin(), set_gene_union.end());
  set_gene_union.erase(std::unique(set_gene_union.begin(), set_gene_union.end()), set_gene_union.end());

  NumericMatrix scores(n_cells, n_sets);
  std::vector<int> rank_by_gene(n_genes);
  std::vector<int> position_by_gene(n_genes);
  std::vector<double> rank_weight_by_gene(n_genes);
  double min_score = R_PosInf;
  double max_score = R_NegInf;

  // Sparse ranking optimization: only sort non-zero genes per cell.
  // Zero-expression genes (the vast majority in scRNA-seq) are all tied
  // at the same value (0.0) and get identical rank/weight, so we handle
  // them as a single block instead of sorting them individually.
  std::vector<int> nonzero_genes;
  nonzero_genes.reserve(n_genes);
  std::vector<double> nonzero_values;
  nonzero_values.reserve(n_genes);
  std::vector<unsigned char> is_nonzero(n_genes, 0);
  std::vector<int> nz_order;
  nz_order.reserve(n_genes);
  const double total_position_sum = static_cast<double>(n_genes) *
    static_cast<double>(n_genes + 1) / 2.0;

  for (int cell = 0; cell < n_cells; ++cell) {
    // Collect non-zero expression values for this cell
    nonzero_genes.clear();
    nonzero_values.clear();

    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (!R_finite(value) || value == 0.0) {
        continue;
      }
      nonzero_genes.push_back(gene);
      nonzero_values.push_back(value);
      is_nonzero[gene] = 1;
    }

    const int n_nonzero = static_cast<int>(nonzero_genes.size());
    const int n_zero = n_genes - n_nonzero;

    int n_negative = 0;
    if (n_nonzero > 0) {
      for (int i = 0; i < n_nonzero; ++i) {
        if (nonzero_values[i] < 0.0) {
          ++n_negative;
        }
      }
    }

    // Zero genes are all tied at value 0.0. Their rank starts after negative
    // non-zero values and uses the original integer-truncated average rank.
    const int zero_rank_int = (n_zero > 0) ?
      n_negative + static_cast<int>((1.0 + static_cast<double>(n_zero)) / 2.0) : 0;
    const double zero_rank_weight = std::pow(std::fabs(static_cast<double>(zero_rank_int)), alpha);

    // Step 1: initialize only genes that can be queried by gene sets. The
    // original dense view assigns the same zero rank to all structural zeros,
    // but downstream scoring only reads genes present in at least one set.
    for (std::vector<int>::const_iterator it = set_gene_union.begin(); it != set_gene_union.end(); ++it) {
      rank_by_gene[*it] = zero_rank_int;
      rank_weight_by_gene[*it] = zero_rank_weight;
    }

    // Step 2: compute non-zero gene ranks (overwrites the zero defaults)
    if (n_nonzero > 0) {
      nz_order.resize(n_nonzero);
      std::iota(nz_order.begin(), nz_order.end(), 0);
      // Sort non-zero genes by value ascending (gene-index tiebreak)
      std::sort(
        nz_order.begin(),
        nz_order.end(),
        [&](int a, int b) {
          if (nonzero_values[a] < nonzero_values[b]) return true;
          if (nonzero_values[a] > nonzero_values[b]) return false;
          return nonzero_genes[a] < nonzero_genes[b];
        }
      );

      // Compute tied ranks. Positive values start after the zero block, while
      // negative values keep their position before zeros in the full ordering.
      int tie_start = 0;
      while (tie_start < n_nonzero) {
        int tie_end = tie_start + 1;
        const double tie_value = nonzero_values[nz_order[tie_start]];
        while (
          tie_end < n_nonzero &&
          nonzero_values[nz_order[tie_end]] == tie_value
        ) {
          ++tie_end;
        }
        const double zero_offset = (tie_value > 0.0) ? static_cast<double>(n_zero) : 0.0;
        const double rank_val = zero_offset +
          (static_cast<double>(tie_start + 1) + static_cast<double>(tie_end)) / 2.0;
        const int rank_int = static_cast<int>(rank_val);
        const double rank_weight = std::pow(std::fabs(static_cast<double>(rank_int)), alpha);
        for (int pos = tie_start; pos < tie_end; ++pos) {
          const int gene = nonzero_genes[nz_order[pos]];
          rank_by_gene[gene] = rank_int;
          rank_weight_by_gene[gene] = rank_weight;
        }
        tie_start = tie_end;
      }

      // Step 3: build positions by descending rank. Insert the zero block
      // between positive and negative non-zero ranks to match a full sort.
      std::sort(
        nz_order.begin(),
        nz_order.end(),
        [&](int a, int b) {
          const int ra = rank_by_gene[nonzero_genes[a]];
          const int rb = rank_by_gene[nonzero_genes[b]];
          if (ra > rb) return true;
          if (ra < rb) return false;
          return nonzero_genes[a] < nonzero_genes[b];
        }
      );
      int pos = 0;
      int nz_pos = 0;
      while (
        nz_pos < n_nonzero &&
        (n_zero == 0 || rank_by_gene[nonzero_genes[nz_order[nz_pos]]] > zero_rank_int)
      ) {
        position_by_gene[nonzero_genes[nz_order[nz_pos]]] = ++pos;
        ++nz_pos;
      }

      if (n_zero > 0) {
        const int zero_start_pos = pos;
        int nonzero_leq = 0;
        for (std::vector<int>::const_iterator it = set_gene_union.begin(); it != set_gene_union.end(); ++it) {
          const int gene = *it;
          while (nonzero_leq < n_nonzero && nonzero_genes[nonzero_leq] <= gene) {
            ++nonzero_leq;
          }
          if (is_nonzero[gene] == 0) {
            position_by_gene[gene] = zero_start_pos + gene + 1 - nonzero_leq;
          }
        }
        pos += n_zero;
      }

      while (nz_pos < n_nonzero) {
        position_by_gene[nonzero_genes[nz_order[nz_pos]]] = ++pos;
        ++nz_pos;
      }
    } else {
      for (std::vector<int>::const_iterator it = set_gene_union.begin(); it != set_gene_union.end(); ++it) {
        position_by_gene[*it] = *it + 1;
      }
    }

    // Step 5: compute ES scores for each gene set
    for (int set_i = 0; set_i < n_sets; ++set_i) {
      const std::vector<int>& set = sets[set_i];
      const int set_size = set_sizes[set_i];
      if (!valid_sets[set_i]) {
        scores(cell, set_i) = R_NaN;
        continue;
      }

      double in_weight_sum = 0.0;
      double in_weighted_position_sum = 0.0;
      double in_position_sum = 0.0;
      for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
        const int gene = *it;
        const double inverse_position = static_cast<double>(n_genes - position_by_gene[gene] + 1);
        const double gene_weight = rank_weight_by_gene[gene];
        in_weight_sum += gene_weight;
        in_weighted_position_sum += gene_weight * inverse_position;
        in_position_sum += inverse_position;
      }

      if (in_weight_sum <= 0.0 || n_genes == set_size) {
        scores(cell, set_i) = R_NaN;
        continue;
      }

      const double in_step = in_weighted_position_sum / in_weight_sum;
      const double out_step = (total_position_sum - in_position_sum) /
        static_cast<double>(n_genes - set_size);
      const double score = in_step - out_step;
      scores(cell, set_i) = score;
      if (R_finite(score)) {
        if (score < min_score) {
          min_score = score;
        }
        if (score > max_score) {
          max_score = score;
        }
      }
    }

    // Cleanup: reset is_nonzero markers for next cell
    for (int i = 0; i < n_nonzero; ++i) {
      is_nonzero[nonzero_genes[i]] = 0;
    }
  }

  if (normalize) {
    const double range = max_score - min_score;
    for (int set_i = 0; set_i < n_sets; ++set_i) {
      for (int cell = 0; cell < n_cells; ++cell) {
        scores(cell, set_i) = scores(cell, set_i) / range;
      }
    }
  }

  return scores;
}

// [[Rcpp::export]]
NumericMatrix zscore_dense(
  S4 expr,
  List gene_sets,
  int min_size = 1,
  int max_size = 2147483647
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  // Build gene -> set adjacency first so sparse passes can ignore genes that
  // are never queried by any retained gene set.
  std::vector<std::vector<int> > gene_to_sets(n_genes);
  std::vector<unsigned char> gene_needed(n_genes, 0);
  std::vector<int> set_sizes(n_sets, 0);
  std::vector<double> inv_sqrt_size(n_sets, 1.0);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    std::vector<int> seen;
    seen.reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        seen.push_back(gene);
      }
    }
    std::sort(seen.begin(), seen.end());
    seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
    for (std::vector<int>::const_iterator it = seen.begin(); it != seen.end(); ++it) {
      const int gene = *it;
      gene_to_sets[gene].push_back(set_i);
      gene_needed[gene] = 1;
      ++set_sizes[set_i];
    }
    if (set_sizes[set_i] > 0) {
      inv_sqrt_size[set_i] = 1.0 / std::sqrt(static_cast<double>(set_sizes[set_i]));
    }
  }

  std::vector<double> gene_sums(n_genes, 0.0);
  std::vector<double> gene_sq_sums(n_genes, 0.0);

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (gene_needed[gene] && R_finite(value)) {
        gene_sums[gene] += value;
        gene_sq_sums[gene] += value * value;
      }
    }
  }

  std::vector<double> gene_means(n_genes, 0.0);
  std::vector<double> gene_sds(n_genes, 1.0);
  if (n_cells > 1) {
    for (int gene = 0; gene < n_genes; ++gene) {
      if (!gene_needed[gene]) {
        gene_means[gene] = R_NaN;
        gene_sds[gene] = R_NaN;
        continue;
      }
      gene_means[gene] = gene_sums[gene] / static_cast<double>(n_cells);
      double var_val = (gene_sq_sums[gene] -
        static_cast<double>(n_cells) * gene_means[gene] * gene_means[gene]) /
        static_cast<double>(n_cells - 1);
      if (var_val < 0.0 && var_val > -1e-12) {
        var_val = 0.0;
      }
      if (R_finite(var_val) && var_val > 0.0) {
        gene_sds[gene] = std::sqrt(var_val);
      } else {
        gene_means[gene] = R_NaN;
        gene_sds[gene] = R_NaN;
      }
    }
  }
  std::vector<double>().swap(gene_sums);
  std::vector<double>().swap(gene_sq_sums);

  // Compute z-score per gene set over the dense expression view used by
  // GSVA::zscoreParam in the R benchmark. Structural zeros contribute
  // (0 - mean) / sd, then stored values add value / sd.
  NumericMatrix scores(n_cells, n_sets);
  std::vector<double> zero_sum(n_sets, 0.0);
  std::vector<int> valid_set_sizes(n_sets, 0);
  for (int gene = 0; gene < n_genes; ++gene) {
    if (!R_finite(gene_means[gene]) || !R_finite(gene_sds[gene]) || gene_sds[gene] <= 0.0) {
      continue;
    }
    const std::vector<int>& sets_for_gene = gene_to_sets[gene];
    if (sets_for_gene.empty()) continue;
    const double zero_contrib = -gene_means[gene] / gene_sds[gene];
    for (std::vector<int>::const_iterator it = sets_for_gene.begin(); it != sets_for_gene.end(); ++it) {
      zero_sum[*it] += zero_contrib;
      ++valid_set_sizes[*it];
    }
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int set_i = 0; set_i < n_sets; ++set_i) {
      scores(cell, set_i) = (valid_set_sizes[set_i] >= min_size && valid_set_sizes[set_i] <= max_size)
        ? zero_sum[set_i] : R_NaN;
    }
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (!R_finite(value) || !R_finite(gene_sds[gene]) || gene_sds[gene] <= 0.0) continue;
      const std::vector<int>& sets_for_gene = gene_to_sets[gene];
      if (sets_for_gene.empty()) continue;
      const double contrib = value / gene_sds[gene];
      for (std::vector<int>::const_iterator it = sets_for_gene.begin(); it != sets_for_gene.end(); ++it) {
        scores(cell, *it) += contrib;
      }
    }
  }
  // Final scaling by 1/√n
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int set_i = 0; set_i < n_sets; ++set_i) {
      if (R_finite(scores(cell, set_i))) {
        scores(cell, set_i) *= (valid_set_sizes[set_i] > 0)
          ? 1.0 / std::sqrt(static_cast<double>(valid_set_sizes[set_i]))
          : inv_sqrt_size[set_i];
      }
    }
  }

  return scores;
}

// [[Rcpp::export]]
NumericMatrix plage_dense(
  S4 expr,
  List gene_sets,
  int min_size = 1,
  int max_size = 2147483647
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<std::vector<int> > sets(n_sets);
  std::vector<unsigned char> row_needed(n_genes, 0);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
    for (std::vector<int>::const_iterator it = sets[set_i].begin(); it != sets[set_i].end(); ++it) {
      row_needed[*it] = 1;
    }
  }

  std::vector<double> row_sums(n_genes, 0.0);
  std::vector<double> row_sq_sums(n_genes, 0.0);

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (row_needed[gene] && R_finite(value) && value != 0.0) {
        row_sums[gene] += value;
        row_sq_sums[gene] += value * value;
      }
    }
  }

  std::vector<int> row_ptr(n_genes + 1, 0);
  std::vector<int> row_counts(n_genes, 0);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (row_needed[gene] && R_finite(value) && value != 0.0) {
        ++row_counts[gene];
      }
    }
  }
  for (int gene = 0; gene < n_genes; ++gene) {
    row_ptr[gene + 1] = row_ptr[gene] + row_counts[gene];
  }
  const int nnz = row_ptr[n_genes];
  std::vector<int> row_cursor(row_ptr);
  std::vector<int> row_cells(nnz);
  std::vector<double> row_values(nnz);

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (row_needed[gene] && R_finite(value) && value != 0.0) {
        const int out = row_cursor[gene]++;
        row_cells[out] = cell;
        row_values[out] = value;
      }
    }
  }
  std::vector<int>().swap(row_cursor);

  std::vector<double> row_means(n_genes, R_NaN);
  std::vector<double> row_sds(n_genes, R_NaN);
  std::vector<unsigned char> row_valid(n_genes, 0);
  std::vector<double> orient_row_means(n_genes, R_NaN);
  std::vector<double> orient_row_sds(n_genes, R_NaN);
  std::vector<unsigned char> orient_row_valid(n_genes, 0);
  if (n_cells > 1) {
    for (int gene = 0; gene < n_genes; ++gene) {
      if (!row_needed[gene]) {
        continue;
      }
      const double mean = row_sums[gene] / static_cast<double>(n_cells);
      double var = (row_sq_sums[gene] - static_cast<double>(n_cells) * mean * mean) /
        static_cast<double>(n_cells - 1);
      if (var < 0.0 && var > -1e-12) {
        var = 0.0;
      }
      if (R_finite(var) && var > 0.0) {
        row_means[gene] = mean;
        row_sds[gene] = std::sqrt(var);
        row_valid[gene] = 1;
        orient_row_means[gene] = mean;
        orient_row_sds[gene] = row_sds[gene];
        orient_row_valid[gene] = 1;
      }
    }
  }
  std::vector<int>().swap(row_counts);

  NumericMatrix scores(n_cells, n_sets);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    const std::vector<int>& set = sets[set_i];
    std::vector<int> valid_genes;
    valid_genes.reserve(set.size());
    for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
      if (row_valid[*it]) {
        valid_genes.push_back(*it);
      }
    }

    const int effective_size = static_cast<int>(valid_genes.size());
    if (effective_size < min_size || effective_size > max_size || n_cells <= 1) {
      for (int cell = 0; cell < n_cells; ++cell) {
        scores(cell, set_i) = R_NaN;
      }
      continue;
    }

    arma::mat z(
      static_cast<arma::uword>(effective_size),
      static_cast<arma::uword>(n_cells),
      arma::fill::zeros
    );
    // Match the dense expression view used by GSVA::plageParam in the R
    // benchmark. Structural zeros contribute (0 - mean) / sd.
    for (int row = 0; row < effective_size; ++row) {
      const int gene = valid_genes[row];
      z.row(static_cast<arma::uword>(row)).fill(-row_means[gene] / row_sds[gene]);
      const int row_start = row_ptr[gene];
      const int row_end = row_ptr[gene + 1];
      for (int ptr = row_start; ptr < row_end; ++ptr) {
        const int cell = row_cells[ptr];
        z(static_cast<arma::uword>(row), static_cast<arma::uword>(cell)) +=
          row_values[ptr] / row_sds[gene];
      }
    }

    arma::vec first_v;
    bool ok = false;
    if (effective_size < n_cells && effective_size <= 1024) {
      // Gene×gene covariance path: when gene set is small relative to cell count,
      // eigendecompose ZZ' (effective_size × effective_size) instead of Z'Z (n_cells × n_cells).
      // This is O(effective_size³ + effective_size² × n_cells) vs O(n_cells³).
      arma::mat gene_cov = z * z.t();
      arma::vec eigval;
      arma::mat eigvec;
      ok = arma::eig_sym(eigval, eigvec, gene_cov);
      if (ok && eigvec.n_cols > 0) {
        // u is the first left singular vector (eigenvector of ZZ')
        arma::vec u = eigvec.col(eigvec.n_cols - 1);
        // Convert to right singular vector: v = Z'u, then normalize
        first_v = z.t() * u;
        double norm_val = arma::norm(first_v, 2);
        if (norm_val > 0.0) {
          first_v /= norm_val;
        }
      }
    } else if (n_cells <= 256) {
      arma::mat crossprod = z.t() * z;
      arma::vec eigval;
      arma::mat eigvec;
      ok = arma::eig_sym(eigval, eigvec, crossprod);
      if (ok && eigvec.n_cols > 0) {
        first_v = eigvec.col(eigvec.n_cols - 1);
      }
    } else {
      arma::mat u;
      arma::vec s;
      arma::mat v;
      ok = arma::svd_econ(
        u,
        s,
        v,
        z,
        "right",
        "dc"
      );
      if (ok && v.n_cols > 0) {
        first_v = v.col(0);
      }
    }

    if (!ok || first_v.n_elem == 0) {
      for (int cell = 0; cell < n_cells; ++cell) {
        scores(cell, set_i) = R_NaN;
      }
      continue;
    }

    double direction = 1.0;
    double dot = 0.0;
    std::vector<double> mean_z_by_cell(n_cells, 0.0);
    int orient_size = 0;
    for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
      const int gene = *it;
      if (!row_valid[gene]) {
        continue;
      }
      ++orient_size;
      const int row_start = row_ptr[gene];
      const int row_end = row_ptr[gene + 1];
      for (int ptr = row_start; ptr < row_end; ++ptr) {
        mean_z_by_cell[row_cells[ptr]] +=
          (row_values[ptr] - row_means[gene]) / row_sds[gene];
      }
    }
    if (orient_size > 0) {
      for (int cell = 0; cell < n_cells; ++cell) {
        const double mean_z = mean_z_by_cell[cell] / static_cast<double>(orient_size);
        dot += first_v(static_cast<arma::uword>(cell)) * mean_z;
      }
    }
    if (R_finite(dot) && dot < 0.0) {
      direction = -1.0;
    }

    for (int cell = 0; cell < n_cells; ++cell) {
      scores(cell, set_i) = direction * first_v(static_cast<arma::uword>(cell));
    }
  }

  return scores;
}

static double gsva_sample_sd_from_freq(
  const std::map<double, int>& freq,
  int n
) {
  if (n < 2) {
    return R_NaN;
  }

  long double sum = 0.0;
  for (std::map<double, int>::const_iterator it = freq.begin(); it != freq.end(); ++it) {
    sum += static_cast<long double>(it->first) * static_cast<long double>(it->second);
  }

  long double mean = sum / static_cast<long double>(n);
  if (R_finite(static_cast<double>(mean))) {
    sum = 0.0;
    for (std::map<double, int>::const_iterator it = freq.begin(); it != freq.end(); ++it) {
      sum += (static_cast<long double>(it->first) - mean) * static_cast<long double>(it->second);
    }
    mean += sum / static_cast<long double>(n);
  }

  sum = 0.0;
  for (std::map<double, int>::const_iterator it = freq.begin(); it != freq.end(); ++it) {
    const long double diff = static_cast<long double>(it->first) - mean;
    sum += diff * diff * static_cast<long double>(it->second);
  }

  return std::sqrt(static_cast<double>(sum / static_cast<long double>(n - 1)));
}

struct GsvaRankEntry {
  int gene;
  double value;
};

static bool gsva_rank_entry_before(const GsvaRankEntry& a, const GsvaRankEntry& b) {
  if (a.value < b.value) {
    return true;
  }
  if (a.value > b.value) {
    return false;
  }
  return a.gene > b.gene;
}

static double gsva_sample_sd_from_values(
  const std::vector<double>& values,
  int start,
  int end
) {
  const int n = end - start;
  if (n < 2) {
    return R_NaN;
  }

  long double mean = 0.0;
  for (int i = start; i < end; ++i) {
    mean += static_cast<long double>(values[i]);
  }
  mean /= static_cast<long double>(n);

  long double ss = 0.0;
  for (int i = start; i < end; ++i) {
    const long double diff = static_cast<long double>(values[i]) - mean;
    ss += diff * diff;
  }
  return std::sqrt(static_cast<double>(ss / static_cast<long double>(n - 1)));
}

static std::vector<std::vector<int> > gsva_sets_from_list(List gene_sets, int n_genes) {
  const int n_sets = gene_sets.size();
  std::vector<std::vector<int> > sets(n_sets);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
  }
  return sets;
}

static void gsva_sparse_rows_from_dgc(
  S4 expr,
  int n_genes,
  int n_cells,
  std::vector<int>& row_ptr,
  std::vector<int>& row_cells,
  std::vector<double>& row_values
) {
  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<int> row_counts(n_genes, 0);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        ++row_counts[gene];
      }
    }
  }

  row_ptr.assign(n_genes + 1, 0);
  for (int gene = 0; gene < n_genes; ++gene) {
    row_ptr[gene + 1] = row_ptr[gene] + row_counts[gene];
  }
  const int nnz = row_ptr[n_genes];
  row_cells.assign(nnz, 0);
  row_values.assign(nnz, 0.0);

  std::vector<int> row_cursor(row_ptr);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        const int out = row_cursor[gene]++;
        row_cells[out] = cell;
        row_values[out] = value;
      }
    }
  }
}

// Restored from commit c9324f14 (paper-benchmark version): pairwise Gaussian
// KDE CDF per gene on non-zero values.  This matches the GSVA R package's
// density() → ecdf() pipeline well enough to give Spearman cor ≈ 0.81.
static std::vector<double> gsva_sparse_kcdf_values(
  const std::vector<int>& row_ptr,
  const std::vector<double>& row_values,
  int n_genes,
  int n_cells,
  bool gaussian
) {
  (void) n_cells;  // kept for API compatibility
  std::vector<double> row_kcdf(row_values.size(), R_NaN);

  for (int gene = 0; gene < n_genes; ++gene) {
    const int row_start = row_ptr[gene];
    const int row_end = row_ptr[gene + 1];
    const int nv = row_end - row_start;
    if (nv <= 0) {
      continue;
    }

    double bw = gaussian ? gsva_sample_sd_from_values(row_values, row_start, row_end) / 4.0 : 0.5;
    if (!R_finite(bw) || bw == 0.0) {
      bw = 0.001;
    }

    std::map<double, double> cdf_by_value;
    for (int y_i = row_start; y_i < row_end; ++y_i) {
      const double y = row_values[y_i];
      if (cdf_by_value.find(y) != cdf_by_value.end()) {
        row_kcdf[y_i] = cdf_by_value[y];
        continue;
      }
      double left_tail = 0.0;
      for (int x_i = row_start; x_i < row_end; ++x_i) {
        const double x = row_values[x_i];
        left_tail += gaussian ?
          R::pnorm(y - x, 0.0, bw, true, false) :
          R::ppois(y, x + bw, true, false);
      }
      left_tail /= static_cast<double>(nv);
      cdf_by_value[y] = left_tail;
      row_kcdf[y_i] = left_tail;
    }
  }

  return row_kcdf;
}

static NumericMatrix gsva_score_sparse_kcdf(
  int n_genes,
  int n_cells,
  const std::vector<std::vector<int> >& sets,
  const std::vector<int>& row_ptr,
  const std::vector<int>& row_cells,
  const std::vector<double>& row_kcdf,
  bool max_diff,
  bool abs_ranking,
  double tau
) {
  NumericMatrix scores(n_cells, sets.size());

  std::vector<int> col_counts(n_cells, 0);
  for (std::size_t ptr = 0; ptr < row_cells.size(); ++ptr) {
    ++col_counts[row_cells[ptr]];
  }
  std::vector<std::vector<GsvaRankEntry> > columns(n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    columns[cell].reserve(col_counts[cell]);
  }
  for (int gene = 0; gene < n_genes; ++gene) {
    for (int ptr = row_ptr[gene]; ptr < row_ptr[gene + 1]; ++ptr) {
      const double value = row_kcdf[ptr];
      if (R_finite(value)) {
        columns[row_cells[ptr]].push_back(GsvaRankEntry{gene, value});
      }
    }
  }

  std::vector<int> rank_by_gene(n_genes, 0);
  std::vector<int> zero_rank_by_gene(n_genes, 0);
  std::vector<double> weight_by_pos(n_genes, 0.0);
  std::vector<unsigned char> in_pos(n_genes, 0);
  std::vector<int> touched_genes;
  std::vector<int> touched_pos;

  for (int cell = 0; cell < n_cells; ++cell) {
    std::vector<GsvaRankEntry>& entries = columns[cell];
    std::sort(entries.begin(), entries.end(), gsva_rank_entry_before);
    touched_genes.clear();
    touched_genes.reserve(entries.size());
    for (std::size_t i = 0; i < entries.size(); ++i) {
      const int gene = entries[i].gene;
      rank_by_gene[gene] = static_cast<int>(i) + 1;
      touched_genes.push_back(gene);
    }

    const int nonzero_count = static_cast<int>(entries.size());
    const int zero_count = n_genes - nonzero_count;
    if (zero_count > 0) {
      int zero_rank = 0;
      for (int gene = 0; gene < n_genes; ++gene) {
        if (rank_by_gene[gene] == 0) {
          zero_rank_by_gene[gene] = ++zero_rank;
        }
      }
    }
    const double sparse_mid = (static_cast<double>(nonzero_count) + 1.0) / 2.0;
    const double zero_sym_rank = std::fabs(sparse_mid - 1.0);
    const double dense_mid = static_cast<double>(n_genes) / 2.0;

    for (int set_i = 0; set_i < static_cast<int>(sets.size()); ++set_i) {
      const std::vector<int>& set = sets[set_i];
      const int set_size = static_cast<int>(set.size());
      if (set_size <= 0 || set_size >= n_genes) {
        scores(cell, set_i) = R_NaN;
        continue;
      }

      touched_pos.clear();
      touched_pos.reserve(set.size());
      double sum_in = 0.0;
      for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
        const int gene = *it;
        const int rank = rank_by_gene[gene];
        int pos = 0;
        double sym_rank = 0.0;
        if (zero_count > 0) {
          if (rank > 0) {
            pos = nonzero_count - rank + 1;
            sym_rank = std::fabs(sparse_mid - static_cast<double>(rank + 1));
          } else {
            pos = n_genes - zero_rank_by_gene[gene] + 1;
            sym_rank = zero_sym_rank;
          }
        } else {
          pos = n_genes - rank + 1;
          sym_rank = std::fabs(dense_mid - static_cast<double>(rank));
        }
        if (pos < 1 || pos > n_genes) {
          continue;
        }
        const int pos_idx = pos - 1;
        const double weight = (tau == 1.0) ? sym_rank : std::pow(sym_rank, tau);
        in_pos[pos_idx] = 1;
        weight_by_pos[pos_idx] = weight;
        touched_pos.push_back(pos_idx);
        sum_in += weight;
      }

      double walk_pos = 0.0;
      double walk_neg = 0.0;
      double in_cdf = 0.0;
      double out_cdf = 0.0;
      const double out_step = 1.0 / static_cast<double>(n_genes - set_size);
      const bool can_scale_in = sum_in > 0.0;
      for (int pos = 0; pos < n_genes; ++pos) {
        if (in_pos[pos]) {
          if (can_scale_in) {
            in_cdf += weight_by_pos[pos] / sum_in;
          }
        } else {
          out_cdf += out_step;
        }
        const double walk = in_cdf - out_cdf;
        if (walk > walk_pos) {
          walk_pos = walk;
        }
        if (walk < walk_neg) {
          walk_neg = walk;
        }
      }
      if (max_diff) {
        scores(cell, set_i) = abs_ranking ? (walk_pos - walk_neg) : (walk_pos + walk_neg);
      } else {
        scores(cell, set_i) = (walk_pos > std::fabs(walk_neg)) ? walk_pos : walk_neg;
      }

      for (std::vector<int>::const_iterator it = touched_pos.begin(); it != touched_pos.end(); ++it) {
        weight_by_pos[*it] = 0.0;
        in_pos[*it] = 0;
      }
    }

    for (std::vector<int>::const_iterator it = touched_genes.begin(); it != touched_genes.end(); ++it) {
      rank_by_gene[*it] = 0;
    }
    if (zero_count > 0) {
      for (int gene = 0; gene < n_genes; ++gene) {
        if (zero_rank_by_gene[gene] != 0) {
          zero_rank_by_gene[gene] = 0;
        }
      }
    }
  }

  return scores;
}

static NumericMatrix gsva_sparse_exact(
  S4 expr,
  List gene_sets,
  bool gaussian,
  bool max_diff,
  bool abs_ranking,
  double tau
) {
  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];

  std::vector<int> row_ptr;
  std::vector<int> row_cells;
  std::vector<double> row_values;
  gsva_sparse_rows_from_dgc(expr, n_genes, n_cells, row_ptr, row_cells, row_values);
  std::vector<std::vector<int> > sets = gsva_sets_from_list(gene_sets, n_genes);
  std::vector<double> row_kcdf = gsva_sparse_kcdf_values(row_ptr, row_values, n_genes, n_cells, gaussian);

  return gsva_score_sparse_kcdf(
    n_genes,
    n_cells,
    sets,
    row_ptr,
    row_cells,
    row_kcdf,
    max_diff,
    abs_ranking,
    tau
  );
}

static void gsva_score_z_chunk(
  const std::vector<double>& z,
  int chunk_start,
  int chunk_len,
  int n_genes,
  const std::vector<std::vector<int> >& sets,
  bool max_diff,
  bool abs_ranking,
  double tau,
  NumericMatrix& scores
) {
  std::vector<int> order(n_genes);
  std::vector<int> gene_at_desc_pos(n_genes);
  std::vector<double> sym_rank_stat(n_genes);
  std::vector<int> in_set(n_genes, 0);

  for (int local_cell = 0; local_cell < chunk_len; ++local_cell) {
    std::iota(order.begin(), order.end(), 0);
    const std::size_t col_offset = static_cast<std::size_t>(local_cell) * n_genes;
    std::sort(
      order.begin(),
      order.end(),
      [&](int a, int b) {
        const double za = z[col_offset + a];
        const double zb = z[col_offset + b];
        if (za < zb) {
          return true;
        }
        if (za > zb) {
          return false;
        }
        return a > b;
      }
    );

    for (int rank_i = 0; rank_i < n_genes; ++rank_i) {
      const int gene = order[rank_i];
      const int rank = rank_i + 1;
      const int desc_pos = n_genes - rank;
      gene_at_desc_pos[desc_pos] = gene;
      sym_rank_stat[gene] = std::fabs(static_cast<double>(n_genes) / 2.0 - static_cast<double>(rank));
    }

    const int cell = chunk_start + local_cell;
    for (int set_i = 0; set_i < static_cast<int>(sets.size()); ++set_i) {
      const std::vector<int>& set = sets[set_i];
      const int set_size = static_cast<int>(set.size());
      if (set_size <= 0 || set_size >= n_genes) {
        scores(cell, set_i) = R_NaN;
        continue;
      }

      double sum_in = 0.0;
      for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
        in_set[*it] = 1;
        sum_in += (tau == 1.0) ? sym_rank_stat[*it] : std::pow(sym_rank_stat[*it], tau);
      }

      double walk_pos = 0.0;
      double walk_neg = 0.0;
      double in_cdf = 0.0;
      double out_cdf = 0.0;
      const double out_step = 1.0 / static_cast<double>(n_genes - set_size);
      if (sum_in > 0.0) {
        for (int pos = 0; pos < n_genes; ++pos) {
          const int gene = gene_at_desc_pos[pos];
          if (in_set[gene]) {
            const double weight = (tau == 1.0) ? sym_rank_stat[gene] : std::pow(sym_rank_stat[gene], tau);
            in_cdf += weight / sum_in;
          } else {
            out_cdf += out_step;
          }
          const double walk = in_cdf - out_cdf;
          if (walk > walk_pos) {
            walk_pos = walk;
          }
          if (walk < walk_neg) {
            walk_neg = walk;
          }
        }
        if (max_diff) {
          scores(cell, set_i) = abs_ranking ? (walk_pos - walk_neg) : (walk_pos + walk_neg);
        } else {
          scores(cell, set_i) = (walk_pos > std::fabs(walk_neg)) ? walk_pos : walk_neg;
        }
      } else {
        scores(cell, set_i) = R_NaN;
      }

      for (std::vector<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
        in_set[*it] = 0;
      }
    }
  }
}

static NumericMatrix gsva_score_transformed_rows(
  int n_genes,
  int n_cells,
  const std::vector<std::vector<int> >& sets,
  const std::vector<int>& row_ptr,
  const std::vector<int>& row_cells,
  const std::vector<double>& row_z_values,
  const std::vector<double>& row_z_zero,
  bool max_diff,
  bool abs_ranking,
  double tau,
  int chunk_size
) {
  NumericMatrix scores(n_cells, sets.size());
  const int chunk_n = (chunk_size > 0 && chunk_size < n_cells) ? chunk_size : n_cells;
  std::vector<double> z(static_cast<std::size_t>(n_genes) * chunk_n);

  for (int chunk_start = 0; chunk_start < n_cells; chunk_start += chunk_n) {
    const int chunk_end = std::min(chunk_start + chunk_n, n_cells);
    const int chunk_len = chunk_end - chunk_start;

    for (int local_cell = 0; local_cell < chunk_len; ++local_cell) {
      const std::size_t col_offset = static_cast<std::size_t>(local_cell) * n_genes;
      for (int gene = 0; gene < n_genes; ++gene) {
        z[col_offset + gene] = row_z_zero[gene];
      }
    }

    for (int gene = 0; gene < n_genes; ++gene) {
      const int row_start = row_ptr[gene];
      const int row_end = row_ptr[gene + 1];
      for (int ptr = row_start; ptr < row_end; ++ptr) {
        const int cell = row_cells[ptr];
        if (cell < chunk_start) {
          continue;
        }
        if (cell >= chunk_end) {
          break;
        }
        const int local_cell = cell - chunk_start;
        z[static_cast<std::size_t>(local_cell) * n_genes + gene] = row_z_values[ptr];
      }
    }

    gsva_score_z_chunk(
      z,
      chunk_start,
      chunk_len,
      n_genes,
      sets,
      max_diff,
      abs_ranking,
      tau,
      scores
    );
  }

  return scores;
}

// [[Rcpp::export]]
NumericMatrix gsva_gaussian_dense(
  S4 expr,
  List gene_sets,
  bool max_diff = true,
  bool abs_ranking = false,
  double tau = 1.0,
  int chunk_size = 0
) {
  // Enable z-score KDE path (frequency-based with zeros, log-odds transform,
  // unified gene ranking) — closer to GSVA R's density() → ecdf() → qnorm().
  // The sparse_exact path is kept as fallback in gsva_sparse_exact.
  // return gsva_sparse_exact(
  //   expr,
  //   gene_sets,
  //   true,
  //   max_diff,
  //   abs_ranking,
  //   tau
  // );

  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<int> row_counts(n_genes, 0);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        ++row_counts[gene];
      }
    }
  }

  std::vector<int> row_ptr(n_genes + 1, 0);
  for (int gene = 0; gene < n_genes; ++gene) {
    row_ptr[gene + 1] = row_ptr[gene] + row_counts[gene];
  }
  const int nnz = row_ptr[n_genes];
  std::vector<int> row_cursor(row_ptr);
  std::vector<int> row_cells(nnz);
  std::vector<double> row_values(nnz);

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        const int out = row_cursor[gene]++;
        row_cells[out] = cell;
        row_values[out] = value;
      }
    }
  }
  std::vector<int>().swap(row_counts);
  std::vector<int>().swap(row_cursor);

  std::vector<std::vector<int> > sets(n_sets);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
  }

  std::vector<double> row_z_zero(n_genes, R_NaN);
  std::vector<double> row_z_values(nnz, R_NaN);
  for (int gene = 0; gene < n_genes; ++gene) {
    std::map<double, int> freq;
    const int row_start = row_ptr[gene];
    const int row_end = row_ptr[gene + 1];
    freq[0.0] = n_cells - (row_end - row_start);
    for (int ptr = row_start; ptr < row_end; ++ptr) {
      ++freq[row_values[ptr]];
    }

    double bw = gsva_sample_sd_from_freq(freq, n_cells) / 4.0;
    if (!R_finite(bw) || bw == 0.0) {
      bw = 0.001;
    }

    std::map<double, double> z_by_value;
    for (std::map<double, int>::const_iterator y_it = freq.begin(); y_it != freq.end(); ++y_it) {
      const double y = y_it->first;
      double left_tail = 0.0;
      for (std::map<double, int>::const_iterator x_it = freq.begin(); x_it != freq.end(); ++x_it) {
        left_tail += static_cast<double>(x_it->second) *
          R::pnorm(y - x_it->first, 0.0, bw, true, false);
      }
      left_tail /= static_cast<double>(n_cells);
      z_by_value[y] = -std::log((1.0 - left_tail) / left_tail);
    }

    row_z_zero[gene] = z_by_value[0.0];
    for (int ptr = row_start; ptr < row_end; ++ptr) {
      row_z_values[ptr] = z_by_value[row_values[ptr]];
    }
  }
  std::vector<double>().swap(row_values);

  return gsva_score_transformed_rows(
    n_genes,
    n_cells,
    sets,
    row_ptr,
    row_cells,
    row_z_values,
    row_z_zero,
    max_diff,
    abs_ranking,
    tau,
    chunk_size
  );
}

// [[Rcpp::export]]
NumericMatrix gsva_poisson_dense(
  S4 expr,
  List gene_sets,
  bool max_diff = true,
  bool abs_ranking = false,
  double tau = 1.0,
  int chunk_size = 0
) {
  // Enable z-score KDE path (frequency-based with zeros, log-odds transform,
  // unified gene ranking) — same rationale as gsva_gaussian_dense.
  // return gsva_sparse_exact(
  //   expr,
  //   gene_sets,
  //   false,
  //   max_diff,
  //   abs_ranking,
  //   tau
  // );

  IntegerVector dims = expr.slot("Dim");
  const int n_genes = dims[0];
  const int n_cells = dims[1];
  const int n_sets = gene_sets.size();

  IntegerVector row_idx = expr.slot("i");
  IntegerVector col_ptr = expr.slot("p");
  NumericVector values = expr.slot("x");

  std::vector<int> row_counts(n_genes, 0);
  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        ++row_counts[gene];
      }
    }
  }

  std::vector<int> row_ptr(n_genes + 1, 0);
  for (int gene = 0; gene < n_genes; ++gene) {
    row_ptr[gene + 1] = row_ptr[gene] + row_counts[gene];
  }
  const int nnz = row_ptr[n_genes];
  std::vector<int> row_cursor(row_ptr);
  std::vector<int> row_cells(nnz);
  std::vector<double> row_values(nnz);

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int ptr = col_ptr[cell]; ptr < col_ptr[cell + 1]; ++ptr) {
      const int gene = row_idx[ptr];
      const double value = values[ptr];
      if (R_finite(value) && value != 0.0) {
        const int out = row_cursor[gene]++;
        row_cells[out] = cell;
        row_values[out] = value;
      }
    }
  }
  std::vector<int>().swap(row_counts);
  std::vector<int>().swap(row_cursor);

  std::vector<std::vector<int> > sets(n_sets);
  for (int set_i = 0; set_i < n_sets; ++set_i) {
    IntegerVector genes = gene_sets[set_i];
    sets[set_i].reserve(genes.size());
    for (int gene_i = 0; gene_i < genes.size(); ++gene_i) {
      const int gene = genes[gene_i] - 1;
      if (gene >= 0 && gene < n_genes) {
        sets[set_i].push_back(gene);
      }
    }
    std::sort(sets[set_i].begin(), sets[set_i].end());
    sets[set_i].erase(std::unique(sets[set_i].begin(), sets[set_i].end()), sets[set_i].end());
  }

  std::vector<double> row_z_zero(n_genes, R_NaN);
  std::vector<double> row_z_values(nnz, R_NaN);
  for (int gene = 0; gene < n_genes; ++gene) {
    std::map<double, int> freq;
    const int row_start = row_ptr[gene];
    const int row_end = row_ptr[gene + 1];
    freq[0.0] = n_cells - (row_end - row_start);
    for (int ptr = row_start; ptr < row_end; ++ptr) {
      ++freq[row_values[ptr]];
    }

    std::map<double, double> z_by_value;
    for (std::map<double, int>::const_iterator y_it = freq.begin(); y_it != freq.end(); ++y_it) {
      const double y = y_it->first;
      double left_tail = 0.0;
      for (std::map<double, int>::const_iterator x_it = freq.begin(); x_it != freq.end(); ++x_it) {
        left_tail += static_cast<double>(x_it->second) * R::ppois(y, x_it->first + 0.5, true, false);
      }
      left_tail /= static_cast<double>(n_cells);
      z_by_value[y] = -std::log((1.0 - left_tail) / left_tail);
    }

    row_z_zero[gene] = z_by_value[0.0];
    for (int ptr = row_start; ptr < row_end; ++ptr) {
      row_z_values[ptr] = z_by_value[row_values[ptr]];
    }
  }
  std::vector<double>().swap(row_values);

  return gsva_score_transformed_rows(
    n_genes,
    n_cells,
    sets,
    row_ptr,
    row_cells,
    row_z_values,
    row_z_zero,
    max_diff,
    abs_ranking,
    tau,
    chunk_size
  );
}
