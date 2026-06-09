// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <random>
#include <utility>
#include <thread>
#include <atomic>

using namespace Rcpp;

int cytotrace2_worker_count(int requested, int n_tasks) {
  if (n_tasks <= 1) return 1;
  int n_workers = std::max(1, std::min(requested, n_tasks));
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    n_workers = std::min(n_workers, static_cast<int>(hardware));
  }
  return std::max(1, n_workers);
}

// ============================================================================
// Utility functions
// ============================================================================

// [[Rcpp::export]]
List cytotrace2_preprocess_numeric(const arma::mat& expression_mapped) {
  int n_genes = expression_mapped.n_rows;
  int n_cells = expression_mapped.n_cols;

  arma::mat ranked_data(n_cells, n_genes);
  arma::mat log2_data(n_cells, n_genes);
  int count_cells_few_genes = 0;

  std::vector<std::pair<double, int>> values;
  values.reserve(n_genes);

  for (int cell = 0; cell < n_cells; cell++) {
    values.clear();
    double col_sum = 0.0;
    int expressed = 0;
    for (int gene = 0; gene < n_genes; gene++) {
      double value = expression_mapped(gene, cell);
      values.push_back(std::make_pair(value, gene));
      col_sum += value;
      if (value > 0.0) {
        expressed++;
      }
    }
    if (expressed < 500) {
      count_cells_few_genes++;
    }

    std::sort(
      values.begin(),
      values.end(),
      [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        if (a.first > b.first) return true;
        if (a.first < b.first) return false;
        return a.second < b.second;
      }
    );

    int start = 0;
    while (start < n_genes) {
      int end = start + 1;
      while (end < n_genes && values[end].first == values[start].first) {
        end++;
      }
      double avg_rank = (static_cast<double>(start + 1) + static_cast<double>(end)) / 2.0;
      for (int pos = start; pos < end; pos++) {
        ranked_data(cell, values[pos].second) = avg_rank;
      }
      start = end;
    }

    if (col_sum > 0.0) {
      double scale = 1000000.0 / col_sum;
      for (int gene = 0; gene < n_genes; gene++) {
        log2_data(cell, gene) = std::log2(expression_mapped(gene, cell) * scale + 1.0);
      }
    } else {
      log2_data.row(cell).zeros();
    }
  }

  return List::create(
    Named("ranked_data") = ranked_data,
    Named("log2_data") = log2_data,
    Named("count_cells_few_genes") = count_cells_few_genes
  );
}

// Softmax over columns (each row independently)
inline arma::mat softmax_rows(const arma::mat& X) {
  arma::mat Y = arma::exp(X);
  arma::vec rowsum = arma::sum(Y, 1);
  Y.each_col() /= rowsum;
  return Y;
}

// ============================================================================
// Stage 2: Binary Module forward pass (single layer)
// ============================================================================

// Pre-extracted layer parameters (pure C++ struct, no R API dependency)
struct LayerParams {
  arma::sp_mat weight;
  arma::rowvec n_vec;
  double maxrank_norm;
  arma::mat gs_backgrounds;
  arma::rowvec scale_factors;
  arma::rowvec running_mean;
  arma::rowvec running_var;
  arma::rowvec out_weight;
  double out_bias;
};

// Extract layer params from R list (called BEFORE parallel region)
LayerParams extract_layer_params(const List& layer_dict) {
  LayerParams p;
  p.weight = as<arma::sp_mat>(layer_dict["weight"]);
  p.n_vec = as<arma::rowvec>(layer_dict["n"]);
  p.maxrank_norm = as<double>(layer_dict["maxrank_norm"]);
  p.gs_backgrounds = arma::mat(as<arma::sp_mat>(layer_dict["gs.backgrounds"]));
  p.scale_factors = as<arma::rowvec>(layer_dict["scale.factors"]);
  p.running_mean = as<arma::rowvec>(layer_dict["running_mean"]);
  p.running_var = as<arma::rowvec>(layer_dict["running_var"]);
  p.out_weight = as<arma::rowvec>(layer_dict["out.weight"]);
  p.out_bias = as<double>(layer_dict["out.bias"]);
  return p;
}

// Pure C++ forward pass (no R API calls, thread-safe)
arma::vec binary_module_forward(
    const arma::mat& rank_data,     // [n_cells x n_genes]
    const arma::mat& log2_data,     // [n_cells x n_genes]
    const LayerParams& p,
    double eps = 1e-5
) {
  int n_cells = rank_data.n_rows;
  int hidden_size = p.n_vec.n_elem;

  // --- UCell score ---
  arma::mat R_UCell(n_cells, hidden_size, arma::fill::zeros);
  for (int j = 0; j < hidden_size; j++) {
    for (arma::sp_mat::const_col_iterator it = p.weight.begin_col(j);
         it != p.weight.end_col(j);
         ++it) {
      arma::uword gene = it.row();
      double weight = *it;
      for (int cell = 0; cell < n_cells; cell++) {
        double rank_value = rank_data(cell, gene);
        R_UCell(cell, j) += (rank_value > p.maxrank_norm ? p.maxrank_norm : rank_value) * weight;
      }
    }
  }

  for (int j = 0; j < hidden_size; j++) {
    double nj = p.n_vec(j);
    double n_offset = nj * (nj + 1.0) / 2.0;
    double denom = nj * p.maxrank_norm;
    R_UCell.col(j) = 1.0 - (R_UCell.col(j) - n_offset) / denom;
  }

  // --- AMS score ---
  arma::mat raw_score = log2_data * p.weight;
  for (int j = 0; j < hidden_size; j++) {
    raw_score.col(j) /= p.n_vec(j);
  }

  arma::mat bg_score = log2_data * p.gs_backgrounds;
  for (int j = 0; j < hidden_size; j++) {
    bg_score.col(j) /= p.scale_factors(j);
  }

  arma::mat R_AMS = raw_score - bg_score;

  // --- Batch Normalization + Linear output ---
  arma::mat R_all = arma::join_rows(R_UCell, R_AMS);
  for (int j = 0; j < 2 * hidden_size; j++) {
    double inv_std = 1.0 / std::sqrt(p.running_var(j) + eps);
    R_all.col(j) = (R_all.col(j) - p.running_mean(j)) * inv_std;
  }

  return R_all * p.out_weight.t() + p.out_bias;
}

// ============================================================================
// Stage 2 (cont.): Binary Encoder (6 layers) -> one model
// ============================================================================

// Pre-extracted model: vector of 6 LayerParams
typedef std::vector<LayerParams> ModelParams;

void extract_model_params(const List& model_dict, ModelParams& mp) {
  int n_layers = model_dict.size();
  mp.resize(n_layers);
  for (int layer = 0; layer < n_layers; layer++) {
    mp[layer] = extract_layer_params(model_dict[layer]);
  }
}

// Pure C++ binary encoder (no R API calls, thread-safe)
void binary_encoder_forward(
    const arma::mat& rank_data,
    const arma::mat& log2_data,
    const ModelParams& model_params,
    arma::mat& prob_pred_out,
    arma::vec& prob_order_out
) {
  int n_layers = model_params.size();
  int n_cells = rank_data.n_rows;

  arma::mat preds(n_cells, n_layers);
  for (int layer = 0; layer < n_layers; layer++) {
    preds.col(layer) = binary_module_forward(rank_data, log2_data, model_params[layer]);
  }

  prob_pred_out = softmax_rows(preds);

  prob_order_out.zeros(n_cells);
  for (int j = 0; j < n_layers; j++) {
    prob_order_out += prob_pred_out.col(j) * (static_cast<double>(j) / 5.0);
  }
}

// ============================================================================
// Stage 2 (cont.): Ensemble Prediction (19 models)
// ============================================================================

List cytotrace2_ensemble_predict(
    const arma::mat& rank_data,      // [n_cells x n_genes]
    const arma::mat& log2_data,      // [n_cells x n_genes]
    const List& parameter_dict,      // list of 19 model dicts
    int cores = 1
) {
  int n_models = parameter_dict.size();
  int n_cells = rank_data.n_rows;
  int n_workers = cytotrace2_worker_count(cores, n_models);

  arma::mat sum_probs(n_cells, 6, arma::fill::zeros);
  arma::vec sum_order(n_cells, arma::fill::zeros);

  if (n_workers == 1) {
    arma::mat prob_pred;
    arma::vec prob_order;
    for (int m = 0; m < n_models; m++) {
      ModelParams mp;
      extract_model_params(parameter_dict[m], mp);
      binary_encoder_forward(rank_data, log2_data, mp, prob_pred, prob_order);
      sum_probs += prob_pred;
      sum_order += prob_order;
    }

    arma::vec predicted_score = sum_order / static_cast<double>(n_models);
    arma::mat avg_probs = sum_probs / static_cast<double>(n_models);

    arma::uvec cat_idx = arma::index_max(avg_probs, 1);
    IntegerVector predicted_category(n_cells);
    for (int i = 0; i < n_cells; i++) {
      predicted_category[i] = static_cast<int>(cat_idx(i)) + 1;
    }

    return List::create(
      Named("score") = predicted_score,
      Named("category") = predicted_category,
      Named("probabilities") = avg_probs
    );
  }

  std::vector<ModelParams> model_params(n_models);
  for (int m = 0; m < n_models; m++) {
    extract_model_params(parameter_dict[m], model_params[m]);
  }

  std::vector<arma::mat> model_probs(n_models);
  std::vector<arma::vec> model_orders(n_models);
  std::atomic<int> next_model(0);
  std::vector<std::thread> workers;
  workers.reserve(n_workers);

  for (int worker = 0; worker < n_workers; worker++) {
    workers.emplace_back([&]() {
      while (true) {
        int m = next_model.fetch_add(1);
        if (m >= n_models) break;
        binary_encoder_forward(
          rank_data, log2_data, model_params[m], model_probs[m], model_orders[m]
        );
      }
    });
  }

  for (std::thread& worker : workers) {
    worker.join();
  }

  for (int m = 0; m < n_models; m++) {
    sum_probs += model_probs[m];
    sum_order += model_orders[m];
  }

  arma::vec predicted_score = sum_order / static_cast<double>(n_models);
  arma::mat avg_probs = sum_probs / static_cast<double>(n_models);

  arma::uvec cat_idx = arma::index_max(avg_probs, 1);
  IntegerVector predicted_category(n_cells);
  for (int i = 0; i < n_cells; i++) {
    predicted_category[i] = static_cast<int>(cat_idx(i)) + 1;
  }

  return List::create(
    Named("score") = predicted_score,
    Named("category") = predicted_category,
    Named("probabilities") = avg_probs
  );
}

// ============================================================================
// Stage 3: Diffusion Smoothing
// ============================================================================

// Compute dispersion (var/mean) for each gene column
arma::vec compute_dispersion(const arma::mat& X) {
  int n_genes = X.n_cols;
  arma::vec disp(n_genes);
  for (int j = 0; j < n_genes; j++) {
    double m = arma::mean(X.col(j));
    if (m == 0.0) {
      disp(j) = 0.0;
    } else {
      double v = arma::var(X.col(j));
      disp(j) = v / m;
    }
  }
  return disp;
}

// Build Markov transition matrix from Pearson correlation
arma::mat build_markov_matrix(const arma::mat& sub_mat) {
  // sub_mat: [n_cells x n_top_genes]
  int n_cells = sub_mat.n_rows;

  // Pearson correlation matrix
  arma::mat D = arma::cor(sub_mat.t());  // [n_cells x n_cells]

  // Zero out diagonal
  D.diag().zeros();

  // Replace NaN with 0
  D.replace(arma::datum::nan, 0.0);

  // Threshold at max(mean(D), 0)
  double cutoff = std::max(arma::mean(arma::mean(D)), 0.0);
  D.elem(arma::find(D < cutoff)).zeros();

  // Row-normalize to Markov matrix
  arma::vec row_sums = arma::sum(D, 1) + 1e-5;
  arma::mat A = D.each_col() / row_sums;

  return A;
}

// Iterative smoothing with restart
arma::vec smoothing_by_diffusion(
    const arma::vec& init_score,
    const arma::mat& markov_mat,
    int maxiter = 10000
) {
  arma::vec prev_score = init_score;
  double init_mean = arma::mean(init_score) + 1e-6;

  for (int iter = 0; iter < maxiter; iter++) {
    arma::vec cur_score = 0.9 * (markov_mat * prev_score) + 0.1 * init_score;

    double change = arma::mean(arma::abs(cur_score - prev_score)) / init_mean;
    if (change < 1e-6) {
      return cur_score;
    }
    prev_score = cur_score;
  }

  return prev_score;
}

arma::vec cytotrace2_diffusion_smooth(
  const arma::mat& log2_data,      // [n_cells x n_genes]
  const arma::vec& raw_scores,     // [n_cells]
  const IntegerVector& smooth_groups
) {
  int n_cells = log2_data.n_rows;
  int n_genes = log2_data.n_cols;
  if (smooth_groups.size() != n_cells) {
    stop("smooth_groups must have one entry per cell");
  }

  // Select top 1000 most variable genes by dispersion
  arma::vec dispersion = compute_dispersion(log2_data);
  int n_top = std::min(1000, n_genes);
  arma::uvec top_idx = arma::sort_index(dispersion, "descend");
  top_idx = top_idx.head(n_top);

  arma::mat sub_mat = log2_data.cols(top_idx);

  arma::vec smoothed_scores(n_cells);

  IntegerVector group_levels = sort_unique(smooth_groups);
  for (int chunk = 0; chunk < group_levels.size(); chunk++) {
    int group = group_levels[chunk];
    std::vector<arma::uword> idx_vec;
    idx_vec.reserve(n_cells);
    for (int i = 0; i < n_cells; i++) {
      if (smooth_groups[i] == group) {
        idx_vec.push_back(static_cast<arma::uword>(i));
      }
    }
    if (idx_vec.empty()) continue;

    arma::uvec chunk_idx(idx_vec);
    arma::mat chunk_log2 = sub_mat.rows(chunk_idx);
    arma::vec chunk_scores = raw_scores(chunk_idx);

    arma::mat markov = build_markov_matrix(chunk_log2);
    arma::vec smoothed = smoothing_by_diffusion(chunk_scores, markov);

    for (int i = 0; i < static_cast<int>(chunk_idx.n_elem); i++) {
      smoothed_scores(chunk_idx(i)) = smoothed(i);
    }
  }

  return smoothed_scores;
}

// ============================================================================
// Stage 4: Binning
// ============================================================================

List cytotrace2_bin_data(
    const arma::vec& smoothed_scores,    // [n_cells]
    const IntegerVector& categories,     // [n_cells] values 1-6
    const StringVector& category_labels  // length 6
) {
  int n_cells = smoothed_scores.n_elem;
  arma::vec binned_scores(n_cells);
  StringVector binned_categories(n_cells);

  // Unit interval bounds for each of 6 categories
  arma::vec limits = arma::linspace(0.0, 1.0, 7);

  for (int cat = 0; cat < 6; cat++) {
    double lower = limits(cat);
    double upper = limits(cat + 1);

    // Find cells in this category (categories are 1-based)
    std::vector<int> cell_indices;
    std::vector<double> cell_scores;
    for (int i = 0; i < n_cells; i++) {
      if (categories[i] == cat + 1) {
        cell_indices.push_back(i);
        cell_scores.push_back(smoothed_scores(i));
      }
    }

    int n_in_cat = cell_indices.size();
    if (n_in_cat == 0) continue;

    // Sort cells by score and get ranks
    std::vector<std::pair<double, int>> score_idx;
    for (int i = 0; i < n_in_cat; i++) {
      score_idx.push_back({cell_scores[i], cell_indices[i]});
    }
    std::sort(score_idx.begin(), score_idx.end());

    double feat_low = lower + 1e-8;
    double feat_high = upper - 1e-8;

    for (int i = 0; i < n_in_cat; i++) {
      int cell_i = score_idx[i].second;
      if (n_in_cat == 1) {
        binned_scores(cell_i) = (feat_low + feat_high) / 2.0;
      } else {
        double frac = static_cast<double>(i) / static_cast<double>(n_in_cat - 1);
        binned_scores(cell_i) = feat_low + frac * (feat_high - feat_low);
      }
      binned_categories[cell_i] = category_labels[cat];
    }
  }

  return List::create(
    Named("preKNN_score") = binned_scores,
    Named("preKNN_potency") = binned_categories
  );
}

// ============================================================================
// Stage 5: Adaptive kNN Smoothing
// ============================================================================

// Maps a score to potency category label index (0-5)
int map_score_to_potency(double score) {
  if (score <= 1.0 / 6.0) return 0;      // Differentiated
  else if (score <= 2.0 / 6.0) return 1; // Unipotent
  else if (score <= 3.0 / 6.0) return 2; // Oligopotent
  else if (score <= 4.0 / 6.0) return 3; // Multipotent
  else if (score <= 5.0 / 6.0) return 4; // Pluripotent
  else return 5;                          // Totipotent
}

// Find shortest consensus segment
int shortest_consensus(const arma::vec& neighbor_scores) {
  int len = neighbor_scores.n_elem;
  int idx_use = 2;
  bool last_part = false;
  int max_i = len / 2;

  for (int i = 2; i <= max_i; i++) {
    double mean1 = arma::mean(neighbor_scores.subvec(0, i - 1));
    double mean2 = arma::mean(neighbor_scores.subvec(i, std::min(2 * i, len) - 1));

    int pot1 = map_score_to_potency(mean1);
    int pot2 = map_score_to_potency(mean2);

    if (pot1 == pot2 && !last_part) {
      idx_use = i;
      last_part = true;
    }
  }

  return 2 * idx_use;
}

List cytotrace2_knn_smooth(
    const arma::mat& pca_coords,         // [n_cells x n_pcs]
    const arma::vec& preKNN_scores,      // [n_cells]
    const StringVector& preKNN_potency,  // [n_cells]
    int cores = 1,
    int seed = 14
) {
  int n_cells = pca_coords.n_rows;

  if (n_cells < 100) {
    // Skip kNN smoothing for small datasets
    StringVector final_potency(n_cells);
    for (int i = 0; i < n_cells; i++) {
      final_potency[i] = preKNN_potency[i];
    }
    return List::create(
      Named("CytoTRACE2_Score") = preKNN_scores,
      Named("CytoTRACE2_Potency") = final_potency
    );
  }

  // kNN smoothing
  arma::vec final_scores(n_cells);

  // Pre-compute row squared norms for vectorized distance computation
  // ||x_i - x_j||² = ||x_i||² + ||x_j||² - 2 * x_i · x_j
  arma::vec row_norms_sq = arma::sum(arma::square(pca_coords), 1);

  for (int i = 0; i < n_cells; i++) {
    // Vectorized Euclidean distances to all other cells in PCA space
    arma::vec dists_sq = row_norms_sq + row_norms_sq(i) -
      2.0 * (pca_coords * pca_coords.row(i).t());
    // Clamp tiny negative values from floating-point rounding
    dists_sq.elem(arma::find(dists_sq < 0.0)).zeros();
    arma::vec dists = arma::sqrt(dists_sq);

    // Normalize distances by max
    double max_dist = dists.max();
    if (max_dist > 0) {
      dists /= max_dist;
    }

    int n_neighbors = std::min(30, n_cells);
    std::vector<std::pair<double, int>> distance_idx;
    distance_idx.reserve(n_cells);
    for (int j = 0; j < n_cells; j++) {
      distance_idx.push_back(std::make_pair(dists(j), j));
    }
    std::partial_sort(
      distance_idx.begin(),
      distance_idx.begin() + n_neighbors,
      distance_idx.end(),
      [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        if (a.first < b.first) return true;
        if (a.first > b.first) return false;
        return a.second < b.second;
      }
    );

    // Extract neighbor scores (sorted by distance)
    arma::vec neighbor_scores(n_neighbors);
    arma::vec neighbor_dists(n_neighbors);
    for (int k = 0; k < n_neighbors; k++) {
      int neighbor_idx = distance_idx[k].second;
      neighbor_scores(k) = preKNN_scores(neighbor_idx);
      neighbor_dists(k) = dists(neighbor_idx);
    }

    int num_neighbors_keep = shortest_consensus(neighbor_scores);

    if (num_neighbors_keep > 1) {
      // Distance-weighted average
      double weight_sum = 0.0;
      double score_sum = 0.0;
      for (int k = 0; k < num_neighbors_keep; k++) {
        double w = std::pow(1.0 - neighbor_dists(k), 2.0);
        weight_sum += w;
        score_sum += w * neighbor_scores(k);
      }
      final_scores(i) = score_sum / weight_sum;
    } else {
      final_scores(i) = preKNN_scores(i);
    }
  }

  // Map scores to potency categories
  StringVector final_potency(n_cells);
  StringVector labels = StringVector::create(
    "Differentiated", "Unipotent", "Oligopotent",
    "Multipotent", "Pluripotent", "Totipotent"
  );

  for (int i = 0; i < n_cells; i++) {
    int cat_idx = map_score_to_potency(final_scores(i));
    final_potency[i] = labels[cat_idx];
  }

  return List::create(
    Named("CytoTRACE2_Score") = final_scores,
    Named("CytoTRACE2_Potency") = final_potency
  );
}

// ============================================================================
// Main entry point: full pipeline in C++
// ============================================================================

// [[Rcpp::export]]
List cytotrace2_main(
    const arma::mat& rank_data,          // [n_cells x n_genes]
    const arma::mat& log2_data,          // [n_cells x n_genes]
    const List& parameter_dict,          // 19 model dicts
    const IntegerVector& smooth_groups,
    int cores,
    int seed,
    const arma::mat& pca_coords              // [n_cells x n_pcs] pre-computed PCA
) {
  int n_cells = rank_data.n_rows;

  // Stage 2: Ensemble prediction
  List pred_result = cytotrace2_ensemble_predict(rank_data, log2_data, parameter_dict, cores);
  arma::vec raw_scores = as<arma::vec>(pred_result["score"]);
  IntegerVector raw_categories = as<IntegerVector>(pred_result["category"]);

  // Map integer categories to labels
  StringVector cat_labels = StringVector::create(
    "Differentiated", "Unipotent", "Oligopotent",
    "Multipotent", "Pluripotent", "Totipotent"
  );
  StringVector preKNN_potency(n_cells);
  for (int i = 0; i < n_cells; i++) {
    preKNN_potency[i] = cat_labels[raw_categories[i] - 1];
  }

  // Stage 3: Diffusion smoothing
  arma::vec smoothed = cytotrace2_diffusion_smooth(
    log2_data, raw_scores, smooth_groups
  );

  if (n_cells <= 10) {
    return List::create(
      Named("CytoTRACE2_Score") = smoothed,
      Named("CytoTRACE2_Potency") = preKNN_potency,
      Named("CytoTRACE2_Relative") = arma::vec(n_cells, arma::fill::zeros),
      Named("preKNN_CytoTRACE2_Score") = smoothed,
      Named("preKNN_CytoTRACE2_Potency") = preKNN_potency
    );
  }

  // Stage 4: Binning
  List bin_result = cytotrace2_bin_data(smoothed, raw_categories, cat_labels);
  arma::vec binned_scores = as<arma::vec>(bin_result["preKNN_score"]);
  StringVector binned_potency = as<StringVector>(bin_result["preKNN_potency"]);

  // Stage 5: kNN smoothing
  List knn_result;
  if (n_cells > 100 && pca_coords.n_rows == n_cells && pca_coords.n_cols > 0) {
    knn_result = cytotrace2_knn_smooth(pca_coords, binned_scores, binned_potency, cores, seed);
  } else if (n_cells > 100) {
    // Fallback: skip kNN if PCA not provided
    StringVector final_potency = binned_potency;
    knn_result = List::create(
      Named("CytoTRACE2_Score") = binned_scores,
      Named("CytoTRACE2_Potency") = final_potency
    );
  } else {
    knn_result = List::create(
      Named("CytoTRACE2_Score") = binned_scores,
      Named("CytoTRACE2_Potency") = binned_potency
    );
  }

  arma::vec final_scores = as<arma::vec>(knn_result["CytoTRACE2_Score"]);
  StringVector final_potency = as<StringVector>(knn_result["CytoTRACE2_Potency"]);

  // Compute relative score
  double min_score = final_scores.min();
  double max_score = final_scores.max();
  double range = max_score - min_score;
  arma::vec relative_scores(n_cells);
  if (range > 0) {
    relative_scores = (final_scores - min_score) / range;
  } else {
    relative_scores.zeros();
  }

  return List::create(
    Named("CytoTRACE2_Score") = final_scores,
    Named("CytoTRACE2_Potency") = final_potency,
    Named("CytoTRACE2_Relative") = relative_scores,
    Named("preKNN_CytoTRACE2_Score") = binned_scores,
    Named("preKNN_CytoTRACE2_Potency") = binned_potency
  );
}
