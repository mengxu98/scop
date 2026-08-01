#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix vector_weighted_arrows_cpp(
    const NumericMatrix& centers,
    const NumericVector& scores,
    const NumericVector& embedding_range,
    double p,
    double arrow_length) {
  const int n = centers.nrow();
  if (centers.ncol() != 2 || scores.size() != n ||
      embedding_range.size() != 4) {
    stop("invalid VECTOR kernel dimensions");
  }
  if (n < 2) {
    return NumericMatrix(0, 8);
  }

  double minimum_distance = R_PosInf;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dx = centers(j, 0) - centers(i, 0);
      const double dy = centers(j, 1) - centers(i, 1);
      const double distance = std::sqrt(dx * dx + dy * dy);
      if (distance > 0.0 && distance < minimum_distance) {
        minimum_distance = distance;
      }
    }
  }
  if (!R_finite(minimum_distance)) {
    return NumericMatrix(0, 8);
  }
  const double one = minimum_distance * arrow_length;

  std::vector<double> values;
  values.reserve(static_cast<std::size_t>(n) * 8);
  std::vector<double> distances(n);
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const double dx = centers(j, 0) - centers(i, 0);
      const double dy = centers(j, 1) - centers(i, 1);
      distances[j] = std::sqrt(dx * dx + dy * dy);
      order[j] = j;
    }
    std::stable_sort(order.begin(), order.end(), [&](int lhs, int rhs) {
      return distances[lhs] < distances[rhs];
    });
    std::vector<int> first_rank(n);
    for (int rank = 0; rank < n; ++rank) {
      first_rank[order[rank]] = rank;
    }

    double denominator = 0.0;
    double numerator_x = 0.0;
    double numerator_y = 0.0;
    for (int j = 0; j < n; ++j) {
      const double score_weight = scores[i] - scores[j];
      const double weight =
          std::pow(p, static_cast<double>(first_rank[j])) * score_weight;
      if (!R_finite(weight)) {
        continue;
      }
      denominator += std::abs(weight);
      if (distances[j] > 0.0 && R_finite(distances[j])) {
        numerator_x +=
            ((centers(j, 0) - centers(i, 0)) / distances[j] * one) * weight;
        numerator_y +=
            ((centers(j, 1) - centers(i, 1)) / distances[j] * one) * weight;
      }
    }
    if (!R_finite(denominator) || denominator == 0.0) {
      continue;
    }
    const double dx = numerator_x / denominator;
    const double dy = numerator_y / denominator;
    const double length = std::sqrt(dx * dx + dy * dy);
    if (!R_finite(dx) || !R_finite(dy) || length == 0.0) {
      continue;
    }
    const double x_end =
        std::max(embedding_range[0],
                 std::min(embedding_range[1], centers(i, 0) + dx));
    const double y_end =
        std::max(embedding_range[2],
                 std::min(embedding_range[3], centers(i, 1) + dy));
    values.push_back(i + 1);
    values.push_back(centers(i, 0));
    values.push_back(centers(i, 1));
    values.push_back(dx);
    values.push_back(dy);
    values.push_back(x_end);
    values.push_back(y_end);
    values.push_back(length);
  }

  NumericMatrix out(values.size() / 8, 8);
  for (int row = 0; row < out.nrow(); ++row) {
    for (int column = 0; column < 8; ++column) {
      out(row, column) = values[static_cast<std::size_t>(row) * 8 + column];
    }
  }
  return out;
}

inline double log_cpm_value(double value, double total) {
  if (!R_finite(value) || !R_finite(total) || total <= 0.0) {
    return 0.0;
  }
  return std::log(value / total * 1e6 + 1.0) / std::log(2.0);
}

// [[Rcpp::export]]
NumericVector cytospace_estimate_fractions_cpp(
    const NumericMatrix& st_expr,
    const NumericMatrix& ref_expr,
    const IntegerVector& labels,
    int n_types,
    const NumericVector& spot_weights) {
  const int genes = st_expr.nrow();
  const int spots = st_expr.ncol();
  const int cells = ref_expr.ncol();
  if (ref_expr.nrow() != genes || labels.size() != cells ||
      spot_weights.size() != spots || n_types < 1) {
    stop("invalid CytoSPACE fraction kernel dimensions");
  }

  std::vector<double> st_totals(spots, 0.0);
  std::vector<double> ref_totals(cells, 0.0);
  for (int j = 0; j < spots; ++j) {
    for (int g = 0; g < genes; ++g) {
      st_totals[j] += st_expr(g, j);
    }
  }
  for (int j = 0; j < cells; ++j) {
    if (labels[j] < 1 || labels[j] > n_types) {
      stop("labels contains an out-of-range cell type code");
    }
    for (int g = 0; g < genes; ++g) {
      ref_totals[j] += ref_expr(g, j);
    }
  }

  std::vector<int> type_counts(n_types, 0);
  NumericMatrix profiles(genes, n_types);
  for (int j = 0; j < cells; ++j) {
    const int type = labels[j] - 1;
    ++type_counts[type];
    for (int g = 0; g < genes; ++g) {
      profiles(g, type) += log_cpm_value(ref_expr(g, j), ref_totals[j]);
    }
  }
  for (int type = 0; type < n_types; ++type) {
    if (type_counts[type] == 0) {
      stop("every cell type must contain at least one reference cell");
    }
    for (int g = 0; g < genes; ++g) {
      profiles(g, type) /= static_cast<double>(type_counts[type]);
    }
  }

  std::vector<double> profile_mean(n_types, 0.0);
  std::vector<double> profile_ss(n_types, 0.0);
  for (int type = 0; type < n_types; ++type) {
    for (int g = 0; g < genes; ++g) {
      profile_mean[type] += profiles(g, type);
    }
    profile_mean[type] /= static_cast<double>(genes);
    for (int g = 0; g < genes; ++g) {
      const double centered = profiles(g, type) - profile_mean[type];
      profile_ss[type] += centered * centered;
    }
  }

  NumericVector fractions(n_types);
  double total_weight = 0.0;
  std::vector<double> spot_values(genes);
  std::vector<double> score(n_types);
  for (int spot = 0; spot < spots; ++spot) {
    double spot_mean = 0.0;
    for (int g = 0; g < genes; ++g) {
      spot_values[g] =
          log_cpm_value(st_expr(g, spot), st_totals[spot]);
      spot_mean += spot_values[g];
    }
    spot_mean /= static_cast<double>(genes);
    double spot_ss = 0.0;
    for (int g = 0; g < genes; ++g) {
      const double centered = spot_values[g] - spot_mean;
      spot_ss += centered * centered;
    }

    double score_sum = 0.0;
    for (int type = 0; type < n_types; ++type) {
      double covariance = 0.0;
      for (int g = 0; g < genes; ++g) {
        covariance += (spot_values[g] - spot_mean) *
                      (profiles(g, type) - profile_mean[type]);
      }
      const double denominator = std::sqrt(spot_ss * profile_ss[type]);
      const double correlation =
          denominator > 0.0 ? covariance / denominator : 0.0;
      score[type] =
          R_finite(correlation) ? std::max(correlation, 0.0) : 0.0;
      score_sum += score[type];
    }
    if (score_sum <= 0.0) {
      std::fill(score.begin(), score.end(), 1.0 / n_types);
    } else {
      for (int type = 0; type < n_types; ++type) {
        score[type] /= score_sum;
      }
    }
    const double weight =
        R_finite(spot_weights[spot]) ? spot_weights[spot] : 0.0;
    total_weight += weight;
    for (int type = 0; type < n_types; ++type) {
      fractions[type] += score[type] * weight;
    }
  }
  if (total_weight <= 0.0) {
    stop("spot_weights must have a positive finite sum");
  }
  for (int type = 0; type < n_types; ++type) {
    fractions[type] /= total_weight;
  }
  return fractions;
}

// [[Rcpp::export]]
DataFrame spatial_pair_count_cpp(
    const CharacterVector& sample,
    const CharacterVector& condition,
    const CharacterVector& subject,
    const CharacterVector& from,
    const CharacterVector& to) {
  const R_xlen_t n = sample.size();
  if (condition.size() != n || subject.size() != n || from.size() != n ||
      to.size() != n) {
    stop("spatial pair columns must have equal length");
  }
  typedef std::tuple<std::string, std::string, std::string, std::string,
                     std::string>
      Key;
  std::map<Key, int> counts;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (CharacterVector::is_na(sample[i]) ||
        CharacterVector::is_na(condition[i]) ||
        CharacterVector::is_na(subject[i]) ||
        CharacterVector::is_na(from[i]) ||
        CharacterVector::is_na(to[i])) {
      continue;
    }
    Key key = std::make_tuple(
        as<std::string>(sample[i]), as<std::string>(condition[i]),
        as<std::string>(subject[i]), as<std::string>(from[i]),
        as<std::string>(to[i]));
    ++counts[key];
  }
  const int size = counts.size();
  CharacterVector sample_out(size), condition_out(size), subject_out(size),
      from_out(size), to_out(size);
  IntegerVector count_out(size);
  int i = 0;
  for (std::map<Key, int>::const_iterator it = counts.begin();
       it != counts.end(); ++it, ++i) {
    sample_out[i] = std::get<0>(it->first);
    condition_out[i] = std::get<1>(it->first);
    subject_out[i] = std::get<2>(it->first);
    from_out[i] = std::get<3>(it->first);
    to_out[i] = std::get<4>(it->first);
    count_out[i] = it->second;
  }
  return DataFrame::create(
      _["sample"] = sample_out,
      _["condition"] = condition_out,
      _["subject"] = subject_out,
      _["from"] = from_out,
      _["to"] = to_out,
      _["count"] = count_out,
      _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
List bbknn_fuzzy_membership_cpp(
    const IntegerMatrix& index,
    const NumericMatrix& distance,
    double local_connectivity = 1.0,
    double bandwidth = 1.0) {
  const int n = index.nrow();
  const int k = index.ncol();
  if (distance.nrow() != n || distance.ncol() != k || n < 2 || k < 2 ||
      local_connectivity < 0.0 || bandwidth <= 0.0) {
    stop("invalid BBKNN fuzzy-membership inputs");
  }
  IntegerMatrix sorted_index(n, k);
  NumericMatrix sorted_distance(n, k);
  NumericMatrix membership(n, k);
  std::vector<std::pair<double, int> > neighbors(k);
  std::vector<double> nonzero;
  nonzero.reserve(k);
  const double target = std::log2(static_cast<double>(k)) * bandwidth;

  for (int row = 0; row < n; ++row) {
    for (int column = 0; column < k; ++column) {
      const int neighbor = index(row, column);
      const double d = distance(row, column);
      if (neighbor < 1 || neighbor > n || !R_finite(d) || d < 0.0) {
        stop("BBKNN neighbors contain invalid indices or distances");
      }
      neighbors[column] = std::make_pair(d, neighbor);
    }
    std::stable_sort(
        neighbors.begin(), neighbors.end(),
        [](const std::pair<double, int>& left,
           const std::pair<double, int>& right) {
          if (left.first < right.first) return true;
          if (left.first > right.first) return false;
          return left.second < right.second;
        });

    nonzero.clear();
    for (int column = 0; column < k; ++column) {
      sorted_distance(row, column) = neighbors[column].first;
      sorted_index(row, column) = neighbors[column].second;
      if (neighbors[column].first > 0.0) {
        nonzero.push_back(neighbors[column].first);
      }
    }

    double rho = 0.0;
    if (!nonzero.empty()) {
      const int integer_connectivity =
          static_cast<int>(std::floor(local_connectivity));
      const double interpolation =
          local_connectivity - integer_connectivity;
      if (integer_connectivity > 0) {
        const int rho_index =
            std::min(integer_connectivity - 1,
                     static_cast<int>(nonzero.size()) - 1);
        rho = nonzero[rho_index];
        if (interpolation > 1e-5 &&
            rho_index + 1 < static_cast<int>(nonzero.size())) {
          rho += interpolation *
                 (nonzero[rho_index + 1] - nonzero[rho_index]);
        }
      } else {
        rho = interpolation * nonzero[0];
      }
    }

    double lower = 0.0;
    double upper = R_PosInf;
    double sigma = 1.0;
    for (int iteration = 0; iteration < 64; ++iteration) {
      double probability_sum = 0.0;
      for (int column = 1; column < k; ++column) {
        const double shifted = sorted_distance(row, column) - rho;
        probability_sum += shifted > 0.0 ? std::exp(-shifted / sigma) : 1.0;
      }
      if (std::abs(probability_sum - target) < 1e-5) {
        break;
      }
      if (probability_sum > target) {
        upper = sigma;
        sigma = (lower + upper) / 2.0;
      } else {
        lower = sigma;
        sigma = R_finite(upper) ? (lower + upper) / 2.0 : sigma * 2.0;
      }
    }
    double mean_distance = 0.0;
    for (int column = 0; column < k; ++column) {
      mean_distance += sorted_distance(row, column);
    }
    mean_distance /= static_cast<double>(k);
    sigma = std::max(sigma, 1e-3 * mean_distance);

    for (int column = 0; column < k; ++column) {
      if (sorted_index(row, column) == row + 1) {
        membership(row, column) = 0.0;
      } else {
        const double shifted = sorted_distance(row, column) - rho;
        membership(row, column) =
            shifted <= 0.0 ? 1.0 : std::exp(-shifted / sigma);
      }
    }
  }
  return List::create(
      _["idx"] = sorted_index,
      _["distance"] = sorted_distance,
      _["membership"] = membership
  );
}
