#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct Pair {
  int left;
  int right;
};

struct Triplet {
  int anchor;
  int near;
  int far;
  double weight;
};

inline double row_distance(
    const double* matrix,
    int rows,
    int cols,
    int left,
    int right,
    int metric) {
  double value = 0.0;
  double left_norm = 0.0;
  double right_norm = 0.0;
  double dot = 0.0;
  for (int col = 0; col < cols; ++col) {
    const double x = matrix[left + rows * col];
    const double y = matrix[right + rows * col];
    if (metric == 1) {
      const double difference = x - y;
      value += difference * difference;
    } else if (metric == 2) {
      value += std::abs(x - y);
    } else if (metric == 3) {
      dot += x * y;
      left_norm += x * x;
      right_norm += y * y;
    } else {
      value += x == y ? 0.0 : 1.0;
    }
  }
  if (metric == 1) {
    return std::sqrt(std::max(0.0, value));
  }
  if (metric == 3) {
    const double denominator = std::sqrt(left_norm * right_norm);
    return denominator > 0.0 ? 1.0 - dot / denominator : 1.0;
  }
  if (metric == 4) {
    return value / static_cast<double>(cols);
  }
  return value;
}

inline double embedding_distance_squared(
    const std::vector<double>& embedding,
    int rows,
    int dimensions,
    int left,
    int right) {
  double distance = 0.0;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const double difference =
      embedding[left + rows * dimension] -
      embedding[right + rows * dimension];
    distance += difference * difference;
  }
  return distance;
}

inline void center_embedding(
    std::vector<double>& embedding,
    int rows,
    int dimensions) {
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    double mean = 0.0;
    for (int row = 0; row < rows; ++row) {
      mean += embedding[row + rows * dimension];
    }
    mean /= static_cast<double>(rows);
    for (int row = 0; row < rows; ++row) {
      embedding[row + rows * dimension] -= mean;
    }
  }
}

inline void add_pair_gradient(
    const Pair& pair,
    double weight,
    double scale,
    bool repulsive,
    const std::vector<double>& embedding,
    std::vector<double>& gradient,
    int rows,
    int dimensions) {
  const double distance = embedding_distance_squared(
    embedding,
    rows,
    dimensions,
    pair.left,
    pair.right
  );
  const double denominator = scale + distance;
  double coefficient = 2.0 * weight * (scale - 1.0) /
    (denominator * denominator);
  if (repulsive) {
    coefficient = -coefficient;
  }
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const int left_index = pair.left + rows * dimension;
    const int right_index = pair.right + rows * dimension;
    const double difference =
      embedding[left_index] - embedding[right_index];
    const double contribution = coefficient * difference;
    gradient[left_index] += contribution;
    gradient[right_index] -= contribution;
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List manifold_exact_knn_cpp(
    Rcpp::NumericMatrix data,
    int k,
    int metric = 1) {
  const int rows = data.nrow();
  const int cols = data.ncol();
  if (rows < 2 || cols < 1 || k < 1 || k > rows || metric < 1 || metric > 4) {
    Rcpp::stop("Invalid native manifold KNN inputs.");
  }
  const double* matrix = data.begin();
  Rcpp::IntegerMatrix index(rows, k);
  Rcpp::NumericMatrix distance(rows, k);
  int* index_ptr = index.begin();
  double* distance_ptr = distance.begin();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector<std::pair<double, int> > candidates(rows);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int query = 0; query < rows; ++query) {
      for (int reference = 0; reference < rows; ++reference) {
        candidates[reference] = std::make_pair(
          row_distance(matrix, rows, cols, query, reference, metric),
          reference
        );
      }
      std::partial_sort(
        candidates.begin(),
        candidates.begin() + k,
        candidates.end(),
        [](const std::pair<double, int>& left,
           const std::pair<double, int>& right) {
          if (left.first < right.first) {
            return true;
          }
          if (left.first > right.first) {
            return false;
          }
          return left.second < right.second;
        }
      );
      for (int neighbor = 0; neighbor < k; ++neighbor) {
        index_ptr[query + rows * neighbor] =
          candidates[neighbor].second + 1;
        distance_ptr[query + rows * neighbor] =
          candidates[neighbor].first;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::_["idx"] = index,
    Rcpp::_["distance"] = distance
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pacmap_optimize_cpp(
    Rcpp::NumericMatrix data,
    Rcpp::NumericMatrix initial,
    Rcpp::IntegerMatrix knn_index,
    int n_mid,
    int n_far,
    double learning_rate,
    int iterations,
    int seed,
    int metric = 1) {
  const int rows = data.nrow();
  const int input_dimensions = data.ncol();
  const int dimensions = initial.ncol();
  if (
    initial.nrow() != rows ||
    knn_index.nrow() != rows ||
    rows < 2 ||
    dimensions < 1 ||
    iterations < 1
  ) {
    Rcpp::stop("Invalid PaCMAP native optimizer inputs.");
  }
  const double* input = data.begin();
  std::mt19937 generator(static_cast<unsigned int>(seed));
  std::uniform_int_distribution<int> random_row(0, rows - 1);
  std::vector<Pair> near_pairs;
  std::vector<Pair> mid_pairs;
  std::vector<Pair> far_pairs;
  std::vector<std::unordered_set<int> > neighbor_sets(rows);

  for (int row = 0; row < rows; ++row) {
    for (int column = 0; column < knn_index.ncol(); ++column) {
      const int neighbor = knn_index(row, column) - 1;
      if (neighbor >= 0 && neighbor < rows && neighbor != row) {
        if (neighbor_sets[row].insert(neighbor).second) {
          near_pairs.push_back(Pair{row, neighbor});
        }
      }
    }
    for (int pair_index = 0; pair_index < n_mid; ++pair_index) {
      std::vector<std::pair<double, int> > candidates;
      candidates.reserve(6);
      for (int attempt = 0; attempt < 6; ++attempt) {
        int candidate = random_row(generator);
        if (candidate == row) {
          candidate = (candidate + 1) % rows;
        }
        candidates.push_back(std::make_pair(
          row_distance(
            input,
            rows,
            input_dimensions,
            row,
            candidate,
            metric
          ),
          candidate
        ));
      }
      std::sort(candidates.begin(), candidates.end());
      mid_pairs.push_back(Pair{
        row,
        candidates[std::min<std::size_t>(1, candidates.size() - 1)].second
      });
    }
    for (int pair_index = 0; pair_index < n_far; ++pair_index) {
      int candidate = random_row(generator);
      int attempts = 0;
      while (
        (candidate == row || neighbor_sets[row].count(candidate) > 0) &&
        attempts < 50
      ) {
        candidate = random_row(generator);
        ++attempts;
      }
      if (candidate != row) {
        far_pairs.push_back(Pair{row, candidate});
      }
    }
  }

  std::vector<double> embedding(initial.begin(), initial.end());
  std::vector<double> gradient(embedding.size(), 0.0);
  std::vector<double> first_moment(embedding.size(), 0.0);
  std::vector<double> second_moment(embedding.size(), 0.0);
  center_embedding(embedding, rows, dimensions);
  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::fill(gradient.begin(), gradient.end(), 0.0);
    double near_weight = 1.0;
    double mid_weight = 0.0;
    double far_weight = 1.0;
    if (iteration < 100) {
      const double progress = static_cast<double>(iteration) / 100.0;
      near_weight = 2.0;
      mid_weight = 1000.0 * (1.0 - progress) + 3.0 * progress;
    } else if (iteration < 200) {
      near_weight = 3.0;
      mid_weight = 3.0;
    }
    for (std::vector<Pair>::const_iterator it = near_pairs.begin();
         it != near_pairs.end(); ++it) {
      add_pair_gradient(
        *it, near_weight, 11.0, false,
        embedding, gradient, rows, dimensions
      );
    }
    for (std::vector<Pair>::const_iterator it = mid_pairs.begin();
         it != mid_pairs.end(); ++it) {
      add_pair_gradient(
        *it, mid_weight, 10001.0, false,
        embedding, gradient, rows, dimensions
      );
    }
    for (std::vector<Pair>::const_iterator it = far_pairs.begin();
         it != far_pairs.end(); ++it) {
      add_pair_gradient(
        *it, far_weight, 2.0, true,
        embedding, gradient, rows, dimensions
      );
    }
    const double beta1 = 0.9;
    const double beta2 = 0.999;
    const double corrected_rate =
      learning_rate * std::sqrt(1.0 - std::pow(beta2, iteration + 1)) /
      (1.0 - std::pow(beta1, iteration + 1));
    for (std::size_t index = 0; index < embedding.size(); ++index) {
      first_moment[index] +=
        (1.0 - beta1) * (gradient[index] - first_moment[index]);
      second_moment[index] +=
        (1.0 - beta2) *
        (gradient[index] * gradient[index] - second_moment[index]);
      embedding[index] -= corrected_rate * first_moment[index] /
        (std::sqrt(second_moment[index]) + 1e-7);
    }
    center_embedding(embedding, rows, dimensions);
  }

  Rcpp::NumericMatrix output(rows, dimensions);
  std::copy(embedding.begin(), embedding.end(), output.begin());
  return output;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix trimap_optimize_cpp(
    Rcpp::NumericMatrix data,
    Rcpp::NumericMatrix initial,
    Rcpp::IntegerMatrix knn_index,
    int n_outliers,
    int n_random,
    double learning_rate,
    int iterations,
    int optimizer,
    int seed,
    int metric = 1) {
  const int rows = data.nrow();
  const int input_dimensions = data.ncol();
  const int dimensions = initial.ncol();
  if (
    initial.nrow() != rows ||
    knn_index.nrow() != rows ||
    rows < 3 ||
    dimensions < 1 ||
    iterations < 1 ||
    optimizer < 1 ||
    optimizer > 3
  ) {
    Rcpp::stop("Invalid TriMap native optimizer inputs.");
  }
  const double* input = data.begin();
  std::mt19937 generator(static_cast<unsigned int>(seed));
  std::uniform_int_distribution<int> random_row(0, rows - 1);
  std::vector<Triplet> triplets;
  std::vector<double> sigma(rows, 1e-10);
  for (int row = 0; row < rows; ++row) {
    const int first = std::min(3, knn_index.ncol() - 1);
    const int last = std::min(5, knn_index.ncol() - 1);
    double total = 0.0;
    int count = 0;
    for (int column = first; column <= last; ++column) {
      const int neighbor = knn_index(row, column) - 1;
      if (neighbor >= 0 && neighbor < rows && neighbor != row) {
        total += row_distance(
          input, rows, input_dimensions, row, neighbor, metric
        );
        ++count;
      }
    }
    if (count > 0) {
      sigma[row] = std::max(total / static_cast<double>(count), 1e-10);
    }
  }

  for (int anchor = 0; anchor < rows; ++anchor) {
    for (int column = 0; column < knn_index.ncol(); ++column) {
      const int near = knn_index(anchor, column) - 1;
      if (near < 0 || near >= rows || near == anchor) {
        continue;
      }
      const double near_distance = row_distance(
        input, rows, input_dimensions, anchor, near, metric
      );
      for (int outlier_index = 0; outlier_index < n_outliers;
           ++outlier_index) {
        int far = random_row(generator);
        while (far == anchor || far == near) {
          far = random_row(generator);
        }
        const double far_distance = row_distance(
          input, rows, input_dimensions, anchor, far, metric
        );
        const double near_similarity =
          -(near_distance * near_distance) /
          (sigma[anchor] * sigma[near]);
        const double far_similarity =
          -(far_distance * far_distance) /
          (sigma[anchor] * sigma[far]);
        triplets.push_back(Triplet{
          anchor,
          near,
          far,
          near_similarity - far_similarity
        });
      }
    }
    for (int random_index = 0; random_index < n_random; ++random_index) {
      int first = random_row(generator);
      int second = random_row(generator);
      while (first == anchor) {
        first = random_row(generator);
      }
      while (second == anchor || second == first) {
        second = random_row(generator);
      }
      const double first_distance = row_distance(
        input, rows, input_dimensions, anchor, first, metric
      );
      const double second_distance = row_distance(
        input, rows, input_dimensions, anchor, second, metric
      );
      const int near = first_distance <= second_distance ? first : second;
      const int far = first_distance <= second_distance ? second : first;
      const double near_similarity =
        -(std::min(first_distance, second_distance) *
          std::min(first_distance, second_distance)) /
        (sigma[anchor] * sigma[near]);
      const double far_similarity =
        -(std::max(first_distance, second_distance) *
          std::max(first_distance, second_distance)) /
        (sigma[anchor] * sigma[far]);
      triplets.push_back(Triplet{
        anchor, near, far, 0.1 * (near_similarity - far_similarity)
      });
    }
  }

  double minimum_weight = R_PosInf;
  if (triplets.empty()) {
    Rcpp::stop("TriMap requires at least one sampled triplet.");
  }
  for (std::vector<Triplet>::const_iterator it = triplets.begin();
       it != triplets.end(); ++it) {
    if (std::isfinite(it->weight)) {
      minimum_weight = std::min(minimum_weight, it->weight);
    }
  }
  for (std::vector<Triplet>::iterator it = triplets.begin();
       it != triplets.end(); ++it) {
    const double shifted = std::max(0.0, it->weight - minimum_weight);
    it->weight = 2.0 * (std::sqrt(1.0 + shifted) - 1.0);
  }

  std::vector<double> embedding(initial.begin(), initial.end());
  std::vector<double> gradient(embedding.size(), 0.0);
  std::vector<double> velocity(embedding.size(), 0.0);
  std::vector<double> gains(embedding.size(), 1.0);
  center_embedding(embedding, rows, dimensions);
  double current_learning_rate = learning_rate;
  double previous_loss = R_PosInf;

  for (int iteration = 0; iteration < iterations; ++iteration) {
    std::fill(gradient.begin(), gradient.end(), 0.0);
    const double momentum = iteration > 250 ? 0.8 : 0.5;
    std::vector<double> lookahead;
    const std::vector<double>* evaluated = &embedding;
    if (optimizer != 2) {
      lookahead.resize(embedding.size());
      for (std::size_t index = 0; index < embedding.size(); ++index) {
        lookahead[index] = embedding[index] + momentum * velocity[index];
      }
      evaluated = &lookahead;
    }
    double loss = 0.0;
    for (std::vector<Triplet>::const_iterator it = triplets.begin();
         it != triplets.end(); ++it) {
      const double near_distance = 1.0 + embedding_distance_squared(
        *evaluated, rows, dimensions, it->anchor, it->near
      );
      const double far_distance = 1.0 + embedding_distance_squared(
        *evaluated, rows, dimensions, it->anchor, it->far
      );
      const double denominator = near_distance + far_distance;
      loss += it->weight * near_distance / denominator;
      const double near_coefficient =
        it->weight * far_distance / (denominator * denominator);
      const double far_coefficient =
        -it->weight * near_distance / (denominator * denominator);
      for (int dimension = 0; dimension < dimensions; ++dimension) {
        const int anchor_index = it->anchor + rows * dimension;
        const int near_index = it->near + rows * dimension;
        const int far_index = it->far + rows * dimension;
        const double near_difference =
          (*evaluated)[anchor_index] - (*evaluated)[near_index];
        const double far_difference =
          (*evaluated)[anchor_index] - (*evaluated)[far_index];
        gradient[anchor_index] +=
          near_coefficient * near_difference +
          far_coefficient * far_difference;
        gradient[near_index] -= near_coefficient * near_difference;
        gradient[far_index] -= far_coefficient * far_difference;
      }
    }
    for (std::size_t index = 0; index < embedding.size(); ++index) {
      const double current_gradient = gradient[index];
      if (optimizer == 1) {
        const int velocity_sign =
          velocity[index] > 0.0 ? 1 : (velocity[index] < 0.0 ? -1 : 0);
        const int gradient_sign =
          current_gradient > 0.0 ? 1 : (current_gradient < 0.0 ? -1 : 0);
        const bool changed_sign =
          velocity_sign != gradient_sign;
        gains[index] = changed_sign ?
          gains[index] + 0.2 : std::max(0.01, gains[index] * 0.8);
        velocity[index] =
          momentum * velocity[index] -
          current_learning_rate * gains[index] * current_gradient;
        embedding[index] += velocity[index];
      } else if (optimizer == 2) {
        embedding[index] -= current_learning_rate * current_gradient;
      } else {
        velocity[index] =
          momentum * velocity[index] -
          current_learning_rate * current_gradient;
        embedding[index] += velocity[index];
      }
    }
    const double normalized_loss =
      loss / static_cast<double>(triplets.size());
    if (optimizer != 1) {
      current_learning_rate *=
        previous_loss > normalized_loss + 1e-7 ? 1.01 : 0.9;
    }
    previous_loss = normalized_loss;
  }

  Rcpp::NumericMatrix output(rows, dimensions);
  std::copy(embedding.begin(), embedding.end(), output.begin());
  return output;
}
