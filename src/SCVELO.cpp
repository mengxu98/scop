#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List scvelo_stochastic_embedding_cpp(
  NumericMatrix spliced,
  NumericMatrix unspliced,
  IntegerMatrix knn_idx,
  NumericMatrix embedding
) {
  const int n_genes = spliced.nrow();
  const int n_cells = spliced.ncol();
  const int n_neighbors = knn_idx.ncol();
  const int n_dims = embedding.ncol();

  if (unspliced.nrow() != n_genes || unspliced.ncol() != n_cells) {
    stop("spliced and unspliced must have identical dimensions");
  }
  if (knn_idx.nrow() != n_cells) {
    stop("knn_idx rows must match the number of cells");
  }
  if (embedding.nrow() != n_cells) {
    stop("embedding rows must match the number of cells");
  }

  const double* spliced_ptr = spliced.begin();
  const double* unspliced_ptr = unspliced.begin();
  const double* embedding_ptr = embedding.begin();
  const int* knn_ptr = knn_idx.begin();

  NumericMatrix ms(n_genes, n_cells);
  NumericMatrix mu(n_genes, n_cells);
  double* ms_ptr = ms.begin();
  double* mu_ptr = mu.begin();

  for (int cell = 0; cell < n_cells; ++cell) {
    int count = 1;
    const int cell_offset = cell * n_genes;
    for (int gene = 0; gene < n_genes; ++gene) {
      const int idx = cell_offset + gene;
      ms_ptr[idx] = spliced_ptr[idx];
      mu_ptr[idx] = unspliced_ptr[idx];
    }
    for (int col = 0; col < n_neighbors; ++col) {
      const int neighbor = knn_ptr[cell + col * n_cells];
      if (neighbor == NA_INTEGER) {
        continue;
      }
      const int nb = neighbor - 1;
      if (nb < 0 || nb >= n_cells || nb == cell) {
        continue;
      }
      ++count;
      const int nb_offset = nb * n_genes;
      for (int gene = 0; gene < n_genes; ++gene) {
        const int cell_idx = cell_offset + gene;
        const int nb_idx = nb_offset + gene;
        ms_ptr[cell_idx] += spliced_ptr[nb_idx];
        mu_ptr[cell_idx] += unspliced_ptr[nb_idx];
      }
    }
    const double inv_count = 1.0 / static_cast<double>(count);
    for (int gene = 0; gene < n_genes; ++gene) {
      const int idx = cell_offset + gene;
      ms_ptr[idx] *= inv_count;
      mu_ptr[idx] *= inv_count;
    }
  }

  NumericVector gamma(n_genes);
  for (int gene = 0; gene < n_genes; ++gene) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int cell = 0; cell < n_cells; ++cell) {
      const int idx = cell * n_genes + gene;
      const double s = ms_ptr[idx];
      numerator += s * mu_ptr[idx];
      denominator += s * s;
    }
    gamma[gene] = denominator > 1e-12 ? numerator / denominator : 0.0;
  }

  NumericMatrix residual(n_genes, n_cells);
  NumericVector velocity_norm(n_cells);
  double* residual_ptr = residual.begin();
  for (int cell = 0; cell < n_cells; ++cell) {
    double norm_velocity_sq = 0.0;
    const int cell_offset = cell * n_genes;
    for (int gene = 0; gene < n_genes; ++gene) {
      const int idx = cell_offset + gene;
      const double v = mu_ptr[idx] - gamma[gene] * ms_ptr[idx];
      residual_ptr[idx] = v;
      norm_velocity_sq += v * v;
    }
    velocity_norm[cell] = std::sqrt(norm_velocity_sq);
  }

  NumericMatrix velocity_embedding(n_cells, n_dims);
  NumericVector confidence(n_cells);
  NumericVector velocity_length(n_cells);
  double* velocity_embedding_ptr = velocity_embedding.begin();

  for (int cell = 0; cell < n_cells; ++cell) {
    const double norm_velocity = velocity_norm[cell];
    velocity_length[cell] = norm_velocity;
    if (norm_velocity <= 0.0) {
      continue;
    }

    double weight_sum = 0.0;
    int used = 0;
    const int cell_offset = cell * n_genes;
    for (int col = 0; col < n_neighbors; ++col) {
      const int neighbor = knn_ptr[cell + col * n_cells];
      if (neighbor == NA_INTEGER) {
        continue;
      }
      const int nb = neighbor - 1;
      if (nb < 0 || nb >= n_cells || nb == cell) {
        continue;
      }

      double dot = 0.0;
      double norm_delta_sq = 0.0;
      const int nb_offset = nb * n_genes;
      for (int gene = 0; gene < n_genes; ++gene) {
        const int cell_idx = cell_offset + gene;
        const double delta = ms_ptr[nb_offset + gene] - ms_ptr[cell_idx];
        dot += residual_ptr[cell_idx] * delta;
        norm_delta_sq += delta * delta;
      }
      const double norm_delta = std::sqrt(norm_delta_sq);
      if (norm_delta <= 0.0) {
        continue;
      }
      const double cosine = dot / (norm_velocity * norm_delta);
      if (cosine <= 0.0 || !std::isfinite(cosine)) {
        continue;
      }
      for (int dim = 0; dim < n_dims; ++dim) {
        velocity_embedding_ptr[cell + dim * n_cells] +=
          cosine * (
            embedding_ptr[nb + dim * n_cells] -
            embedding_ptr[cell + dim * n_cells]
          );
      }
      weight_sum += cosine;
      ++used;
    }

    if (weight_sum > 0.0) {
      for (int dim = 0; dim < n_dims; ++dim) {
        velocity_embedding_ptr[cell + dim * n_cells] /= weight_sum;
      }
      confidence[cell] = static_cast<double>(used) / static_cast<double>(n_neighbors);
    }
  }

  return List::create(
    _["velocity_embedding"] = velocity_embedding,
    _["confidence"] = confidence,
    _["velocity_length"] = velocity_length,
    _["gamma"] = gamma
  );
}
