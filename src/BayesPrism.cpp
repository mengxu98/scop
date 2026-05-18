// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

using namespace Rcpp;

namespace {

struct BayesPrismSampleSpec {
  std::vector<int> genes;
  std::vector<int> counts;
};

int bayesprism_worker_count(int requested, int n_tasks) {
  if (requested <= 1 || n_tasks <= 1) {
    return 1;
  }

  int n_threads = std::min(requested, n_tasks);
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    n_threads = std::min(n_threads, static_cast<int>(hardware));
  }
  return std::max(1, n_threads);
}

BayesPrismSampleSpec bayesprism_prepare_sample(
  const double* mixture_ptr,
  int n_samples,
  int n_genes,
  int sample_idx
) {
  BayesPrismSampleSpec spec;
  spec.genes.reserve(n_genes);
  spec.counts.reserve(n_genes);

  for (int g = 0; g < n_genes; ++g) {
    const double value = mixture_ptr[sample_idx + n_samples * g];
    const int count = static_cast<int>(std::llround(value));
    if (count > 0) {
      spec.genes.push_back(g);
      spec.counts.push_back(count);
    }
  }

  return spec;
}

void bayesprism_sample_multinomial(
  int total,
  const std::vector<double>& probs,
  std::vector<int>& out,
  std::mt19937_64& rng
) {
  const int k_total = static_cast<int>(probs.size());
  std::fill(out.begin(), out.end(), 0);
  if (total <= 0 || k_total == 0) {
    return;
  }
  if (k_total == 1) {
    out[0] = total;
    return;
  }

  int remaining = total;
  double remaining_prob = 1.0;
  for (int k = 0; k < k_total - 1; ++k) {
    if (remaining <= 0) {
      break;
    }
    double prob = 0.0;
    if (remaining_prob > 0.0) {
      prob = probs[k] / remaining_prob;
    }
    prob = std::max(0.0, std::min(1.0, prob));
    std::binomial_distribution<int> binom(remaining, prob);
    const int draw = binom(rng);
    out[k] = draw;
    remaining -= draw;
    remaining_prob -= probs[k];
  }
  out[k_total - 1] = remaining;
}

void bayesprism_sample_dirichlet(
  const std::vector<double>& alpha,
  std::vector<double>& out,
  std::mt19937_64& rng
) {
  const int k_total = static_cast<int>(alpha.size());
  double total = 0.0;
  for (int k = 0; k < k_total; ++k) {
    const double shape = std::max(alpha[k], 1e-12);
    std::gamma_distribution<double> gamma(shape, 1.0);
    out[k] = gamma(rng);
    total += out[k];
  }

  if (!std::isfinite(total) || total <= 0.0) {
    const double uniform = 1.0 / static_cast<double>(k_total);
    for (int k = 0; k < k_total; ++k) {
      out[k] = uniform;
    }
    return;
  }

  for (int k = 0; k < k_total; ++k) {
    out[k] /= total;
  }
}

} // namespace

// [[Rcpp::export]]
List bayesprism_gibbs_initial_cpp(
  NumericMatrix mixture,
  NumericMatrix phi,
  IntegerVector gibbs_idx,
  double alpha = 1.0,
  int seed = 123,
  int n_threads = 1
) {
  const int n_samples = mixture.nrow();
  const int n_genes = mixture.ncol();
  const int n_states = phi.nrow();
  if (phi.ncol() != n_genes) {
    stop("phi columns must match mixture columns.");
  }
  if (gibbs_idx.size() == 0) {
    stop("gibbs_idx must contain at least one iteration.");
  }

  int max_iter = 0;
  for (int i = 0; i < gibbs_idx.size(); ++i) {
    if (gibbs_idx[i] > max_iter) {
      max_iter = gibbs_idx[i];
    }
  }
  if (max_iter < 1) {
    stop("gibbs_idx must contain positive iteration indices.");
  }

  std::vector<unsigned char> keep_flag(static_cast<std::size_t>(max_iter + 1), 0);
  for (int i = 0; i < gibbs_idx.size(); ++i) {
    const int idx = gibbs_idx[i];
    if (idx >= 1 && idx <= max_iter) {
      keep_flag[static_cast<std::size_t>(idx)] = 1;
    }
  }
  const int n_kept = std::count(keep_flag.begin(), keep_flag.end(), static_cast<unsigned char>(1));
  if (n_kept == 0) {
    stop("No retained Gibbs iterations after applying gibbs_idx.");
  }

  const double* mixture_ptr = REAL(mixture);
  const double* phi_ptr = REAL(phi);

  std::vector<BayesPrismSampleSpec> sample_specs(static_cast<std::size_t>(n_samples));
  for (int n = 0; n < n_samples; ++n) {
    sample_specs[static_cast<std::size_t>(n)] =
      bayesprism_prepare_sample(mixture_ptr, n_samples, n_genes, n);
  }

  NumericVector z_out(static_cast<R_xlen_t>(n_samples) * n_genes * n_states);
  NumericMatrix theta_out(n_samples, n_states);
  NumericMatrix theta_cv_out(n_samples, n_states);
  double* z_ptr = REAL(z_out);
  double* theta_ptr = REAL(theta_out);
  double* theta_cv_ptr = REAL(theta_cv_out);

  const int workers = bayesprism_worker_count(n_threads, n_samples);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(workers));

  auto worker_fun = [&](int worker_id) {
    for (int n = worker_id; n < n_samples; n += workers) {
      const BayesPrismSampleSpec& spec = sample_specs[static_cast<std::size_t>(n)];
      std::mt19937_64 rng(
        static_cast<std::uint64_t>(seed) + 104729ULL * static_cast<std::uint64_t>(n + 1)
      );
      std::vector<double> theta(static_cast<std::size_t>(n_states), 1.0 / static_cast<double>(n_states));
      std::vector<double> theta_sum(static_cast<std::size_t>(n_states), 0.0);
      std::vector<double> theta_sq_sum(static_cast<std::size_t>(n_states), 0.0);
      std::vector<double> dirichlet_alpha(static_cast<std::size_t>(n_states), alpha);
      std::vector<double> probs(static_cast<std::size_t>(n_states), 0.0);
      std::vector<int> z_gene(static_cast<std::size_t>(n_states), 0);
      std::vector<double> z_state_sum(static_cast<std::size_t>(n_states), 0.0);
      std::vector<double> z_sum(static_cast<std::size_t>(n_genes) * n_states, 0.0);

      for (int iter = 1; iter <= max_iter; ++iter) {
        std::fill(z_state_sum.begin(), z_state_sum.end(), 0.0);
        for (std::size_t idx = 0; idx < spec.genes.size(); ++idx) {
          const int g = spec.genes[idx];
          const int total = spec.counts[idx];
          double prob_total = 0.0;
          for (int k = 0; k < n_states; ++k) {
            const double value = phi_ptr[k + n_states * g] * theta[static_cast<std::size_t>(k)];
            probs[static_cast<std::size_t>(k)] = value;
            prob_total += value;
          }
          if (!std::isfinite(prob_total) || prob_total <= 0.0) {
            const double uniform = 1.0 / static_cast<double>(n_states);
            std::fill(probs.begin(), probs.end(), uniform);
          } else {
            for (int k = 0; k < n_states; ++k) {
              probs[static_cast<std::size_t>(k)] /= prob_total;
            }
          }

          bayesprism_sample_multinomial(total, probs, z_gene, rng);
          for (int k = 0; k < n_states; ++k) {
            const double draw = static_cast<double>(z_gene[static_cast<std::size_t>(k)]);
            z_state_sum[static_cast<std::size_t>(k)] += draw;
            if (keep_flag[static_cast<std::size_t>(iter)]) {
              z_sum[static_cast<std::size_t>(g) + static_cast<std::size_t>(n_genes) * k] += draw;
            }
          }
        }

        for (int k = 0; k < n_states; ++k) {
          dirichlet_alpha[static_cast<std::size_t>(k)] =
            z_state_sum[static_cast<std::size_t>(k)] + alpha;
        }
        bayesprism_sample_dirichlet(dirichlet_alpha, theta, rng);

        if (keep_flag[static_cast<std::size_t>(iter)]) {
          for (int k = 0; k < n_states; ++k) {
            const double value = theta[static_cast<std::size_t>(k)];
            theta_sum[static_cast<std::size_t>(k)] += value;
            theta_sq_sum[static_cast<std::size_t>(k)] += value * value;
          }
        }
      }

      for (int k = 0; k < n_states; ++k) {
        const double mean = theta_sum[static_cast<std::size_t>(k)] / static_cast<double>(n_kept);
        const double second = theta_sq_sum[static_cast<std::size_t>(k)] / static_cast<double>(n_kept);
        const double variance = std::max(0.0, second - mean * mean);
        theta_ptr[n + n_samples * k] = mean;
        theta_cv_ptr[n + n_samples * k] =
          (mean > 0.0) ? std::sqrt(variance) / mean : NA_REAL;
      }

      for (int g = 0; g < n_genes; ++g) {
        for (int k = 0; k < n_states; ++k) {
          z_ptr[n + n_samples * g + static_cast<R_xlen_t>(n_samples) * n_genes * k] =
            z_sum[static_cast<std::size_t>(g) + static_cast<std::size_t>(n_genes) * k] /
            static_cast<double>(n_kept);
        }
      }
    }
  };

  for (int worker = 0; worker < workers; ++worker) {
    pool.emplace_back(worker_fun, worker);
  }
  for (std::size_t i = 0; i < pool.size(); ++i) {
    pool[i].join();
  }

  z_out.attr("dim") = IntegerVector::create(n_samples, n_genes, n_states);
  theta_out.attr("dimnames") = List::create(R_NilValue, R_NilValue);
  theta_cv_out.attr("dimnames") = List::create(R_NilValue, R_NilValue);

  return List::create(
    _["Z"] = z_out,
    _["theta"] = theta_out,
    _["theta_cv"] = theta_cv_out,
    _["constant"] = 0.0
  );
}

// [[Rcpp::export]]
List bayesprism_gibbs_final_cpp(
  NumericMatrix mixture,
  NumericMatrix phi,
  IntegerVector gibbs_idx,
  double alpha = 1.0,
  int seed = 123,
  int n_threads = 1
) {
  const int n_samples = mixture.nrow();
  const int n_genes = mixture.ncol();
  const int n_types = phi.nrow();
  if (phi.ncol() != n_genes) {
    stop("phi columns must match mixture columns.");
  }
  if (gibbs_idx.size() == 0) {
    stop("gibbs_idx must contain at least one iteration.");
  }

  int max_iter = 0;
  for (int i = 0; i < gibbs_idx.size(); ++i) {
    if (gibbs_idx[i] > max_iter) {
      max_iter = gibbs_idx[i];
    }
  }
  if (max_iter < 1) {
    stop("gibbs_idx must contain positive iteration indices.");
  }

  std::vector<unsigned char> keep_flag(static_cast<std::size_t>(max_iter + 1), 0);
  for (int i = 0; i < gibbs_idx.size(); ++i) {
    const int idx = gibbs_idx[i];
    if (idx >= 1 && idx <= max_iter) {
      keep_flag[static_cast<std::size_t>(idx)] = 1;
    }
  }
  const int n_kept = std::count(keep_flag.begin(), keep_flag.end(), static_cast<unsigned char>(1));
  if (n_kept == 0) {
    stop("No retained Gibbs iterations after applying gibbs_idx.");
  }

  const double* mixture_ptr = REAL(mixture);
  const double* phi_ptr = REAL(phi);

  std::vector<BayesPrismSampleSpec> sample_specs(static_cast<std::size_t>(n_samples));
  for (int n = 0; n < n_samples; ++n) {
    sample_specs[static_cast<std::size_t>(n)] =
      bayesprism_prepare_sample(mixture_ptr, n_samples, n_genes, n);
  }

  NumericMatrix theta_out(n_samples, n_types);
  NumericMatrix theta_cv_out(n_samples, n_types);
  double* theta_ptr = REAL(theta_out);
  double* theta_cv_ptr = REAL(theta_cv_out);

  const int workers = bayesprism_worker_count(n_threads, n_samples);
  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(workers));

  auto worker_fun = [&](int worker_id) {
    for (int n = worker_id; n < n_samples; n += workers) {
      const BayesPrismSampleSpec& spec = sample_specs[static_cast<std::size_t>(n)];
      std::mt19937_64 rng(
        static_cast<std::uint64_t>(seed) + 104729ULL * static_cast<std::uint64_t>(n + 1)
      );
      std::vector<double> theta(static_cast<std::size_t>(n_types), 1.0 / static_cast<double>(n_types));
      std::vector<double> theta_sum(static_cast<std::size_t>(n_types), 0.0);
      std::vector<double> theta_sq_sum(static_cast<std::size_t>(n_types), 0.0);
      std::vector<double> dirichlet_alpha(static_cast<std::size_t>(n_types), alpha);
      std::vector<double> probs(static_cast<std::size_t>(n_types), 0.0);
      std::vector<int> z_gene(static_cast<std::size_t>(n_types), 0);
      std::vector<double> z_type_sum(static_cast<std::size_t>(n_types), 0.0);

      for (int iter = 1; iter <= max_iter; ++iter) {
        std::fill(z_type_sum.begin(), z_type_sum.end(), 0.0);
        for (std::size_t idx = 0; idx < spec.genes.size(); ++idx) {
          const int g = spec.genes[idx];
          const int total = spec.counts[idx];
          double prob_total = 0.0;
          for (int k = 0; k < n_types; ++k) {
            const double value = phi_ptr[k + n_types * g] * theta[static_cast<std::size_t>(k)];
            probs[static_cast<std::size_t>(k)] = value;
            prob_total += value;
          }
          if (!std::isfinite(prob_total) || prob_total <= 0.0) {
            const double uniform = 1.0 / static_cast<double>(n_types);
            std::fill(probs.begin(), probs.end(), uniform);
          } else {
            for (int k = 0; k < n_types; ++k) {
              probs[static_cast<std::size_t>(k)] /= prob_total;
            }
          }

          bayesprism_sample_multinomial(total, probs, z_gene, rng);
          for (int k = 0; k < n_types; ++k) {
            z_type_sum[static_cast<std::size_t>(k)] +=
              static_cast<double>(z_gene[static_cast<std::size_t>(k)]);
          }
        }

        for (int k = 0; k < n_types; ++k) {
          dirichlet_alpha[static_cast<std::size_t>(k)] =
            z_type_sum[static_cast<std::size_t>(k)] + alpha;
        }
        bayesprism_sample_dirichlet(dirichlet_alpha, theta, rng);

        if (keep_flag[static_cast<std::size_t>(iter)]) {
          for (int k = 0; k < n_types; ++k) {
            const double value = theta[static_cast<std::size_t>(k)];
            theta_sum[static_cast<std::size_t>(k)] += value;
            theta_sq_sum[static_cast<std::size_t>(k)] += value * value;
          }
        }
      }

      for (int k = 0; k < n_types; ++k) {
        const double mean = theta_sum[static_cast<std::size_t>(k)] / static_cast<double>(n_kept);
        const double second = theta_sq_sum[static_cast<std::size_t>(k)] / static_cast<double>(n_kept);
        const double variance = std::max(0.0, second - mean * mean);
        theta_ptr[n + n_samples * k] = mean;
        theta_cv_ptr[n + n_samples * k] =
          (mean > 0.0) ? std::sqrt(variance) / mean : NA_REAL;
      }
    }
  };

  for (int worker = 0; worker < workers; ++worker) {
    pool.emplace_back(worker_fun, worker);
  }
  for (std::size_t i = 0; i < pool.size(); ++i) {
    pool[i].join();
  }

  theta_out.attr("dimnames") = List::create(R_NilValue, R_NilValue);
  theta_cv_out.attr("dimnames") = List::create(R_NilValue, R_NilValue);

  return List::create(
    _["theta"] = theta_out,
    _["theta_cv"] = theta_cv_out
  );
}
