// [[Rcpp::depends(RcppArmadillo, cli)]]
#include <RcppArmadillo.h>
#include <thisutils/cli_progress.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <random>
#include <vector>

using namespace Rcpp;

namespace {

struct Edge {
  int to;
  int rev;
  int cap;
  long long cost;
};

void add_edge(std::vector<std::vector<Edge> >& graph, int from, int to, int cap, long long cost) {
  Edge a = {to, static_cast<int>(graph[to].size()), cap, cost};
  Edge b = {from, static_cast<int>(graph[from].size()), 0, -cost};
  graph[from].push_back(a);
  graph[to].push_back(b);
}

long long min_cost_flow(
  std::vector<std::vector<Edge> >& graph,
  int source,
  int sink,
  int flow,
  bool verbose
) {
  const long long INF = std::numeric_limits<long long>::max() / 4;
  const int n = graph.size();
  std::vector<long long> potential(n, 0), dist(n);
  std::vector<int> prev_node(n), prev_edge(n);
  long long total_cost = 0;
  thisutils::cli_progress progress(
    flow,
    verbose,
    "Assign sampled reference cells to spatial spots"
  );

  for (int pushed = 0; pushed < flow; ++pushed) {
    if (thisutils::should_check_interrupt(pushed, flow, verbose)) {
      Rcpp::checkUserInterrupt();
    }
    progress.set(pushed);
    std::fill(dist.begin(), dist.end(), INF);
    dist[source] = 0;
    typedef std::pair<long long, int> QueueNode;
    std::priority_queue<QueueNode, std::vector<QueueNode>, std::greater<QueueNode> > pq;
    pq.push(QueueNode(0, source));

    while (!pq.empty()) {
      QueueNode current = pq.top();
      pq.pop();
      long long d = current.first;
      int v = current.second;
      if (d != dist[v]) {
        continue;
      }
      for (int i = 0; i < static_cast<int>(graph[v].size()); ++i) {
        const Edge& e = graph[v][i];
        if (e.cap <= 0) {
          continue;
        }
        long long nd = d + e.cost + potential[v] - potential[e.to];
        if (nd < dist[e.to]) {
          dist[e.to] = nd;
          prev_node[e.to] = v;
          prev_edge[e.to] = i;
          pq.push(QueueNode(nd, e.to));
        }
      }
    }

    if (dist[sink] == INF) {
      stop("Unable to find a complete CytoSPACE assignment for the current target counts.");
    }

    for (int v = 0; v < n; ++v) {
      if (dist[v] < INF) {
        potential[v] += dist[v];
      }
    }

    int aug = 1;
    for (int v = sink; v != source; v = prev_node[v]) {
      aug = std::min(aug, graph[prev_node[v]][prev_edge[v]].cap);
    }
    for (int v = sink; v != source; v = prev_node[v]) {
      Edge& e = graph[prev_node[v]][prev_edge[v]];
      e.cap -= aug;
      graph[v][e.rev].cap += aug;
      total_cost += static_cast<long long>(aug) * e.cost;
    }
  }
  progress.set(flow, true);

  return total_cost;
}

arma::mat normalize_log_cpm(const NumericMatrix& input) {
  arma::mat out(input.nrow(), input.ncol());
  for (int j = 0; j < input.ncol(); ++j) {
    for (int i = 0; i < input.nrow(); ++i) {
      out(i, j) = input(i, j);
    }
  }
  for (arma::uword j = 0; j < out.n_cols; ++j) {
    double total = arma::accu(out.col(j));
    if (!std::isfinite(total) || total <= 0.0) {
      out.col(j).zeros();
      continue;
    }
    out.col(j) = arma::log2((out.col(j) * (1000000.0 / total)) + 1.0);
    out.col(j).replace(arma::datum::nan, 0.0);
    out.col(j).replace(arma::datum::inf, 0.0);
  }
  return out;
}

void standardize_columns_inplace(arma::mat& x) {
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double mean = arma::mean(x.col(j));
    x.col(j) -= mean;
    double ss = arma::dot(x.col(j), x.col(j));
    if (!std::isfinite(ss) || ss <= 0.0) {
      x.col(j).zeros();
    } else {
      x.col(j) /= std::sqrt(ss);
    }
  }
}

arma::mat pearson_cor_matrix(arma::mat x, arma::mat y) {
  standardize_columns_inplace(x);
  standardize_columns_inplace(y);
  arma::mat out = x.t() * y;
  out.transform([](double value) {
    if (!std::isfinite(value)) {
      return 0.0;
    }
    if (value > 1.0) {
      return 1.0;
    }
    if (value < -1.0) {
      return -1.0;
    }
    return value;
  });
  return out;
}

long long scaled_cost(double value, int distance_code, double jitter = 0.0) {
  double cost = distance_code == 1 ? value : 1.0 - value;
  if (!std::isfinite(cost) || cost < 0.0) {
    cost = 0.0;
  }
  return static_cast<long long>(std::floor(cost * 1000000.0 + jitter + 1.0));
}

double numpy_random_sample(std::mt19937& rng) {
  unsigned long a = static_cast<unsigned long>(rng()) >> 5;
  unsigned long b = static_cast<unsigned long>(rng()) >> 6;
  return (a * 67108864.0 + b) / 9007199254740992.0;
}

} // namespace

// [[Rcpp::export]]
List cytospace_assign(
  NumericMatrix sc_expr,
  NumericMatrix st_expr,
  IntegerVector spot_capacities,
  int seed = 1,
  bool upstream_tie_break = true,
  bool verbose = false
) {
  if (sc_expr.nrow() != st_expr.nrow()) {
    stop("sc_expr and st_expr must have the same number of genes.");
  }
  if (spot_capacities.size() != st_expr.ncol()) {
    stop("spot_capacities must have one value per spatial spot.");
  }

  arma::mat sc_norm = normalize_log_cpm(sc_expr);
  arma::mat st_norm = normalize_log_cpm(st_expr);
  arma::mat cor_mat = pearson_cor_matrix(sc_norm, st_norm);

  const int n_cells = sc_expr.ncol();
  const int n_spots = st_expr.ncol();
  int total_capacity = 0;
  for (int spot = 0; spot < n_spots; ++spot) {
    if (spot_capacities[spot] < 0) {
      stop("spot_capacities cannot contain negative values.");
    }
    total_capacity += spot_capacities[spot];
  }
  if (total_capacity != n_cells) {
    stop("spot_capacities must sum to the number of sampled reference cells.");
  }

  std::vector<int> out_cell(n_cells);
  std::vector<int> out_spot(n_cells);
  std::vector<double> out_score;
  long long total_cost = 0;

  const int source = 0;
  const int cell_offset = 1;
  const int spot_offset = cell_offset + n_cells;
  const int sink = spot_offset + n_spots;
  std::vector<std::vector<Edge> > graph(sink + 1);

  for (int cell = 0; cell < n_cells; ++cell) {
    add_edge(graph, source, cell_offset + cell, 1, 0);
  }
  for (int spot = 0; spot < n_spots; ++spot) {
    add_edge(graph, spot_offset + spot, sink, spot_capacities[spot], 0);
  }

  std::vector<unsigned char> jitter;
  if (upstream_tie_break) {
    jitter.resize(static_cast<std::size_t>(n_spots) * static_cast<std::size_t>(n_cells));
    std::mt19937 rng(static_cast<unsigned int>(seed));
    for (int spot = 0; spot < n_spots; ++spot) {
      for (int cell = 0; cell < n_cells; ++cell) {
        jitter[static_cast<std::size_t>(spot) * n_cells + cell] =
          static_cast<unsigned char>(std::floor(numpy_random_sample(rng) * 10.0));
      }
    }
  }

  for (int cell = 0; cell < n_cells; ++cell) {
    for (int spot = 0; spot < n_spots; ++spot) {
      double metric = cor_mat(cell, spot);
      double edge_jitter = upstream_tie_break ?
        static_cast<double>(jitter[static_cast<std::size_t>(spot) * n_cells + cell]) :
        0.0;
      add_edge(
        graph,
        cell_offset + cell,
        spot_offset + spot,
        1,
        scaled_cost(metric, 0, edge_jitter)
      );
    }
  }

  total_cost += min_cost_flow(graph, source, sink, n_cells, verbose);

  for (int cell = 0; cell < n_cells; ++cell) {
    int assigned_spot = -1;
    double score = NA_REAL;
    int node = cell_offset + cell;
    for (std::size_t e_idx = 0; e_idx < graph[node].size(); ++e_idx) {
      const Edge& e = graph[node][e_idx];
      if (e.to >= spot_offset && e.to < spot_offset + n_spots && e.cap == 0) {
        assigned_spot = e.to - spot_offset;
        score = cor_mat(cell, assigned_spot);
        break;
      }
    }
    if (assigned_spot < 0) {
      stop("Internal CytoSPACE assignment recovery failed.");
    }
    out_cell[cell] = cell + 1;
    out_spot[cell] = assigned_spot + 1;
    out_score.push_back(score);
  }

  IntegerVector cell_index(out_cell.begin(), out_cell.end());
  IntegerVector spot_index(out_spot.begin(), out_spot.end());
  NumericVector score(out_score.begin(), out_score.end());

  return List::create(
    _["cell_index"] = cell_index,
    _["spot_index"] = spot_index,
    _["score"] = score,
    _["total_cost"] = static_cast<double>(total_cost)
  );
}
