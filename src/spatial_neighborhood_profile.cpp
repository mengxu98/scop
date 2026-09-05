#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>

namespace {
struct Tile {
  long long x, y;
  bool operator==(const Tile& b) const { return x == b.x && y == b.y; }
};
struct TileHash {
  size_t operator()(const Tile& a) const {
    size_t h = std::hash<long long>{}(a.x);
    return h ^ (std::hash<long long>{}(a.y) + 0x9e3779b9 + (h << 6) + (h >> 2));
  }
};
}

// Validated R entry point supplies one sample, 1-based IDs and sorted radii.
// Direct reduction avoids materializing a potentially large cell-cell edge list.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector spatial_neighborhood_profile_cpp(
    Rcpp::NumericMatrix xy, Rcpp::IntegerVector group, Rcpp::IntegerVector query,
    Rcpp::NumericVector radii, int ng) {
  const int n = xy.nrow(), nq = query.size(), nr = radii.size();
  if (xy.ncol() != 2 || group.size() != n || nr < 1 || ng < 1)
    Rcpp::stop("Invalid neighborhood profile dimensions");
  for (int g : group) if (g < 1 || g > ng) Rcpp::stop("Invalid group index");
  for (int q : query) if (q < 1 || q > n) Rcpp::stop("Invalid target index");
  Rcpp::IntegerVector out(static_cast<R_xlen_t>(nq) * ng * nr);
  out.attr("dim") = Rcpp::IntegerVector::create(nq, ng, nr);
  const double width = radii[nr - 1];
  std::vector<double> r2(nr);
  for (int b = 0; b < nr; ++b) r2[b] = radii[b] * radii[b];
  std::unordered_map<Tile, std::vector<int>, TileHash> tiles;
  tiles.reserve(n);
  for (int j = 0; j < n; ++j) {
    if (j % 16384 == 0) Rcpp::checkUserInterrupt();
    tiles[{static_cast<long long>(std::floor(xy(j, 0) / width)),
           static_cast<long long>(std::floor(xy(j, 1) / width))}].push_back(j);
  }
  for (int i = 0; i < nq; ++i) {
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
    const int q = query[i] - 1;
    const auto tx = static_cast<long long>(std::floor(xy(q, 0) / width));
    const auto ty = static_cast<long long>(std::floor(xy(q, 1) / width));
    for (int dx = -1; dx <= 1; ++dx) for (int dy = -1; dy <= 1; ++dy) {
      auto tile = tiles.find({tx + dx, ty + dy});
      if (tile == tiles.end()) continue;
      for (int j : tile->second) {
        if (q == j) continue;
        const double x = xy(q, 0) - xy(j, 0), y = xy(q, 1) - xy(j, 1);
        const double d2 = x * x + y * y;
        if (d2 > r2.back()) continue;
        const int b = std::lower_bound(r2.begin(), r2.end(), d2) - r2.begin();
        out[i + static_cast<R_xlen_t>(nq) * (group[j] - 1 + ng * b)]++;
      }
    }
  }
  return out;
}
