#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct CccAgg {
  std::string sender;
  std::string receiver;
  double sum = 0.0;
  double max = 0.0;
  double count = 0.0;
  int n = 0;
};

struct LianaAgg {
  std::string source;
  std::string target;
  std::string ligand_complex;
  std::string receptor_complex;
  std::string sample;
  double sum_score = 0.0;
  double min_pvalue = R_PosInf;
  bool has_pvalue = false;
  std::string classification;
  std::vector<std::string> methods;
  std::vector<std::string> liana_methods;
  std::vector<std::string> resources;
  std::unordered_set<std::string> methods_seen;
  std::unordered_set<std::string> liana_methods_seen;
  std::unordered_set<std::string> resources_seen;
  bool classification_set = false;
  bool use_sample = false;
  int n = 0;
};

std::string ccc_as_string(const CharacterVector& x, int i) {
  if (CharacterVector::is_na(x[i])) {
    return "NA";
  }
  return as<std::string>(x[i]);
}

double ccc_score_or_zero(double x) {
  if (NumericVector::is_na(x) || R_IsNaN(x)) {
    return 0.0;
  }
  return x;
}

double ccc_significant_or_zero(double x) {
  if (NumericVector::is_na(x) || R_IsNaN(x) || !R_finite(x)) {
    return 0.0;
  }
  return x;
}

}  // namespace

// [[Rcpp::export]]
DataFrame ccc_aggregate_long_cpp(
    CharacterVector sender,
    CharacterVector receiver,
    NumericVector score,
    NumericVector significant
) {
  const int n_in = sender.size();
  if (n_in == 0) {
    return DataFrame::create(
      _["sender"] = CharacterVector(),
      _["receiver"] = CharacterVector(),
      _["sum"] = NumericVector(),
      _["mean"] = NumericVector(),
      _["max"] = NumericVector(),
      _["count"] = NumericVector()
    );
  }
  if (receiver.size() != n_in || score.size() != n_in || significant.size() != n_in) {
    stop("ccc_aggregate_long_cpp: sender, receiver, score, and significant must have the same length");
  }

  std::unordered_map<std::string, int> index;
  std::vector<CccAgg> groups;
  index.reserve(static_cast<std::size_t>(n_in));
  groups.reserve(static_cast<std::size_t>(n_in));

  int n = n_in;
  Rcpp::CharacterVector sender_in = sender;
  Rcpp::CharacterVector receiver_in = receiver;

  for (int i = 0; i < n; ++i) {
    const std::string sender_i = ccc_as_string(sender_in, i);
    const std::string receiver_i = ccc_as_string(receiver_in, i);
    if (sender_i.empty() && receiver_i.empty()) continue;
    const std::string key = sender_i + "\n" + receiver_i;
    int pos;
    const auto found = index.find(key);
    if (found == index.end()) {
      pos = static_cast<int>(groups.size());
      index.emplace(key, pos);
      CccAgg agg;
      agg.sender = sender_i;
      agg.receiver = receiver_i;
      groups.push_back(agg);
    } else {
      pos = found->second;
    }

    const double score_i = ccc_score_or_zero(score[i]);
    const double significant_i = ccc_significant_or_zero(significant[i]);
    CccAgg& agg = groups[static_cast<std::size_t>(pos)];
    if (agg.n == 0 || score_i > agg.max) {
      agg.max = score_i;
    }
    agg.sum += score_i;
    agg.count += significant_i;
    ++agg.n;
  }

  std::vector<int> order(groups.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = static_cast<int>(i);
  }
  std::sort(order.begin(), order.end(), [&groups](int lhs, int rhs) {
    const CccAgg& a = groups[static_cast<std::size_t>(lhs)];
    const CccAgg& b = groups[static_cast<std::size_t>(rhs)];
    if (a.sender == b.sender) {
      return a.receiver < b.receiver;
    }
    return a.sender < b.sender;
  });

  const int n_groups = static_cast<int>(groups.size());
  CharacterVector out_sender(n_groups);
  CharacterVector out_receiver(n_groups);
  NumericVector out_sum(n_groups);
  NumericVector out_mean(n_groups);
  NumericVector out_max(n_groups);
  NumericVector out_count(n_groups);

  for (int i = 0; i < n_groups; ++i) {
    const CccAgg& agg = groups[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])];
    out_sender[i] = agg.sender;
    out_receiver[i] = agg.receiver;
    out_sum[i] = agg.sum;
    out_mean[i] = agg.n > 0 ? agg.sum / static_cast<double>(agg.n) : R_NaN;
    out_max[i] = agg.max;
    out_count[i] = agg.count;
  }

  return DataFrame::create(
    _["sender"] = out_sender,
    _["receiver"] = out_receiver,
    _["sum"] = out_sum,
    _["mean"] = out_mean,
    _["max"] = out_max,
    _["count"] = out_count
  );
}

// [[Rcpp::export]]
DataFrame ccc_aggregate_liana_table_cpp(
    CharacterVector source,
    CharacterVector target,
    CharacterVector ligand_complex,
    CharacterVector receptor_complex,
    CharacterVector sample,
    NumericVector score,
    NumericVector pvalue,
    CharacterVector classification,
    CharacterVector method,
    CharacterVector liana_method,
    CharacterVector resource,
    bool has_sample
) {
  const int n = source.size();
  if (n == 0) {
    return DataFrame::create(
      _["source"] = CharacterVector(),
      _["target"] = CharacterVector(),
      _["ligand_complex"] = CharacterVector(),
      _["receptor_complex"] = CharacterVector(),
      _["score"] = NumericVector(),
      _["pvalue"] = NumericVector(),
      _["classification"] = CharacterVector(),
      _["method"] = CharacterVector(),
      _["liana_method"] = CharacterVector(),
      _["resource"] = CharacterVector()
    );
  }
  if (target.size() != n || ligand_complex.size() != n || receptor_complex.size() != n ||
      score.size() != n || pvalue.size() != n) {
    stop("ccc_aggregate_liana_table_cpp: all vectors must have the same length");
  }

  auto make_key = [&](int i) -> std::string {
    std::string k = as<std::string>(source[i]) + std::string("\r") +
                    as<std::string>(target[i]) + std::string("\r") +
                    as<std::string>(ligand_complex[i]) + std::string("\r") +
                    as<std::string>(receptor_complex[i]);
    if (has_sample) {
      k += std::string("\r") + as<std::string>(sample[i]);
    }
    return k;
  };

  std::unordered_map<std::string, int> index;
  std::vector<LianaAgg> groups;
  index.reserve(static_cast<std::size_t>(n));
  groups.reserve(static_cast<std::size_t>(n));

  bool has_classification = classification.size() == n;
  bool has_method = method.size() == n;
  bool has_liana_method = liana_method.size() == n;
  bool has_resource = resource.size() == n;

  for (int i = 0; i < n; ++i) {
    const std::string key = make_key(i);
    int pos;
    const auto found = index.find(key);
    if (found == index.end()) {
      pos = static_cast<int>(groups.size());
      index.emplace(key, pos);
      LianaAgg agg;
      agg.source = as<std::string>(source[i]);
      agg.target = as<std::string>(target[i]);
      agg.ligand_complex = as<std::string>(ligand_complex[i]);
      agg.receptor_complex = as<std::string>(receptor_complex[i]);
      agg.n = 0;
      agg.use_sample = has_sample;
      if (has_sample) {
        agg.sample = as<std::string>(sample[i]);
      }
      groups.push_back(agg);
    } else {
      pos = found->second;
    }

    LianaAgg& agg = groups[static_cast<std::size_t>(pos)];

    if (has_classification && !agg.classification_set) {
      std::string cls = as<std::string>(classification[i]);
      if (!cls.empty() && cls != "NA") {
        agg.classification = cls;
        agg.classification_set = true;
      }
    }

    double s = score[i];
    if (NumericVector::is_na(s) || R_IsNaN(s)) {
      s = 0.0;
    }
    agg.sum_score += s;

    double p = pvalue[i];
    if (!NumericVector::is_na(p) && R_finite(p)) {
      if (!agg.has_pvalue || p < agg.min_pvalue) {
        agg.min_pvalue = p;
      }
      agg.has_pvalue = true;
    }

    if (has_method) {
      std::string m = as<std::string>(method[i]);
      if (!m.empty() && m != "NA") {
        if (agg.methods_seen.find(m) == agg.methods_seen.end()) {
          agg.methods_seen.insert(m);
          agg.methods.push_back(m);
        }
      }
    }
    if (has_liana_method) {
      std::string lm = as<std::string>(liana_method[i]);
      if (!lm.empty() && lm != "NA") {
        if (agg.liana_methods_seen.find(lm) == agg.liana_methods_seen.end()) {
          agg.liana_methods_seen.insert(lm);
          agg.liana_methods.push_back(lm);
        }
      }
    }
    if (has_resource) {
      std::string r = as<std::string>(resource[i]);
      if (!r.empty() && r != "NA") {
        if (agg.resources_seen.find(r) == agg.resources_seen.end()) {
          agg.resources_seen.insert(r);
          agg.resources.push_back(r);
        }
      }
    }
    ++agg.n;
  }

  auto join_vec = [](const std::vector<std::string>& v) -> std::string {
    if (v.empty()) return "";
    std::string result;
    for (std::size_t j = 0; j < v.size(); ++j) {
      if (j > 0) result += ";";
      result += v[j];
    }
    return result;
  };

  const int n_groups = static_cast<int>(groups.size());
  CharacterVector out_source(n_groups);
  CharacterVector out_target(n_groups);
  CharacterVector out_ligand(n_groups);
  CharacterVector out_receptor(n_groups);
  CharacterVector out_sample(n_groups);
  NumericVector out_score(n_groups);
  NumericVector out_pvalue(n_groups);
  CharacterVector out_classification(n_groups);
  CharacterVector out_method(n_groups);
  CharacterVector out_liana_method(n_groups);
  CharacterVector out_resource(n_groups);

  for (int i = 0; i < n_groups; ++i) {
    const LianaAgg& agg = groups[i];
    out_source[i] = agg.source;
    out_target[i] = agg.target;
    out_ligand[i] = agg.ligand_complex;
    out_receptor[i] = agg.receptor_complex;
    if (has_sample) out_sample[i] = agg.sample;
    out_score[i] = agg.sum_score;
    out_pvalue[i] = agg.has_pvalue ? agg.min_pvalue : 1.0;
    out_classification[i] = agg.classification_set ? agg.classification : "Unclassified";
out_method[i] = join_vec(agg.methods);
    out_liana_method[i] = join_vec(agg.liana_methods);
    out_resource[i] = join_vec(agg.resources);
  }

  if (has_sample) {
    return DataFrame::create(
      _["source"] = out_source,
      _["target"] = out_target,
      _["ligand_complex"] = out_ligand,
      _["receptor_complex"] = out_receptor,
      _["sample"] = out_sample,
      _["score"] = out_score,
      _["pvalue"] = out_pvalue,
      _["classification"] = out_classification,
      _["method"] = out_method,
      _["liana_method"] = out_liana_method,
      _["resource"] = out_resource
    );
  }
  return DataFrame::create(
    _["source"] = out_source,
    _["target"] = out_target,
    _["ligand_complex"] = out_ligand,
    _["receptor_complex"] = out_receptor,
    _["score"] = out_score,
    _["pvalue"] = out_pvalue,
    _["classification"] = out_classification,
    _["method"] = out_method,
    _["liana_method"] = out_liana_method,
    _["resource"] = out_resource
  );
}
