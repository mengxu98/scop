#ifndef SCOP_LOG_MESSAGE_H
#define SCOP_LOG_MESSAGE_H

#include <Rcpp.h>
#include <cli/progress.h>
#include <algorithm>
#include <ctime>
#include <sstream>
#include <string>

namespace scop {

inline std::string current_timestamp() {
  std::time_t now = std::time(nullptr);
  std::tm tm_buf;
#if defined(_WIN32)
  localtime_s(&tm_buf, &now);
#else
  localtime_r(&now, &tm_buf);
#endif
  char buffer[20];
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &tm_buf);
  return std::string(buffer);
}

class CliProgress {
 public:
  CliProgress(int total, bool enabled, const std::string& message) :
    bar_(R_NilValue),
    enabled_(enabled && total > 0) {
    if (enabled_) {
      std::ostringstream label;
      label << "[" << current_timestamp() << "] " << message;
      Rcpp::List config = Rcpp::List::create(
        Rcpp::_["name"] = label.str(),
        Rcpp::_["type"] = "iterator",
        Rcpp::_["format"] = "{cli::symbol$info} {cli::pb_name} {cli::pb_bar} {cli::pb_percent} ETA: {cli::pb_eta}",
        Rcpp::_["clear"] = true
      );
      bar_ = PROTECT(cli_progress_bar(static_cast<double>(total), config));
    }
  }

  ~CliProgress() {
    if (enabled_) {
      cli_progress_done(bar_);
      UNPROTECT(1);
    }
  }

  void set(int value, bool force = false) {
    if (enabled_ && (force || CLI_SHOULD_TICK)) {
      cli_progress_update(
        bar_,
        static_cast<double>(value),
        0,
        force ? 1 : 0
      );
    }
  }

 private:
  SEXP bar_;
  bool enabled_;
};

inline bool should_check_interrupt(int current, int total, bool progress_enabled) {
  const int step = std::max(1, total / 100);
  return current % step == 0 || (progress_enabled && CLI_SHOULD_TICK);
}

} // namespace scop

#endif
