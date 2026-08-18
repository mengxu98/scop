#ifndef SCOP_THREAD_UTILS_H
#define SCOP_THREAD_UTILS_H

#include <algorithm>
#include <thread>

// Shared worker-count heuristic used by the parallel C++ backends.
// Returns the number of threads to spawn for `n_tasks` tasks given a
// requested core count, capped by the number of tasks and the hardware.
inline int worker_count(int requested, int n_tasks) {
  if (requested <= 1 || n_tasks <= 1) {
    return 1;
  }
  int cores = std::min(requested, n_tasks);
  const unsigned int hardware = std::thread::hardware_concurrency();
  if (hardware > 0) {
    cores = std::min(cores, static_cast<int>(hardware));
  }
  return std::max(1, cores);
}

#endif // SCOP_THREAD_UTILS_H
