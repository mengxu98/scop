#ifndef SCOP_THREAD_UTILS_H
#define SCOP_THREAD_UTILS_H

#include <algorithm>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

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

inline int omp_thread_count(int requested, int n_tasks = -1) {
#ifdef _OPENMP
  int available = omp_get_max_threads();
  if (available < 1) {
    available = 1;
  }
  int n = requested > 0 ? std::min(requested, available) : available;
#else
  int n = 1;
#endif
  if (n_tasks >= 0 && n > n_tasks) {
    n = n_tasks;
  }
  if (n < 1) {
    n = 1;
  }
  return n;
}

#endif
