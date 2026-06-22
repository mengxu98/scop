

#ifndef SCOP_KNN_H
#define SCOP_KNN_H

#include <stdint.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#define SCOP_API __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
#define SCOP_API __attribute__((visibility("default")))
#else
#define SCOP_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

SCOP_API int knn_descent_f32(
    int            n_cells,
    int            n_dims,
    int            k,
    const float   *X,
    int            n_trees,
    int            n_iters,
    int            leaf_size,
    uint64_t       seed,
    int32_t       *knn_idx,
    float         *knn_dist,
    int            cores);

#ifdef __cplusplus
}
#endif

#endif
