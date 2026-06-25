

#ifndef SCOP_UMAP_H
#define SCOP_UMAP_H

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

#define UMAP_SMOOTH_K_TOLERANCE 1.0e-5f
#define UMAP_MIN_K_DIST_SCALE   1.0e-3f


SCOP_API void umap_smooth_knn_dist_f32(
    int          n_samples,
    int          n_neighbors,
    const float *distances,
    float        k,
    int          n_iter,
    float        local_connectivity,
    float        bandwidth,
    float       *sigmas,
    float       *rhos
);


SCOP_API int umap_smooth_knn_dist_parallel_f32(
    int          n_samples,
    int          n_neighbors,
    const float *distances,
    float        k,
    int          n_iter,
    float        local_connectivity,
    float        bandwidth,
    float       *sigmas,
    float       *rhos,
    int          cores
);


SCOP_API void umap_membership_strengths_f32(
    int            n_samples,
    int            n_neighbors,
    const int32_t *knn_indices,
    const float   *knn_dists,
    const float   *sigmas,
    const float   *rhos,
    int            return_dists,
    int            bipartite,
    int32_t       *rows,
    int32_t       *cols,
    float         *vals,
    float         *dists
);

SCOP_API int umap_membership_strengths_parallel_f32(
    int            n_samples,
    int            n_neighbors,
    const int32_t *knn_indices,
    const float   *knn_dists,
    const float   *sigmas,
    const float   *rhos,
    int            return_dists,
    int            bipartite,
    int32_t       *rows,
    int32_t       *cols,
    float         *vals,
    float         *dists,
    int            cores
);


SCOP_API void umap_symmetrize_fuzzy_graph_f32(
    int            n_samples,
    int            n_neighbors,
    const int32_t *rows,
    const int32_t *cols,
    const float   *vals,
    float          set_op_mix_ratio,
    int32_t       *out_rows,
    int32_t       *out_cols,
    float         *out_vals
);

SCOP_API int umap_symmetrize_fuzzy_graph_parallel_f32(
    int            n_samples,
    int            n_neighbors,
    const int32_t *rows,
    const int32_t *cols,
    const float   *vals,
    float          set_op_mix_ratio,
    int32_t       *out_rows,
    int32_t       *out_cols,
    float         *out_vals,
    int            cores
);



SCOP_API void umap_optimize_layout_euclidean_f32(
    float         *embedding,
    int            n_samples,
    int            dim,
    const int32_t *head,
    const int32_t *tail,
    int64_t        n_edges,
    const float   *epochs_per_sample,
    int            n_epochs,
    float          a,
    float          b,
    float          gamma,
    float          initial_alpha,
    float          negative_sample_rate,
    uint64_t       seed
);


SCOP_API int umap_optimize_layout_euclidean_parallel_f32(
    float         *embedding,
    int            n_samples,
    int            dim,
    const int32_t *head,
    const int32_t *tail,
    int64_t        n_edges,
    const float   *epochs_per_sample,
    int            n_epochs,
    float          a,
    float          b,
    float          gamma,
    float          initial_alpha,
    float          negative_sample_rate,
    uint64_t       seed,
    int            cores
);

#ifdef __cplusplus
}
#endif

#endif
