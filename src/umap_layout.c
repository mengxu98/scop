#include "umap.h"
#include "umap_layout.h"
#include <stdlib.h>

typedef struct {
  double *positive;
  double *negative;
  double *negative_step;
} EdgeSchedule;

static void release_schedule(EdgeSchedule schedule) {
  free(schedule.positive);
  free(schedule.negative);
  free(schedule.negative_step);
}

static int make_schedule(EdgeSchedule* schedule, const float* epochs,
                         int64_t edges, float negative_sample_rate) {
  schedule->positive = (double*)malloc((size_t)edges * sizeof(double));
  schedule->negative = (double*)malloc((size_t)edges * sizeof(double));
  schedule->negative_step = (double*)malloc((size_t)edges * sizeof(double));
  if (!schedule->positive || !schedule->negative || !schedule->negative_step) {
    release_schedule(*schedule);
    return 0;
  }
  for (int64_t edge = 0; edge < edges; ++edge) {
    double every = (double)epochs[edge];
    schedule->negative_step[edge] = every / (double)negative_sample_rate;
    schedule->positive[edge] = every;
    schedule->negative[edge] = schedule->negative_step[edge];
  }
  return 1;
}

static void apply_positive(float* left, float* right, int dims,
                           float a, float b, float alpha) {
  float scale = layout_pull(layout_distance2(left, right, dims), a, b);
  for (int dim = 0; dim < dims; ++dim) {
    float delta = layout_bound(scale * (left[dim] - right[dim])) * alpha;
    left[dim] += delta;
    right[dim] -= delta;
  }
}

static void apply_negative(float* current, float* other, int dims,
                           float a, float b, float gamma, float alpha) {
  float distance2 = layout_distance2(current, other, dims);
  float scale = layout_push(distance2, a, b, gamma);
  if (scale <= 0.0f) return;
  for (int dim = 0; dim < dims; ++dim) {
    current[dim] += layout_bound(scale * (current[dim] - other[dim])) * alpha;
  }
}

static void update_edge(float* embedding, int vertices, int dims,
                        const int32_t* head, const int32_t* tail,
                        int64_t edge, int epoch, EdgeSchedule* schedule,
                        const float* epochs, float a, float b, float gamma,
                        float alpha, uint64_t seed) {
  const int32_t source = head[edge];
  const int32_t target = tail[edge];
  float* current = embedding + (int64_t)source * dims;
  float* paired = embedding + (int64_t)target * dims;
  apply_positive(current, paired, dims, a, b, alpha);
  schedule->positive[edge] += (double)epochs[edge];

  int draws = (int)(((double)epoch - schedule->negative[edge]) /
                    schedule->negative_step[edge]);
  if (draws <= 0) return;

  uint32_t rng[3];
  layout_rng_state(seed, edge, epoch, rng);
  for (int draw = 0; draw < draws; ++draw) {
    int32_t sampled = (int32_t)(layout_rng_next(rng) % (uint32_t)vertices);
    if (sampled == source) {
      float* self = embedding + (int64_t)sampled * dims;
      if (layout_distance2(current, self, dims) <= 0.0f) continue;
    }
    apply_negative(current, embedding + (int64_t)sampled * dims,
                   dims, a, b, gamma, alpha);
  }
  schedule->negative[edge] += (double)draws * schedule->negative_step[edge];
}

void umap_optimize_layout_euclidean_f32(
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
) {
  if (!embedding || !head || !tail || !epochs_per_sample) return;
  if (n_samples <= 0 || dim <= 0 || n_edges <= 0 || n_epochs <= 0) return;
  if (negative_sample_rate <= 0.0f) return;

  EdgeSchedule schedule = {0, 0, 0};
  if (!make_schedule(&schedule, epochs_per_sample, n_edges,
                     negative_sample_rate)) {
    return;
  }

  for (int epoch = 0; epoch < n_epochs; ++epoch) {
    float alpha = initial_alpha *
      (1.0f - (float)epoch / (float)n_epochs);
    for (int64_t edge = 0; edge < n_edges; ++edge) {
      if (epochs_per_sample[edge] <= 0.0f) continue;
      if (schedule.positive[edge] > (double)epoch) continue;
      update_edge(embedding, n_samples, dim, head, tail, edge, epoch,
                  &schedule, epochs_per_sample, a, b, gamma, alpha, seed);
    }
  }

  release_schedule(schedule);
}
