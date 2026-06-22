
#include "knn.h"

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#define KNN_FAR 1.0e30f

typedef struct {
  const float *x;
  int rows;
  int dims;
  int width;
  int leaf;
  float *distance;
  int32_t *index;
  uint8_t *fresh;
  double *projection;
} ProjectionTree;

typedef struct {
  int first;
  int last;
  int rows;
  int dims;
  int width;
  int leaf;
  const float *x;
  float *distance;
  int32_t *index;
  uint8_t *fresh;
  uint64_t seed;
} TreeTask;

typedef struct {
  int first;
  int last;
  int width;
  int sources;
  float **source_distance;
  int32_t **source_index;
  float *distance;
  int32_t *index;
  uint8_t *fresh;
} MergeTask;

typedef struct {
  int first;
  int last;
  int rows;
  int dims;
  int width;
  const float *x;
  const int32_t *snapshot_index;
  const uint8_t *snapshot_fresh;
  const int32_t *reverse_index;
  const int *reverse_count;
  float *distance;
  int32_t *index;
  uint8_t *fresh;
  long updates;
} DescentTask;

static long machine_cores(void) {
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors > 0 ? (long)si.dwNumberOfProcessors : 1L;
#else
  long n = sysconf(_SC_NPROCESSORS_ONLN);
  return n > 0 ? n : 1L;
#endif
}

static int clamp_workers(int requested, int limit) {
  int workers = requested;
  if (workers <= 0) {
    long detected = machine_cores();
    workers = detected > 0 ? (int)detected : 1;
  }
  if (workers < 1) workers = 1;
  if (limit > 0 && workers > limit) workers = limit;
  return workers;
}

static uint64_t next_random(uint64_t *state) {
  uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

static void swap_i32(int32_t *a, int32_t *b) {
  int32_t tmp = *a;
  *a = *b;
  *b = tmp;
}

static void swap_f64(double *a, double *b) {
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

static float squared_distance(const float *a, const float *b, int dims) {
#if defined(__ARM_NEON)
  float32x4_t sum0 = vdupq_n_f32(0.0f);
  float32x4_t sum1 = sum0;
  int dim = 0;
  for (; dim + 8 <= dims; dim += 8) {
    float32x4_t d0 = vsubq_f32(vld1q_f32(a + dim), vld1q_f32(b + dim));
    float32x4_t d1 = vsubq_f32(vld1q_f32(a + dim + 4), vld1q_f32(b + dim + 4));
    sum0 = vfmaq_f32(sum0, d0, d0);
    sum1 = vfmaq_f32(sum1, d1, d1);
  }
  for (; dim + 4 <= dims; dim += 4) {
    float32x4_t d0 = vsubq_f32(vld1q_f32(a + dim), vld1q_f32(b + dim));
    sum0 = vfmaq_f32(sum0, d0, d0);
  }
  float out = vaddvq_f32(vaddq_f32(sum0, sum1));
  for (; dim < dims; ++dim) {
    float delta = a[dim] - b[dim];
    out += delta * delta;
  }
  return out;
#else
  float out = 0.0f;
  for (int dim = 0; dim < dims; ++dim) {
    float delta = a[dim] - b[dim];
    out += delta * delta;
  }
  return out;
#endif
}

static float dot_product(const float *a, const float *b, int dims) {
#if defined(__ARM_NEON)
  float32x4_t sum0 = vdupq_n_f32(0.0f);
  float32x4_t sum1 = sum0;
  int dim = 0;
  for (; dim + 8 <= dims; dim += 8) {
    sum0 = vfmaq_f32(sum0, vld1q_f32(a + dim), vld1q_f32(b + dim));
    sum1 = vfmaq_f32(sum1, vld1q_f32(a + dim + 4), vld1q_f32(b + dim + 4));
  }
  for (; dim + 4 <= dims; dim += 4) {
    sum0 = vfmaq_f32(sum0, vld1q_f32(a + dim), vld1q_f32(b + dim));
  }
  float out = vaddvq_f32(vaddq_f32(sum0, sum1));
  for (; dim < dims; ++dim) out += a[dim] * b[dim];
  return out;
#else
  float out = 0.0f;
  for (int dim = 0; dim < dims; ++dim) out += a[dim] * b[dim];
  return out;
#endif
}

static void candidate_reset(float *distance,
                            int32_t *index,
                            uint8_t *fresh,
                            int width) {
  for (int i = 0; i < width; ++i) {
    distance[i] = KNN_FAR;
    index[i] = -1;
    fresh[i] = 1;
  }
}

static int candidate_insert(float *distance,
                            int32_t *index,
                            uint8_t *fresh,
                            int width,
                            int32_t value,
                            float metric) {
  if (metric >= distance[0]) return 0;
  for (int pos = 0; pos < width; ++pos) {
    if (index[pos] == value) return 0;
  }
  distance[0] = metric;
  index[0] = value;
  fresh[0] = 1;

  int node = 0;
  for (;;) {
    int left = node * 2 + 1;
    int right = left + 1;
    int largest = node;
    if (left < width && distance[left] > distance[largest]) largest = left;
    if (right < width && distance[right] > distance[largest]) largest = right;
    if (largest == node) break;
    float dtmp = distance[node];
    distance[node] = distance[largest];
    distance[largest] = dtmp;
    int32_t itmp = index[node];
    index[node] = index[largest];
    index[largest] = itmp;
    uint8_t ftmp = fresh[node];
    fresh[node] = fresh[largest];
    fresh[largest] = ftmp;
    node = largest;
  }
  return 1;
}

static void nth_projection(double *projection,
                           int32_t *index,
                           int count,
                           int target,
                           uint64_t *rng) {
  int left = 0;
  int right = count - 1;
  while (left < right) {
    int pivot_pos = left + (int)(next_random(rng) % (uint64_t)(right - left + 1));
    double pivot = projection[pivot_pos];
    swap_f64(&projection[pivot_pos], &projection[left]);
    swap_i32(&index[pivot_pos], &index[left]);
    int store = left;
    for (int pos = left + 1; pos <= right; ++pos) {
      if (projection[pos] < pivot) {
        ++store;
        swap_f64(&projection[store], &projection[pos]);
        swap_i32(&index[store], &index[pos]);
      }
    }
    swap_f64(&projection[left], &projection[store]);
    swap_i32(&index[left], &index[store]);
    if (store == target) return;
    if (store < target) {
      left = store + 1;
    } else {
      right = store - 1;
    }
  }
}

static void add_leaf_pairs(ProjectionTree *tree, int32_t *order, int first, int last) {
  int dims = tree->dims;
  int width = tree->width;
  for (int a = first; a < last; ++a) {
    int32_t left = order[a];
    const float *left_x = tree->x + (size_t)left * dims;
    for (int b = a + 1; b < last; ++b) {
      int32_t right = order[b];
      float metric = squared_distance(left_x, tree->x + (size_t)right * dims, dims);
      candidate_insert(tree->distance + (size_t)left * width,
                       tree->index + (size_t)left * width,
                       tree->fresh + (size_t)left * width,
                       width, right, metric);
      candidate_insert(tree->distance + (size_t)right * width,
                       tree->index + (size_t)right * width,
                       tree->fresh + (size_t)right * width,
                       width, left, metric);
    }
  }
}

static void split_tree(ProjectionTree *tree,
                       int32_t *order,
                       int first,
                       int last,
                       uint64_t *rng) {
  int count = last - first;
  if (count <= tree->leaf) {
    add_leaf_pairs(tree, order, first, last);
    return;
  }

  int pick_a = (int)(next_random(rng) % (uint64_t)count);
  int pick_b = pick_a;
  while (pick_b == pick_a) {
    pick_b = (int)(next_random(rng) % (uint64_t)count);
  }

  const float *a = tree->x + (size_t)order[first + pick_a] * tree->dims;
  const float *b = tree->x + (size_t)order[first + pick_b] * tree->dims;
  float *direction = (float *)malloc((size_t)tree->dims * sizeof(float));
  if (direction == 0) {
    add_leaf_pairs(tree, order, first, last);
    return;
  }
  for (int dim = 0; dim < tree->dims; ++dim) {
    direction[dim] = b[dim] - a[dim];
  }
  for (int pos = 0; pos < count; ++pos) {
    int32_t row = order[first + pos];
    tree->projection[pos] = (double)dot_product(
      tree->x + (size_t)row * tree->dims, direction, tree->dims);
  }
  free(direction);

  nth_projection(tree->projection, order + first, count, count / 2, rng);
  int middle = first + count / 2;
  split_tree(tree, order, first, middle, rng);
  split_tree(tree, order, middle, last, rng);
}

static void shuffle_rows(int32_t *order, int rows, uint64_t *rng) {
  for (int i = 0; i < rows; ++i) order[i] = i;
  for (int i = rows - 1; i > 0; --i) {
    int j = (int)(next_random(rng) % (uint64_t)(i + 1));
    swap_i32(&order[i], &order[j]);
  }
}

static void *build_trees(void *raw) {
  TreeTask *task = (TreeTask *)raw;
  size_t total = (size_t)task->rows * (size_t)task->width;
  for (size_t i = 0; i < total; ++i) {
    task->distance[i] = KNN_FAR;
    task->index[i] = -1;
    task->fresh[i] = 1;
  }

  int32_t *order = (int32_t *)malloc((size_t)task->rows * sizeof(int32_t));
  double *projection = (double *)malloc((size_t)task->rows * sizeof(double));
  if (order == 0 || projection == 0) {
    free(order);
    free(projection);
    return 0;
  }

  ProjectionTree tree;
  tree.x = task->x;
  tree.rows = task->rows;
  tree.dims = task->dims;
  tree.width = task->width;
  tree.leaf = task->leaf;
  tree.distance = task->distance;
  tree.index = task->index;
  tree.fresh = task->fresh;
  tree.projection = projection;

  uint64_t rng = task->seed;
  for (int tree_id = task->first; tree_id < task->last; ++tree_id) {
    (void)tree_id;
    shuffle_rows(order, task->rows, &rng);
    split_tree(&tree, order, 0, task->rows, &rng);
  }
  free(order);
  free(projection);
  return 0;
}

static void *merge_candidates(void *raw) {
  MergeTask *task = (MergeTask *)raw;
  int width = task->width;
  for (int row = task->first; row < task->last; ++row) {
    float *out_distance = task->distance + (size_t)row * width;
    int32_t *out_index = task->index + (size_t)row * width;
    uint8_t *out_fresh = task->fresh + (size_t)row * width;
    candidate_reset(out_distance, out_index, out_fresh, width);
    for (int source = 0; source < task->sources; ++source) {
      const float *src_distance = task->source_distance[source] + (size_t)row * width;
      const int32_t *src_index = task->source_index[source] + (size_t)row * width;
      for (int pos = 0; pos < width; ++pos) {
        if (src_index[pos] >= 0) {
          candidate_insert(out_distance, out_index, out_fresh, width,
                           src_index[pos], src_distance[pos]);
        }
      }
    }
  }
  return 0;
}

static void *refine_neighbors(void *raw) {
  DescentTask *task = (DescentTask *)raw;
  long updates = 0;
  int width = task->width;
  int dims = task->dims;
  const int32_t *snapshot = task->snapshot_index;
  const uint8_t *is_fresh = task->snapshot_fresh;

  for (int row = task->first; row < task->last; ++row) {
    const float *row_x = task->x + (size_t)row * dims;
    float *row_distance = task->distance + (size_t)row * width;
    int32_t *row_index = task->index + (size_t)row * width;
    uint8_t *row_fresh = task->fresh + (size_t)row * width;

    for (int slot = 0; slot < width; ++slot) {
      if (!is_fresh[(size_t)row * width + slot]) continue;
      int32_t neighbor = snapshot[(size_t)row * width + slot];
      if (neighbor < 0 || neighbor >= task->rows) continue;
      const int32_t *neighbor_list = snapshot + (size_t)neighbor * width;
      for (int pos = 0; pos < width; ++pos) {
        int32_t candidate = neighbor_list[pos];
        if (candidate < 0 || candidate >= task->rows || candidate == row) continue;
        float metric = squared_distance(row_x, task->x + (size_t)candidate * dims, dims);
        updates += candidate_insert(row_distance, row_index, row_fresh,
                                    width, candidate, metric);
      }
    }

    const int32_t *reverse = task->reverse_index + (size_t)row * width;
    for (int pos = 0; pos < task->reverse_count[row]; ++pos) {
      int32_t candidate = reverse[pos];
      if (candidate < 0 || candidate >= task->rows || candidate == row) continue;
      float metric = squared_distance(row_x, task->x + (size_t)candidate * dims, dims);
      updates += candidate_insert(row_distance, row_index, row_fresh,
                                  width, candidate, metric);
    }
  }
  task->updates = updates;
  return 0;
}

static void row_chunks(int total, int workers, int worker, int *first, int *last) {
  int base = total / workers;
  int extra = total % workers;
  int offset = worker * base + (worker < extra ? worker : extra);
  int count = base + (worker < extra ? 1 : 0);
  *first = offset;
  *last = offset + count;
}

static int allocate_forest_buffers(int workers,
                                   int rows,
                                   int width,
                                   float ***distances,
                                   int32_t ***indices,
                                   uint8_t ***fresh) {
  *distances = (float **)malloc((size_t)workers * sizeof(float *));
  *indices = (int32_t **)malloc((size_t)workers * sizeof(int32_t *));
  *fresh = (uint8_t **)malloc((size_t)workers * sizeof(uint8_t *));
  if (*distances == 0 || *indices == 0 || *fresh == 0) return 0;
  for (int worker = 0; worker < workers; ++worker) {
    (*distances)[worker] = (float *)malloc((size_t)rows * width * sizeof(float));
    (*indices)[worker] = (int32_t *)malloc((size_t)rows * width * sizeof(int32_t));
    (*fresh)[worker] = (uint8_t *)malloc((size_t)rows * width);
    if ((*distances)[worker] == 0 || (*indices)[worker] == 0 ||
        (*fresh)[worker] == 0) {
      return 0;
    }
  }
  return 1;
}

static void release_forest_buffers(int workers,
                                   float **distances,
                                   int32_t **indices,
                                   uint8_t **fresh) {
  if (distances != 0) {
    for (int worker = 0; worker < workers; ++worker) free(distances[worker]);
  }
  if (indices != 0) {
    for (int worker = 0; worker < workers; ++worker) free(indices[worker]);
  }
  if (fresh != 0) {
    for (int worker = 0; worker < workers; ++worker) free(fresh[worker]);
  }
  free(distances);
  free(indices);
  free(fresh);
}

static int initialize_with_trees(int rows,
                                 int dims,
                                 int width,
                                 const float *x,
                                 int trees,
                                 int leaf,
                                 uint64_t seed,
                                 int workers,
                                 float *distance,
                                 int32_t *index,
                                 uint8_t *fresh) {
  int tree_workers = clamp_workers(workers, trees);
  float **tree_distance = 0;
  int32_t **tree_index = 0;
  uint8_t **tree_fresh = 0;
  pthread_t *handles = (pthread_t *)malloc((size_t)tree_workers * sizeof(pthread_t));
  TreeTask *tasks = (TreeTask *)malloc((size_t)tree_workers * sizeof(TreeTask));
  if (handles == 0 || tasks == 0 ||
      !allocate_forest_buffers(tree_workers, rows, width,
                               &tree_distance, &tree_index, &tree_fresh)) {
    free(handles);
    free(tasks);
    release_forest_buffers(tree_workers, tree_distance, tree_index, tree_fresh);
    return 0;
  }

  int launched = 0;
  int start_tree = 0;
  for (int worker = 0; worker < tree_workers; ++worker) {
    int first = 0;
    int last = 0;
    row_chunks(trees, tree_workers, worker, &first, &last);
    start_tree = last;
    tasks[worker].first = first;
    tasks[worker].last = last;
    tasks[worker].rows = rows;
    tasks[worker].dims = dims;
    tasks[worker].width = width;
    tasks[worker].leaf = leaf;
    tasks[worker].x = x;
    tasks[worker].distance = tree_distance[worker];
    tasks[worker].index = tree_index[worker];
    tasks[worker].fresh = tree_fresh[worker];
    tasks[worker].seed = seed ^ (UINT64_C(0x9E37) * (uint64_t)(worker + 1));
    if (first == last) continue;
    if (pthread_create(&handles[launched], 0, build_trees, &tasks[worker]) == 0) {
      launched += 1;
    } else {
      build_trees(&tasks[worker]);
    }
  }
  (void)start_tree;
  for (int worker = 0; worker < launched; ++worker) pthread_join(handles[worker], 0);

  int merge_workers = clamp_workers(workers, rows);
  pthread_t *merge_handles = (pthread_t *)malloc((size_t)merge_workers * sizeof(pthread_t));
  MergeTask *merge_tasks = (MergeTask *)malloc((size_t)merge_workers * sizeof(MergeTask));
  if (merge_handles == 0 || merge_tasks == 0) {
    free(merge_handles);
    free(merge_tasks);
    free(handles);
    free(tasks);
    release_forest_buffers(tree_workers, tree_distance, tree_index, tree_fresh);
    return 0;
  }
  launched = 0;
  for (int worker = 0; worker < merge_workers; ++worker) {
    int first = 0;
    int last = 0;
    row_chunks(rows, merge_workers, worker, &first, &last);
    merge_tasks[worker].first = first;
    merge_tasks[worker].last = last;
    merge_tasks[worker].width = width;
    merge_tasks[worker].sources = tree_workers;
    merge_tasks[worker].source_distance = tree_distance;
    merge_tasks[worker].source_index = tree_index;
    merge_tasks[worker].distance = distance;
    merge_tasks[worker].index = index;
    merge_tasks[worker].fresh = fresh;
    if (first == last) continue;
    if (pthread_create(&merge_handles[launched], 0, merge_candidates,
                       &merge_tasks[worker]) == 0) {
      launched += 1;
    } else {
      merge_candidates(&merge_tasks[worker]);
    }
  }
  for (int worker = 0; worker < launched; ++worker) pthread_join(merge_handles[worker], 0);

  free(merge_handles);
  free(merge_tasks);
  free(handles);
  free(tasks);
  release_forest_buffers(tree_workers, tree_distance, tree_index, tree_fresh);
  return 1;
}

static void build_reverse_links(int rows,
                                int width,
                                const int32_t *snapshot,
                                const uint8_t *fresh,
                                int32_t *reverse_index,
                                int *reverse_count) {
  memset(reverse_count, 0, (size_t)rows * sizeof(int));
  for (int row = 0; row < rows; ++row) {
    for (int slot = 0; slot < width; ++slot) {
      if (!fresh[(size_t)row * width + slot]) continue;
      int32_t neighbor = snapshot[(size_t)row * width + slot];
      if (neighbor >= 0 && neighbor < rows && reverse_count[neighbor] < width) {
        reverse_index[(size_t)neighbor * width + reverse_count[neighbor]] = row;
        reverse_count[neighbor] += 1;
      }
    }
  }
}

static int run_descent(int rows,
                       int dims,
                       int width,
                       const float *x,
                       int iterations,
                       int workers,
                       float *distance,
                       int32_t *index,
                       uint8_t *fresh,
                       int32_t *snapshot_index,
                       uint8_t *snapshot_fresh,
                       int32_t *reverse_index,
                       int *reverse_count) {
  int descent_workers = clamp_workers(workers, rows);
  pthread_t *handles = (pthread_t *)malloc((size_t)descent_workers * sizeof(pthread_t));
  DescentTask *tasks = (DescentTask *)malloc((size_t)descent_workers * sizeof(DescentTask));
  if (handles == 0 || tasks == 0) {
    free(handles);
    free(tasks);
    return 0;
  }

  for (int iter = 0; iter < iterations; ++iter) {
    memcpy(snapshot_index, index, (size_t)rows * width * sizeof(int32_t));
    memcpy(snapshot_fresh, fresh, (size_t)rows * width);
    memset(fresh, 0, (size_t)rows * width);
    build_reverse_links(rows, width, snapshot_index, snapshot_fresh,
                        reverse_index, reverse_count);

    int launched = 0;
    for (int worker = 0; worker < descent_workers; ++worker) {
      int first = 0;
      int last = 0;
      row_chunks(rows, descent_workers, worker, &first, &last);
      tasks[worker].first = first;
      tasks[worker].last = last;
      tasks[worker].rows = rows;
      tasks[worker].dims = dims;
      tasks[worker].width = width;
      tasks[worker].x = x;
      tasks[worker].snapshot_index = snapshot_index;
      tasks[worker].snapshot_fresh = snapshot_fresh;
      tasks[worker].reverse_index = reverse_index;
      tasks[worker].reverse_count = reverse_count;
      tasks[worker].distance = distance;
      tasks[worker].index = index;
      tasks[worker].fresh = fresh;
      tasks[worker].updates = 0;
      if (first == last) continue;
      if (pthread_create(&handles[launched], 0, refine_neighbors,
                         &tasks[worker]) == 0) {
        launched += 1;
      } else {
        refine_neighbors(&tasks[worker]);
      }
    }
    for (int worker = 0; worker < launched; ++worker) pthread_join(handles[worker], 0);
    long updates = 0;
    for (int worker = 0; worker < descent_workers; ++worker) updates += tasks[worker].updates;
    if (updates < (long)((double)rows * width * 0.0002)) break;
  }

  free(handles);
  free(tasks);
  return 1;
}

SCOP_API int knn_descent_f32(int n_cells,
                             int n_dims,
                             int k,
                             const float *X,
                             int n_trees,
                             int n_iters,
                             int leaf_size,
                             uint64_t seed,
                             int32_t *knn_idx,
                             float *knn_dist,
                             int cores) {
  int rows = n_cells;
  int dims = n_dims;
  if (rows <= 0 || dims <= 0 || k <= 0 || k >= rows) return -1;
  if (n_trees < 1) n_trees = 1;
  if (n_iters < 1) n_iters = 1;
  if (leaf_size < k + 1) leaf_size = k + 1;

  uint8_t *fresh = (uint8_t *)malloc((size_t)rows * k);
  int32_t *snapshot_index = (int32_t *)malloc((size_t)rows * k * sizeof(int32_t));
  uint8_t *snapshot_fresh = (uint8_t *)malloc((size_t)rows * k);
  int32_t *reverse_index = (int32_t *)malloc((size_t)rows * k * sizeof(int32_t));
  int *reverse_count = (int *)malloc((size_t)rows * sizeof(int));
  if (fresh == 0 || snapshot_index == 0 || snapshot_fresh == 0 ||
      reverse_index == 0 || reverse_count == 0) {
    free(fresh);
    free(snapshot_index);
    free(snapshot_fresh);
    free(reverse_index);
    free(reverse_count);
    return -2;
  }

  int workers = clamp_workers(cores, rows);
  int ok = initialize_with_trees(rows, dims, k, X, n_trees, leaf_size, seed,
                                 workers, knn_dist, knn_idx, fresh);
  if (ok) {
    ok = run_descent(rows, dims, k, X, n_iters, workers, knn_dist, knn_idx,
                     fresh, snapshot_index, snapshot_fresh, reverse_index,
                     reverse_count);
  }

  free(fresh);
  free(snapshot_index);
  free(snapshot_fresh);
  free(reverse_index);
  free(reverse_count);
  return ok ? 0 : -2;
}
