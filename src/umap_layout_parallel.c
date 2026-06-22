
#include "umap.h"
#include "umap_layout.h"

int loader_process_is_current(void);

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#define UMAP_SUBSTEP_CAP 64

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int arrived;
    int phase;
    int parties;
    int stopped;
} LayoutGate;

typedef struct {
    float *embedding;
    int samples;
    int dims;
    const int32_t *head;
    const int32_t *tail;
    int64_t edges;
    const float *epochs_per_sample;
    double *edge_next_positive;
    double *edge_next_negative;
    double *edge_negative_interval;
    float *partials;
    size_t embedding_values;
    int epochs;
    int slices;
    float a;
    float b;
    float gamma;
    float initial_alpha;
    float negative_sample_rate;
    uint64_t seed;
    LayoutGate *gate;
    int workers;
} LayoutShared;

typedef struct {
    int id;
    LayoutShared *shared;
} LayoutJob;

typedef struct {
    float *partials;
    double *positive;
    double *negative;
    double *interval;
    pthread_t *handles;
    LayoutJob *jobs;
} LayoutArena;

static long layout_cpu_cores(void) {
#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return si.dwNumberOfProcessors > 0 ? (long)si.dwNumberOfProcessors : 1L;
#else
    long value = sysconf(_SC_NPROCESSORS_ONLN);
    return value > 0 ? value : 1L;
#endif
}

static int layout_core_count(int requested, int64_t edges) {
    int n = requested;
    if (n <= 0) {
        long cpus = layout_cpu_cores();
        n = cpus > 0 ? (int)cpus : 1;
    }
    if (n < 1) n = 1;
    if ((int64_t)n > edges) n = (int)edges;
    return n;
}

static int layout_slice_count(int64_t edges, int samples) {
    const char *env = getenv("UMAP_LAYOUT_SUBSTEPS");
    if (env != 0) {
        long value = strtol(env, 0, 10);
        if (value >= 1) return value > (long)edges ? (int)edges : (int)value;
    }
    if (samples <= 0) return 1;
    long slices = (long)((edges + samples - 1) / samples);
    if (slices < 1) slices = 1;
    if (slices > UMAP_SUBSTEP_CAP) slices = UMAP_SUBSTEP_CAP;
    if (slices > edges) slices = (int)edges;
    return (int)slices;
}

static void gate_open(LayoutGate *gate, int parties) {
    pthread_mutex_init(&gate->mutex, 0);
    pthread_cond_init(&gate->cond, 0);
    gate->arrived = 0;
    gate->phase = 0;
    gate->parties = parties;
    gate->stopped = 0;
}

static void gate_close(LayoutGate *gate) {
    pthread_mutex_destroy(&gate->mutex);
    pthread_cond_destroy(&gate->cond);
}

static int gate_wait(LayoutGate *gate) {
    pthread_mutex_lock(&gate->mutex);
    if (gate->stopped) {
        pthread_mutex_unlock(&gate->mutex);
        return 1;
    }
    int phase = gate->phase;
    gate->arrived += 1;
    if (gate->arrived == gate->parties) {
        gate->arrived = 0;
        gate->phase += 1;
        pthread_cond_broadcast(&gate->cond);
    } else {
        while (phase == gate->phase && !gate->stopped) {
            pthread_cond_wait(&gate->cond, &gate->mutex);
        }
    }
    int stopped = gate->stopped;
    pthread_mutex_unlock(&gate->mutex);
    return stopped;
}

static void gate_stop(LayoutGate *gate) {
    pthread_mutex_lock(&gate->mutex);
    gate->stopped = 1;
    pthread_cond_broadcast(&gate->cond);
    pthread_mutex_unlock(&gate->mutex);
}

static int layout_bad_input(float *embedding,
                            int samples,
                            int dims,
                            const int32_t *head,
                            const int32_t *tail,
                            int64_t edges,
                            const float *epochs_per_sample,
                            int epochs,
                            float negative_sample_rate) {
    return embedding == 0 || head == 0 || tail == 0 || epochs_per_sample == 0 ||
        samples <= 0 || dims <= 0 || edges <= 0 || epochs <= 0 ||
        negative_sample_rate <= 0.0f;
}

static void layout_serial_path(float *embedding,
                            int samples,
                            int dims,
                            const int32_t *head,
                            const int32_t *tail,
                            int64_t edges,
                            const float *epochs_per_sample,
                            int epochs,
                            float a,
                            float b,
                            float gamma,
                            float initial_alpha,
                            float negative_sample_rate,
                            uint64_t seed) {
    umap_optimize_layout_euclidean_f32(
        embedding, samples, dims, head, tail, edges, epochs_per_sample,
        epochs, a, b, gamma, initial_alpha, negative_sample_rate, seed);
}

static void arena_clear(LayoutArena *arena) {
    free(arena->partials);
    free(arena->positive);
    free(arena->negative);
    free(arena->interval);
    free(arena->handles);
    free(arena->jobs);
    arena->partials = 0;
    arena->positive = 0;
    arena->negative = 0;
    arena->interval = 0;
    arena->handles = 0;
    arena->jobs = 0;
}

static int arena_make(LayoutArena *arena,
                      int workers,
                      size_t embedding_values,
                      int64_t edges) {
    arena->partials = (float *)malloc((size_t)workers * embedding_values * sizeof(float));
    arena->positive = (double *)malloc((size_t)edges * sizeof(double));
    arena->negative = (double *)malloc((size_t)edges * sizeof(double));
    arena->interval = (double *)malloc((size_t)edges * sizeof(double));
    arena->handles = (pthread_t *)malloc((size_t)workers * sizeof(pthread_t));
    arena->jobs = (LayoutJob *)malloc((size_t)workers * sizeof(LayoutJob));
    if (arena->partials == 0 || arena->positive == 0 || arena->negative == 0 ||
        arena->interval == 0 || arena->handles == 0 || arena->jobs == 0) {
        arena_clear(arena);
        return 0;
    }
    return 1;
}

static void arena_seed_epochs(LayoutArena *arena,
                              int64_t edges,
                              const float *epochs_per_sample,
                              float negative_sample_rate) {
    for (int64_t edge = 0; edge < edges; ++edge) {
        double every = (double)epochs_per_sample[edge];
        double neg_every = every / (double)negative_sample_rate;
        arena->positive[edge] = every;
        arena->negative[edge] = neg_every;
        arena->interval[edge] = neg_every;
    }
}

static int64_t slice_edges(int64_t edges, int slice, int slices) {
    if (edges <= (int64_t)slice) return 0;
    return (edges - (int64_t)slice + (int64_t)slices - 1) / (int64_t)slices;
}

static void slice_range(int64_t count,
                        int worker,
                        int workers,
                        int64_t *begin,
                        int64_t *end) {
    *begin = (count * (int64_t)worker) / (int64_t)workers;
    *end = (count * (int64_t)(worker + 1)) / (int64_t)workers;
}

static void value_range(size_t count,
                        int worker,
                        int workers,
                        size_t *begin,
                        size_t *end) {
    *begin = (count * (size_t)worker) / (size_t)workers;
    *end = (count * (size_t)(worker + 1)) / (size_t)workers;
}

static void layout_add_pair(float *partial,
                            const float *left,
                            const float *right,
                            int32_t left_id,
                            int32_t right_id,
                            int dims,
                            float coeff) {
    float *left_grad = partial + (int64_t)left_id * dims;
    float *right_grad = partial + (int64_t)right_id * dims;
    for (int dim = 0; dim < dims; ++dim) {
        float value = layout_bound(coeff * (left[dim] - right[dim]));
        left_grad[dim] += value;
        right_grad[dim] -= value;
    }
}

static void layout_add_noise(LayoutShared *state,
                             float *partial,
                             int64_t edge,
                             int epoch,
                             int32_t vertex,
                             const float *origin,
                             int draws) {
    uint32_t rng[3];
    layout_rng_state(state->seed, edge, epoch, rng);
    float *grad = partial + (int64_t)vertex * state->dims;
    for (int draw = 0; draw < draws; ++draw) {
        int32_t other_id = (int32_t)(layout_rng_next(rng) % (uint32_t)state->samples);
        const float *other = state->embedding + (int64_t)other_id * state->dims;
        float distance2 = layout_distance2(origin, other, state->dims);
        if (distance2 <= 0.0f && vertex == other_id) continue;
        float coeff = layout_push(distance2, state->a, state->b, state->gamma);
        if (coeff <= 0.0f) continue;
        for (int dim = 0; dim < state->dims; ++dim) {
            grad[dim] += layout_bound(coeff * (origin[dim] - other[dim]));
        }
    }
}

static void layout_visit_edge(LayoutShared *state,
                              float *partial,
                              int64_t edge,
                              int epoch) {
    if (state->epochs_per_sample[edge] <= 0.0f) return;
    if (state->edge_next_positive[edge] > (double)epoch) return;

    int32_t left_id = state->head[edge];
    int32_t right_id = state->tail[edge];
    const float *left = state->embedding + (int64_t)left_id * state->dims;
    const float *right = state->embedding + (int64_t)right_id * state->dims;
    float distance2 = layout_distance2(left, right, state->dims);
    float coeff = layout_pull(distance2, state->a, state->b);

    layout_add_pair(partial, left, right, left_id, right_id, state->dims, coeff);
    state->edge_next_positive[edge] += (double)state->epochs_per_sample[edge];

    int negatives = (int)(((double)epoch - state->edge_next_negative[edge]) /
                          state->edge_negative_interval[edge]);
    if (negatives > 0) {
        layout_add_noise(state, partial, edge, epoch, left_id, left, negatives);
        state->edge_next_negative[edge] +=
            (double)negatives * state->edge_negative_interval[edge];
    }
}

static void layout_step_edges(LayoutShared *state,
                              int worker,
                              int slice,
                              int epoch,
                              float *partial) {
    int64_t count = slice_edges(state->edges, slice, state->slices);
    int64_t begin = 0;
    int64_t end = 0;
    slice_range(count, worker, state->workers, &begin, &end);
    for (int64_t pos = begin; pos < end; ++pos) {
        int64_t edge = (int64_t)slice + pos * (int64_t)state->slices;
        layout_visit_edge(state, partial, edge, epoch);
    }
}

static void layout_apply_partials(LayoutShared *state,
                                  int worker,
                                  float alpha) {
    size_t begin = 0;
    size_t end = 0;
    value_range(state->embedding_values, worker, state->workers, &begin, &end);
    for (size_t offset = begin; offset < end; ++offset) {
        float total = 0.0f;
        for (int id = 0; id < state->workers; ++id) {
            total += state->partials[(size_t)id * state->embedding_values + offset];
        }
        state->embedding[offset] += alpha * total;
    }
}

static void *layout_run_job(void *raw) {
    LayoutJob *job = (LayoutJob *)raw;
    LayoutShared *state = job->shared;
    float *partial = state->partials + (size_t)job->id * state->embedding_values;

    for (int epoch = 0; epoch < state->epochs; ++epoch) {
        float alpha = epoch == 0
            ? state->initial_alpha
            : state->initial_alpha *
                (1.0f - (float)(epoch - 1) / (float)state->epochs);

        for (int slice = 0; slice < state->slices; ++slice) {
            memset(partial, 0, state->embedding_values * sizeof(float));
            layout_step_edges(state, job->id, slice, epoch, partial);
            if (gate_wait(state->gate)) return 0;
            layout_apply_partials(state, job->id, alpha);
            if (gate_wait(state->gate)) return 0;
        }
    }
    return 0;
}

static int launch_jobs(LayoutShared *state,
                       LayoutArena *arena,
                       LayoutGate *gate) {
    int started = 0;
    for (int id = 0; id < state->workers; ++id) {
        arena->jobs[id].id = id;
        arena->jobs[id].shared = state;
        if (pthread_create(&arena->handles[id], 0, layout_run_job, &arena->jobs[id]) != 0) {
            gate_stop(gate);
            for (int joined = 0; joined < started; ++joined) {
                pthread_join(arena->handles[joined], 0);
            }
            return 0;
        }
        started += 1;
    }
    for (int id = 0; id < state->workers; ++id) {
        pthread_join(arena->handles[id], 0);
    }
    return 1;
}

int umap_optimize_layout_euclidean_parallel_f32(
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
) {
    if (layout_bad_input(embedding, n_samples, dim, head, tail, n_edges,
                         epochs_per_sample, n_epochs, negative_sample_rate)) {
        return -1;
    }

    if (!loader_process_is_current()) {
        layout_serial_path(embedding, n_samples, dim, head, tail, n_edges,
                        epochs_per_sample, n_epochs, a, b, gamma,
                        initial_alpha, negative_sample_rate, seed);
        return -1;
    }

    int workers = layout_core_count(cores, n_edges);
    if (workers == 1) {
        layout_serial_path(embedding, n_samples, dim, head, tail, n_edges,
                        epochs_per_sample, n_epochs, a, b, gamma,
                        initial_alpha, negative_sample_rate, seed);
        return -1;
    }

    LayoutArena arena;
    arena.partials = 0;
    arena.positive = 0;
    arena.negative = 0;
    arena.interval = 0;
    arena.handles = 0;
    arena.jobs = 0;

    size_t embedding_values = (size_t)n_samples * (size_t)dim;
    if (!arena_make(&arena, workers, embedding_values, n_edges)) {
        layout_serial_path(embedding, n_samples, dim, head, tail, n_edges,
                        epochs_per_sample, n_epochs, a, b, gamma,
                        initial_alpha, negative_sample_rate, seed);
        return -1;
    }
    arena_seed_epochs(&arena, n_edges, epochs_per_sample, negative_sample_rate);

    LayoutGate gate;
    gate_open(&gate, workers);

    LayoutShared state;
    state.embedding = embedding;
    state.samples = n_samples;
    state.dims = dim;
    state.head = head;
    state.tail = tail;
    state.edges = n_edges;
    state.epochs_per_sample = epochs_per_sample;
    state.edge_next_positive = arena.positive;
    state.edge_next_negative = arena.negative;
    state.edge_negative_interval = arena.interval;
    state.partials = arena.partials;
    state.embedding_values = embedding_values;
    state.epochs = n_epochs;
    state.slices = layout_slice_count(n_edges, n_samples);
    state.a = a;
    state.b = b;
    state.gamma = gamma;
    state.initial_alpha = initial_alpha;
    state.negative_sample_rate = negative_sample_rate;
    state.seed = seed;
    state.gate = &gate;
    state.workers = workers;

    int ok = launch_jobs(&state, &arena, &gate);
    gate_close(&gate);
    arena_clear(&arena);

    if (!ok) {
        layout_serial_path(embedding, n_samples, dim, head, tail, n_edges,
                        epochs_per_sample, n_epochs, a, b, gamma,
                        initial_alpha, negative_sample_rate, seed);
        return -1;
    }
    return 0;
}
