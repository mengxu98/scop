#ifndef UMAP_LAYOUT_H
#define UMAP_LAYOUT_H

#include <stdint.h>
#include <math.h>

static inline uint64_t layout_hash64(uint64_t state) {
  state += UINT64_C(0x9E3779B97F4A7C15);
  state = (state ^ (state >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  state = (state ^ (state >> 27)) * UINT64_C(0x94D049BB133111EB);
  return state ^ (state >> 31);
}

static inline void layout_rng_state(uint64_t seed, int64_t edge, int epoch,
                                    uint32_t state[3]) {
  uint64_t x = seed ^
    (UINT64_C(0x9E3779B97F4A7C15) * (uint64_t)(edge + 1)) ^
    (UINT64_C(0xD1B54A32D192ED03) * (uint64_t)(epoch + 1));
  x = layout_hash64(x);
  state[0] = (uint32_t)x | 2u;
  x = layout_hash64(x);
  state[1] = (uint32_t)x | 8u;
  x = layout_hash64(x);
  state[2] = (uint32_t)x | 16u;
}

static inline uint32_t layout_rng_next(uint32_t state[3]) {
  uint32_t a = state[0], b = state[1], c = state[2];
  state[0] = ((a & 4294967294u) << 12) ^ (((a << 13) ^ a) >> 19);
  state[1] = ((b & 4294967288u) << 4) ^ (((b << 2) ^ b) >> 25);
  state[2] = ((c & 4294967280u) << 17) ^ (((c << 3) ^ c) >> 11);
  return state[0] ^ state[1] ^ state[2];
}

static inline float layout_bound(float value) {
  return value > 4.0f ? 4.0f : (value < -4.0f ? -4.0f : value);
}

static inline float layout_distance2(const float* left, const float* right,
                                     int dims) {
  float total = 0.0f;
  for (int dim = 0; dim < dims; ++dim) {
    float delta = left[dim] - right[dim];
    total += delta * delta;
  }
  return total;
}

static inline float layout_pull(float distance2, float a, float b) {
  if (distance2 <= 0.0f) return 0.0f;
  float base = powf(distance2, b - 1.0f);
  return (-2.0f * a * b * base) / (a * base * distance2 + 1.0f);
}

static inline float layout_push(float distance2, float a, float b,
                                float gamma) {
  if (distance2 <= 0.0f) return 0.0f;
  return (2.0f * gamma * b) /
    ((0.001f + distance2) * (a * powf(distance2, b) + 1.0f));
}

#endif
