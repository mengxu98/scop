# Tests for scvelo-compatible velocity embedding autoscaling
#
# `velocity_embedding_autoscale` replicates scvelo's
# `V_emb /= 3 * quiver_autoscale(X_emb, V_emb)` step with matplotlib's
# default `plt.subplots()` geometry. The exact scale factor is a fixed point:
# after scaling, re-running `quiver_autoscale` on the result must yield 1/3.

quiver_autoscale_ref <- function(embedding, velocity_embedding) {
  rx <- diff(range(embedding[, 1L]))
  ry <- diff(range(embedding[, 2L]))
  bbox_aspect <- ((0.88 - 0.11) * 4.8) / ((0.9 - 0.125) * 6.4)
  ratio <- bbox_aspect * rx / ry
  mean_pixel_norm <- mean(
    sqrt(velocity_embedding[, 1L]^2 + (velocity_embedding[, 2L] * ratio)^2)
  )
  sn <- max(10, sqrt(nrow(embedding)))
  1.8 * sn * mean_pixel_norm / (1.1 * rx)
}

test_that("velocity_embedding_autoscale preserves direction", {
  set.seed(42)
  X <- cbind(
    UMAP_1 = stats::rnorm(50, sd = 4),
    UMAP_2 = stats::rnorm(50, sd = 3)
  )
  V <- cbind(
    velocity_1 = stats::rnorm(50, sd = 0.5),
    velocity_2 = stats::rnorm(50, sd = 0.5)
  )
  out <- scop:::velocity_embedding_autoscale(X, V)
  cosv <- rowSums(out * V) / (sqrt(rowSums(out^2)) * sqrt(rowSums(V^2)))
  expect_equal(cosv, rep(1, nrow(X)), tolerance = 1e-12)
})

test_that("velocity_embedding_autoscale is scale-equivariant", {
  set.seed(42)
  X <- cbind(UMAP_1 = stats::rnorm(30, sd = 3), UMAP_2 = stats::rnorm(30, sd = 2))
  V <- cbind(velocity_1 = stats::rnorm(30), velocity_2 = stats::rnorm(30))
  out1 <- scop:::velocity_embedding_autoscale(X, V)
  out2 <- scop:::velocity_embedding_autoscale(X, V * 2.5)
  expect_equal(out1, out2, tolerance = 1e-12)
})

test_that("velocity_embedding_autoscale reaches the scvelo fixed point", {
  set.seed(42)
  X <- cbind(UMAP_1 = stats::rnorm(100, sd = 5), UMAP_2 = stats::rnorm(100, sd = 4))
  V <- cbind(velocity_1 = stats::rnorm(100), velocity_2 = stats::rnorm(100))
  out <- scop:::velocity_embedding_autoscale(X, V)
  expect_equal(
    quiver_autoscale_ref(X, out),
    1 / 3,
    tolerance = 1e-12
  )
})
