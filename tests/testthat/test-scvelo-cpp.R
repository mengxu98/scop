make_scanpy_mock <- function(
  n_genes = 10,
  n_cells = 20,
  n_neighbors = 4,
  n_dims = 2,
  seed = 42
) {
  set.seed(seed)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 2, rate = 1),
    nrow = n_genes,
    ncol = n_cells
  )
  beta <- 0.5
  unspliced <- beta * spliced + matrix(
    stats::rnorm(n_genes * n_cells, sd = 0.2),
    nrow = n_genes,
    ncol = n_cells
  )
  unspliced <- pmax(unspliced, 0)

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, n_neighbors))
  }

  embedding <- matrix(
    stats::rnorm(n_cells * n_dims, sd = 2),
    nrow = n_cells,
    ncol = n_dims
  )

  list(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding,
    n_genes = n_genes,
    n_cells = n_cells,
    n_neighbors = n_neighbors,
    n_dims = n_dims,
    beta = beta
  )
}

make_velocity_plot_object <- function(n_cells = 6L) {
  counts <- matrix(
    1,
    nrow = 2,
    ncol = n_cells,
    dimnames = list(paste0("gene", 1:2), paste0("cell", seq_len(n_cells)))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  embedding <- cbind(
    UMAP_1 = seq(-1, 1, length.out = n_cells),
    UMAP_2 = stats::approx(
      seq_len(6), c(-1, -0.4, 0, 0.3, 0.7, 1),
      xout = seq(1, 6, length.out = n_cells)
    )$y
  )
  rownames(embedding) <- colnames(srt)
  velocity <- cbind(
    velocity_1 = rep(0.25, n_cells),
    velocity_2 = rep(0.1, n_cells)
  )
  rownames(velocity) <- colnames(srt)
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = "RNA"
  )
  srt[["stochastic_umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = velocity,
    key = "velocity_",
    assay = "RNA"
  )
  srt
}

test_that("vector projections preserve raw and grid arrow styles", {
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  object$group <- factor(rep(c("A", "B"), each = 6L))
  original <- object
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  default_segment <- ggplot2::ggplot_build(
    ggplot2::ggplot(data.frame(x = 0, y = 0)) +
      ggplot2::geom_segment(ggplot2::aes(x, y, xend = x + 1, yend = y + 1))
  )$data[[1L]]
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    for (type in c("raw", "grid")) {
      for (group in list(NULL, "group")) {
        args <- list(
          object = object, reduction = "umap", plot_type = type,
          group.by = group, density = if (type == "raw") 1 else 0.15,
          n_neighbors = 3, min_mass = 0, return_layer = TRUE, verbose = FALSE
        )
        if (identical(plotter, PseudotimeProjectionPlot)) {
          args$time_key <- "time"
          args$method <- "gradient"
          args$label <- FALSE
        }
        layers <- do.call(plotter, args)
        expect_length(layers$velocity_layer, if (type == "raw" && !is.null(group)) 2L else 1L)
        arrow <- layers$velocity_layer[[1L]]$geom_params$arrow
        expect_identical(grid::unitType(arrow$length), "mm")
        expect_equal(as.numeric(arrow$length), 1.5)
        built <- ggplot2::ggplot_build(ggplot2::ggplot() + layers$velocity_layer)$data[[1L]]
        if (identical(plotter, VelocityPlot)) {
          expect_true(all(built$linewidth == if (type == "raw") 0.25 else 0.3))
          expect_true(all(built$alpha == if (type == "raw") 0.75 else 0.8))
          if (type == "raw") {
            expect_equal(built$xend, unname(embedding[, 1L]) + 0.25)
            expect_equal(built$yend, unname(embedding[, 2L]) + 0.1)
          }
        } else {
          expect_true(all(built$linewidth == default_segment$linewidth))
          expect_true(all(is.na(built$alpha)))
          if (type == "grid") {
            expect_identical(rownames(layers$cell_layer[[1L]]$data), colnames(object))
          }
        }
        if (type == "grid" || is.null(group)) {
          expect_true(all(built$colour == "black"))
        } else {
          colors <- thisplot::palette_colors(object$group, palette = "Chinese")
          expect_equal(built$colour, unname(colors[as.character(object$group)]))
        }
      }
    }
  }
  expect_identical(object, original)
})


test_that("vector projection calls restore caller RNG", {
  skip_if_not_installed("withr")
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  original <- object
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    args <- list(
      object = object, reduction = "umap", plot_type = "raw", density = 0.5,
      n_neighbors = 3, seed = 37, verbose = FALSE
    )
    if (identical(plotter, PseudotimeProjectionPlot)) {
      args$time_key <- "time"
      args$method <- "gradient"
      args$label <- FALSE
    }
    for (has_seed in c(TRUE, FALSE)) {
      for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        for (return_layer in c(FALSE, TRUE)) {
          expect_rng_restored(
            do.call(plotter, c(args, list(return_layer = return_layer))),
            has_seed, kind
          )
        }
        expect_rng_restored(
          expect_error(do.call(plotter, utils::modifyList(args,
            list(plot_type = "invalid")
          )), "arg"),
          has_seed, kind
        )
        expect_rng_restored(
          expect_error(do.call(plotter, c(args,
            list(theme_use = function() stop("theme failure"))
          )), "theme failure"),
          has_seed, kind
        )
      }
    }
  }
  expect_identical(object, original)
})

test_that("velocity and pseudotime stream layers remain composable", {
  skip_if_not_installed("metR")
  object <- make_velocity_plot_object(n_cells = 12L)
  object$time <- seq(0, 1, length.out = ncol(object))
  original <- object
  for (plotter in list(VelocityPlot, PseudotimeProjectionPlot)) {
    for (color in list(NULL, "black")) {
      args <- list(
        object = object, reduction = "umap", plot_type = "stream",
        density = 0.2, n_neighbors = 3, streamline_n = 4,
        streamline_L = 0.6, streamline_minL = 0.1,
        streamline_color = color, return_layer = TRUE, verbose = FALSE
      )
      if (identical(plotter, PseudotimeProjectionPlot)) {
        args$time_key <- "time"
        args$method <- "gradient"
        args$show_cells <- FALSE
      }
      layers <- do.call(plotter, args)
      expect_length(layers$velocity_layer, if (is.null(color)) 5L else 3L)
      expect_s3_class(
        ggplot2::ggplotGrob(ggplot2::ggplot() + layers), "gtable"
      )
    }
  }
  expect_identical(object, original)
})

test_that("scanpy_stochastic_embedding_cpp returns a named list with correct fields", {
  dat <- make_scanpy_mock()
  out <- scanpy_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )

  expect_type(out, "list")
  required <- c("velocity_embedding", "confidence", "velocity_length", "gamma")
  expect_true(all(required %in% names(out)))

  expect_true(is.matrix(out$velocity_embedding))
  expect_equal(dim(out$velocity_embedding), c(dat$n_cells, dat$n_dims))
  expect_true(all(is.finite(out$velocity_embedding)))

  expect_type(out$confidence, "double")
  expect_length(out$confidence, dat$n_cells)
  expect_true(all(out$confidence >= 0, na.rm = TRUE))
  expect_true(all(out$confidence <= 1, na.rm = TRUE))

  expect_type(out$velocity_length, "double")
  expect_length(out$velocity_length, dat$n_cells)
  expect_true(all(out$velocity_length >= 0, na.rm = TRUE))

  expect_type(out$gamma, "double")
  expect_length(out$gamma, dat$n_genes)
  expect_true(all(out$gamma >= 0, na.rm = TRUE))
})


test_that("scvelo handles all-zero expression", {
  n_genes <- 5
  n_cells <- 10
  spliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(0, nrow = n_genes, ncol = n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 3))
  }
  set.seed(1)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scanpy_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_true(all(out$velocity_embedding == 0))
  expect_true(all(out$confidence == 0))
  expect_true(all(out$velocity_length == 0))
  expect_true(all(out$gamma == 0))
})


test_that("scvelo estimates gamma ~ 1 when unspliced = spliced", {
  n_genes <- 20
  n_cells <- 30
  set.seed(1)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 3, rate = 1),
    nrow = n_genes, ncol = n_cells
  )
  unspliced <- spliced
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 4)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 4))
  }
  set.seed(1)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scanpy_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_equal(out$gamma, rep(1, n_genes), tolerance = 0.1)
  expect_true(mean(out$velocity_length) < 0.1)
})


test_that("scvelo handles single cell gracefully", {
  n_genes <- 5
  n_cells <- 1
  spliced <- matrix(runif(n_genes), nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(runif(n_genes), nrow = n_genes, ncol = n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 3)
  embedding <- matrix(runif(2), nrow = 1, ncol = 2)

  expect_error(
    scanpy_stochastic_embedding_cpp(
      spliced = spliced,
      unspliced = unspliced,
      knn_idx = knn_idx,
      embedding = embedding
    ),
    NA
  )
})


test_that("scvelo works with 1 neighbor", {
  n_genes <- 5
  n_cells <- 10
  set.seed(1)
  spliced <- matrix(runif(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  unspliced <- matrix(runif(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 1)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sample(candidates, 1)
  }
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scanpy_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_equal(dim(out$velocity_embedding), c(10, 2))
  expect_true(all(is.finite(out$velocity_embedding)))
})


test_that("scanpy_stochastic_embedding_cpp is deterministic", {
  dat <- make_scanpy_mock(seed = 42)
  out1 <- scanpy_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )
  out2 <- scanpy_stochastic_embedding_cpp(
    spliced = dat$spliced,
    unspliced = dat$unspliced,
    knn_idx = dat$knn_idx,
    embedding = dat$embedding
  )
  expect_equal(out1, out2)
})


test_that("scvelo gamma reflects relationship between spliced and unspliced", {
  n_genes <- 10
  n_cells <- 30
  set.seed(1)
  spliced <- matrix(runif(n_genes * n_cells, 1, 10), nrow = n_genes, ncol = n_cells)
  unspliced <- spliced
  unspliced[1:5, ] <- unspliced[1:5, ] * 2.0
  unspliced[6:10, ] <- unspliced[6:10, ] * 0.1

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 4)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, 4))
  }
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scanpy_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_true(mean(out$gamma[1:5]) > mean(out$gamma[6:10]) + 0.5)
})


test_that("scvelo rejects mismatched spliced/unspliced dimensions", {
  dat <- make_scanpy_mock()
  expect_error(
    scanpy_stochastic_embedding_cpp(
      spliced = dat$spliced[1:5, , drop = FALSE],
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx,
      embedding = dat$embedding
    ),
    "spliced and unspliced must have identical dimensions",
    ignore.case = TRUE
  )
})

test_that("scvelo rejects mismatched knn_idx rows vs cells", {
  dat <- make_scanpy_mock()
  expect_error(
    scanpy_stochastic_embedding_cpp(
      spliced = dat$spliced,
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx[1:5, , drop = FALSE],
      embedding = dat$embedding
    ),
    "knn_idx",
    ignore.case = TRUE
  )
})

test_that("scvelo rejects mismatched embedding rows vs cells", {
  dat <- make_scanpy_mock()
  expect_error(
    scanpy_stochastic_embedding_cpp(
      spliced = dat$spliced,
      unspliced = dat$unspliced,
      knn_idx = dat$knn_idx,
      embedding = dat$embedding[1:5, , drop = FALSE]
    ),
    "embedding",
    ignore.case = TRUE
  )
})


test_that("scvelo embedding output is finite and well-scaled", {
  n_genes <- 30
  n_cells <- 50
  set.seed(1)
  spliced <- matrix(
    stats::rgamma(n_genes * n_cells, shape = 2, rate = 1),
    nrow = n_genes, ncol = n_cells
  )
  unspliced <- 0.3 * spliced + matrix(
    rnorm(n_genes * n_cells, sd = 0.1),
    nrow = n_genes, ncol = n_cells
  )
  unspliced <- pmax(unspliced, 0)

  knn_idx <- matrix(NA_integer_, nrow = n_cells, ncol = 6)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    knn_idx[i, ] <- sort(sample(candidates, min(6, length(candidates))))
  }
  set.seed(99)
  embedding <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)

  out <- scanpy_stochastic_embedding_cpp(
    spliced = spliced,
    unspliced = unspliced,
    knn_idx = knn_idx,
    embedding = embedding
  )

  expect_true(all(is.finite(out$velocity_embedding)))
  expect_true(all(is.finite(out$confidence)))
  expect_true(all(is.finite(out$velocity_length)))
  expect_true(all(is.finite(out$gamma)))

  expect_gt(mean(out$velocity_length), 0)
  expect_gt(max(out$confidence), 0)
})
