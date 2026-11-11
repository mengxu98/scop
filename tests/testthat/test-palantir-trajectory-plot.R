test_that("Palantir branch selection preserves finite-probability handling", {
  dat <- data.frame(
    BranchA = c(0.4, NA_real_, Inf, -Inf, 0.5),
    BranchB = c(0.4, 0.9, -Inf, NA_real_, 0.5),
    BranchC = c(0.1, 0.2, NA_real_, NA_real_, 0.4),
    row.names = paste0("Cell", seq_len(5))
  )
  branches <- colnames(dat)
  legacy <- function(x) {
    probs <- as.matrix(x[, branches, drop = FALSE])
    probs[is.na(probs)] <- -Inf
    selected <- branches[max.col(probs, ties.method = "first")]
    has_prob <- apply(probs, 1L, function(values) any(is.finite(values)))
    selected[!has_prob] <- NA_character_
    factor(selected, levels = branches)
  }

  expect_identical(
    palantir_branch_selection(dat, branches),
    legacy(dat)
  )
})

test_that("Palantir cell colours preserve data and composable layers", {
  skip_if_not_installed("ggnewscale")
  n <- 24L
  counts <- matrix(1L, nrow = 2L, ncol = n,
    dimnames = list(c("gene1", "gene2"), paste0("cell", seq_len(n)))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  time <- seq(0, 1, length.out = n)
  embedding <- cbind(UMAP_1 = time, UMAP_2 = sin(time * pi))
  rownames(embedding) <- colnames(object)
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embedding, key = "UMAP_", assay = "RNA"
  )
  object$palantir_pseudotime <- time
  object$A_diff_potential <- 1 - time
  object$B_diff_potential <- time
  object$score <- time^2
  object$category <- factor(rep(c("A", "B"), length.out = n))
  object$category[[2L]] <- NA
  original <- object
  expected <- list(
    pseudotime = time,
    branch_selection = factor(rep(c("A_diff_potential", "B_diff_potential"), each = n / 2),
      levels = c("A_diff_potential", "B_diff_potential")
    ),
    score = time^2,
    category = object$category,
    none = NULL
  )
  args <- list(
    object = object, reduction = "umap", trajectory_method = "bin",
    smooth = FALSE, n_bins = 4, min_cells_per_bin = 1,
    return_layer = TRUE, verbose = FALSE
  )
  for (mode in names(expected)) {
    layers <- do.call(PalantirTrajectoryPlot, c(args, list(cell_color = mode)))
    built <- ggplot2::ggplot_build(ggplot2::ggplot() + layers)
    expect_s3_class(ggplot2::ggplot_gtable(built), "gtable")
    if (mode == "none") {
      expect_length(layers$cell_layer, 0L)
    } else {
      expect_equal(unname(layers$cell_layer[[1L]]$data$.cell_color), unname(expected[[mode]]))
      expect_identical(rownames(layers$cell_layer[[1L]]$data), colnames(object))
      points <- built$data[[1L]]
      expect_equal(unname(points$x), unname(embedding[, 1L]))
      expect_equal(unname(points$y), unname(embedding[, 2L]))
    }
  }
  expect_error(do.call(PalantirTrajectoryPlot, c(args, list(cell_color = "absent"))),
    "cell_color"
  )
  expect_identical(object, original)
  for (has_seed in c(TRUE, FALSE)) {
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      for (return_layer in c(FALSE, TRUE)) {
        expect_rng_restored(
          do.call(PalantirTrajectoryPlot, utils::modifyList(args,
            list(seed = 37, return_layer = return_layer)
          )),
          has_seed, kind
        )
      }
      expect_rng_restored(
        expect_error(do.call(PalantirTrajectoryPlot, utils::modifyList(args,
          list(seed = 37, trajectory_method = "invalid")
        )), "arg"),
        has_seed, kind
      )
      expect_rng_restored(
        expect_error(do.call(PalantirTrajectoryPlot, utils::modifyList(args,
          list(seed = 37, theme_use = function() stop("theme failure"))
        )), "theme failure"),
        has_seed, kind
      )
    }
  }
})
