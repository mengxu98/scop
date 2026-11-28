test_that("split FeatureDimPlot preserves shared point scales and missing-panel fallback", {
  counts <- rbind(
    gene1 = c(0, 0, 0, 0, 2, 3, 4, 5, 1, 3, 5, 8),
    gene2 = c(1, 4, 9, 5, 3, 4, 3, 4, 3, 5, 11, 21)
  )
  cells <- paste0("cell", seq_len(ncol(counts)))
  colnames(counts) <- cells
  obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  embedding <- cbind(UMAP_1 = seq_along(cells), UMAP_2 = seq_along(cells) %% 5)
  rownames(embedding) <- cells
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "UMAP_", assay = "RNA"
  )
  obj$batch <- factor(rep(c("A", "B", "C"), each = 4L), levels = c("B", "A", "C"))

  for (mode in list("feature", "all", NULL)) {
    plots <- FeatureDimPlot(
      obj, rownames(counts), layer = "counts", reduction = "umap",
      split.by = "batch", keep_scale = mode, upper_quantile = 0.8,
      lower_cutoff = 1.25, raster = FALSE, combine = FALSE, verbose = FALSE
    )
    expect_identical(names(plots), c("B:gene1", "A:gene1", "C:gene1",
      "B:gene2", "A:gene2", "C:gene2"))
    limits <- function(plot) plot$scales$get_scales("colour")$limits
    expect_null(limits(plots[["A:gene1"]]))
    expect_true(all(is.na(plots[["A:gene1"]]$data$value)))
    for (name in setdiff(names(plots), "A:gene1")) {
      rows <- strsplit(name, ":", fixed = TRUE)[[1L]]
      pool <- if (is.null(mode)) {
        counts[rows[[2L]], obj$batch == rows[[1L]]]
      } else if (identical(mode, "feature")) {
        counts[rows[[2L]], ]
      } else {
        counts
      }
      expected <- c(1.25, unname(stats::quantile(pool[pool > 0], .8)) + 0.001)
      expect_equal(limits(plots[[name]]), expected)
      expect_true(all(plots[[name]]$data$value >= expected[1] &
        plots[[name]]$data$value <= expected[2], na.rm = TRUE))
    }
  }
})
