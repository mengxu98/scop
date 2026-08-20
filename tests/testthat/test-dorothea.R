test_that("RunDorothea writes TF activity metadata when storing a new assay", {
  counts <- matrix(
    c(
      5, 0, 4,
      0, 3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  regulons <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    mor = c(1, -1),
    confidence = c("A", "B")
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = regulons,
    confidence = c("A", "B"),
    method = "ulm",
    new_assay = TRUE,
    verbose = FALSE
  )

  expect_true("dorothea" %in% SeuratObject::Assays(out))
  expect_true("regulons" %in% names(out@tools$Dorothea))
  expect_true(all(c("dorothea_TF1", "dorothea_TF2") %in% colnames(out@meta.data)))
  expect_equal(unname(out$dorothea_TF1), c(1, 2, 3))
  expect_equal(unname(out$dorothea_TF2), c(-1, -2, -3))
})

test_that("RunDorothea can skip TF activity metadata", {
  counts <- matrix(
    c(
      5, 0, 4,
      0, 3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  regulons <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    mor = c(1, -1),
    confidence = c("A", "B")
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = regulons,
    confidence = c("A", "B"),
    method = "ulm",
    new_assay = TRUE,
    add_meta = FALSE,
    verbose = FALSE
  )

  expect_true("dorothea" %in% SeuratObject::Assays(out))
  expect_false(any(grepl("^dorothea_", colnames(out@meta.data))))
  expect_true("scores" %in% names(out@tools$Dorothea))
  expect_false(out@tools$Dorothea$parameters$add_meta)
})

test_that("RunDorothea accepts RunGRN-style unsigned networks", {
  counts <- matrix(
    c(
      5, 0, 4,
      0, 3, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  grn <- data.frame(
    TF = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    importance = c(0.8, 0.4)
  )

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = grn,
    method = "ulm",
    verbose = FALSE
  )

  expect_true("dorothea" %in% SeuratObject::Assays(out))
  expect_false(out@tools$Dorothea$network_info$signed)
  expect_equal(out@tools$Dorothea$network_info$source, "custom")
  expect_equal(out@tools$Dorothea$network_info$format, "scop_grn")
  expect_equal(out@tools$Dorothea$regulons$tf, grn$TF)
  expect_equal(out@tools$Dorothea$regulons$mor, grn$importance)
  expect_equal(out@tools$Dorothea$regulons_input, grn)
})

test_that("RunDorothea loads CollecTRI with the SCOP species mapping", {
  counts <- matrix(
    c(5, 0, 4, 0, 3, 2),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  collectri_call <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    get_namespace_fun = function(pkg, fun) {
      if (identical(fun, "get_collectri")) {
        return(function(organism, split_complexes) {
          collectri_call$organism <- organism
          collectri_call$split_complexes <- split_complexes
          data.frame(
            source = c("TF1", "TF2"),
            target = c("Gene1", "Gene2"),
            mor = c(1, -1)
          )
        })
      }
      stop("unexpected namespace function: ", fun)
    },
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = "collectri",
    species = "Mus_musculus",
    method = "ulm",
    verbose = FALSE
  )

  expect_equal(collectri_call$organism, "mouse")
  expect_false(collectri_call$split_complexes)
  expect_equal(out@tools$Dorothea$network_info$source, "collectri")
  expect_equal(out@tools$Dorothea$network_info$format, "collectri")
  expect_true(out@tools$Dorothea$network_info$signed)
})

test_that("RunDorothea rejects invalid unsigned network weights before backend execution", {
  counts <- matrix(
    c(5, 0, 4, 0, 3, 2),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(...) stop("backend should not run")
  )
  expect_error(
    RunDorothea(
      srt,
      regulons = data.frame(
        TF = "TF1",
        target = "Gene1",
        importance = -1
      ),
      verbose = FALSE
    ),
    "non-negative"
  )
})

test_that("RunDorothea accepts source, regulator, and weight aliases", {
  counts <- matrix(
    c(5, 0, 4, 0, 3, 2),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(method) {
      function(...) {
        data.frame(
          source = rep(c("TF1", "TF2"), each = 3),
          condition = rep(colnames(srt), times = 2),
          score = c(1, 2, 3, -1, -2, -3)
        )
      }
    }
  )

  out <- RunDorothea(
    srt,
    regulons = data.frame(
      regulator = c("TF1", "TF2"),
      target = c("Gene1", "Gene2"),
      weight = c(0.7, 0.3)
    ),
    method = "ulm",
    verbose = FALSE
  )

  expect_equal(out@tools$Dorothea$network_info$regulator_column, "regulator")
  expect_equal(out@tools$Dorothea$network_info$weight_column, "weight")
  expect_false(out@tools$Dorothea$network_info$signed)
})

test_that("RunDorothea does not mutate the input when decoupleR fails", {
  counts <- matrix(
    c(5, 0, 4, 0, 3, 2),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:3))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  before <- list(
    assays = SeuratObject::Assays(srt),
    metadata = srt@meta.data,
    tools = srt@tools
  )
  testthat::local_mocked_bindings(
    check_r = function(...) TRUE,
    dorothea_get_run_fun = function(...) {
      function(...) stop("mock decoupleR failure")
    }
  )

  expect_error(
    RunDorothea(
      srt,
      regulons = data.frame(
        tf = "TF1",
        target = "Gene1",
        mor = 1
      ),
      verbose = FALSE
    ),
    "mock decoupleR failure"
  )
  expect_equal(SeuratObject::Assays(srt), before$assays)
  expect_equal(srt@meta.data, before$metadata)
  expect_equal(srt@tools, before$tools)
})

local({
  make_dorothea_srt <- function() {
    counts <- matrix(
      c(
        5, 0, 4, 1,
        0, 3, 2, 6
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("Gene1", "Gene2"), paste0("Cell", 1:4))
    )
    srt <- Seurat::CreateSeuratObject(
      counts = methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
    )
    srt <- Seurat::NormalizeData(srt, verbose = FALSE)
    srt$CellType <- c("A", "A", "B", "B")
    scores <- matrix(
      c(
        1, 2, -1, -2,
        -3, -1, 2, 4
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("TF1", "TF2"), colnames(srt))
    )
    srt <- dorothea_attach_assay(srt, scores, "dorothea")
    embeddings <- matrix(
      c(0, 1, 0, 1, 0, 0, 1, 1),
      ncol = 2,
      dimnames = list(colnames(srt), c("UMAP_1", "UMAP_2"))
    )
    srt[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings,
      key = "UMAP_",
      assay = "RNA"
    )
    srt@tools[["Dorothea"]] <- list(
      scores = scores,
      regulons = data.frame(
        tf = c("TF1", "TF1", "TF2"),
        target = c("Gene1", "Gene2", "Gene2"),
        mor = c(1, -1, 1),
        stringsAsFactors = FALSE
      ),
      parameters = list(assay_name = "dorothea")
    )
    srt
  }

  test_that("DorotheaPlot bar and lollipop return ggplot objects", {
    srt <- make_dorothea_srt()
    p_bar <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      plot_type = "bar",
      verbose = FALSE
    )
    p_lollipop <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      plot_type = "lollipop",
      verbose = FALSE
    )
    p_volcano <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      plot_type = "volcano",
      verbose = FALSE
    )
    expect_s3_class(p_bar, "ggplot")
    expect_s3_class(p_lollipop, "ggplot")
    expect_s3_class(p_volcano, "ggplot")
  })

  test_that("DorotheaPlot volcano plots all TFs and caps underflow p-values", {
    srt <- make_dorothea_srt()
    out <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      plot_type = "volcano",
      top_n = 1,
      nlabel = 0,
      return_data = TRUE,
      verbose = FALSE
    )
    expect_equal(sort(as.character(out$data$TF)), c("TF1", "TF2"))
    expect_equal(nrow(ggplot2::layer_data(out$plot, 3)), 2L)
    capped <- dorothea_volcano_y(data.frame(
      p_val_adj = c(0, 1e-42, 0.2),
      neglog10_p_val_adj = c(307, 42, 0.7)
    ))
    expect_equal(capped, c(50, 42, 0.7))
    labeled <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      plot_type = "volcano",
      padjustCutoff = 1,
      nlabel = 1,
      verbose = FALSE
    )
    geoms <- vapply(labeled$layers, function(layer) class(layer$geom)[[1L]], character(1))
    expect_true("GeomTextRepel" %in% geoms)
    expect_equal(nrow(labeled$layers[[match("GeomTextRepel", geoms)]]$data), 1L)
  })

  test_that("DorotheaPlot heatmap uses GroupHeatmap on the dorothea assay", {
    srt <- make_dorothea_srt()
    ht <- DorotheaPlot(
      srt,
      group.by = "CellType",
      plot_type = "heatmap",
      features = c("TF1", "TF2"),
      verbose = FALSE
    )
    expect_true(is.list(ht))
    expect_true("plot" %in% names(ht))
    expect_true("matrix_list" %in% names(ht))
    expect_equal(sort(rownames(ht$matrix_list[[1]])), c("TF1", "TF2"))
  })

  test_that("DorotheaPlot dim, stat, and targets return ggplot objects", {
    srt <- make_dorothea_srt()
    p_dim <- DorotheaPlot(
      srt,
      group.by = "CellType",
      features = "TF1",
      plot_type = "dim",
      reduction = "umap",
      verbose = FALSE
    )
    p_stat <- DorotheaPlot(
      srt,
      group.by = "CellType",
      features = c("TF1", "TF2"),
      plot_type = "stat",
      verbose = FALSE
    )
    p_targets <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      features = "TF1",
      plot_type = "targets",
      verbose = FALSE
    )
    expect_true(inherits(p_dim, "ggplot") || inherits(p_dim, "patchwork"))
    expect_true(inherits(p_stat, "ggplot") || inherits(p_stat, "patchwork"))
    expect_s3_class(p_targets, "ggplot")
    geoms <- vapply(p_targets$layers, function(layer) class(layer$geom)[[1L]], character(1))
    if ("GeomTextRepel" %in% geoms) {
      expect_true(all(
        p_targets$layers[[match("GeomTextRepel", geoms)]]$data$support != "NS"
      ))
    }
  })

  test_that("DorotheaPlot targets uses unsigned labels for unsigned networks", {
    srt <- make_dorothea_srt()
    srt@tools$Dorothea$network_info <- list(
      source = "custom",
      label = "custom TF network",
      signed = FALSE
    )
    srt@tools$Dorothea$regulons$mor <- abs(srt@tools$Dorothea$regulons$mor)
    out <- DorotheaPlot(
      srt,
      group.by = "CellType",
      group1 = "A",
      group2 = "B",
      features = "TF1",
      plot_type = "targets",
      return_data = TRUE,
      verbose = FALSE
    )
    expect_true(all(out$data$support %in% c("Higher group1", "Higher group2", "NS")))
    expect_false(any(out$data$support %in% c("Support", "Oppose")))
    expect_true(any(vapply(
      out$plot$scales$scales,
      function(scale) identical(scale$name, "Importance"),
      logical(1)
    )))
  })
})
