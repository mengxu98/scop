test_that("the Seurat marker reference follows concrete and inherited classes", {
  dispatch <- environment(seurat_reference_find_all_markers)$FindMarkers
  record_method <- function(method) {
    force(method)
    function(object, ...) list(method = method, object = object, args = list(...))
  }
  testthat::local_mocked_bindings(
    FindMarkers = function(...) stop("The generic must not select the reference"),
    FindMarkers.Seurat = record_method("Seurat"),
    FindMarkers.Assay = record_method("Assay"),
    FindMarkers.StdAssay = record_method("StdAssay"),
    FindMarkers.SCTAssay = record_method("SCTAssay"),
    FindMarkers.DimReduc = record_method("DimReduc"),
    FindMarkers.default = record_method("default"),
    .package = "Seurat"
  )
  counts <- Matrix::Matrix(
    matrix(
      seq_len(12),
      nrow = 3,
      dimnames = list(paste0("g", 1:3), paste0("c", 1:4))
    ),
    sparse = TRUE
  )
  objects <- list(
    Seurat = SeuratObject::CreateSeuratObject(counts),
    Assay = SeuratObject::CreateAssayObject(counts),
    StdAssay = SeuratObject::CreateAssay5Object(counts),
    SCTAssay = methods::as(SeuratObject::CreateAssayObject(counts), "SCTAssay"),
    DimReduc = SeuratObject::CreateDimReducObject(
      embeddings = matrix(
        seq_len(8),
        nrow = 4,
        dimnames = list(colnames(counts), c("PC_1", "PC_2"))
      ),
      assay = "RNA",
      key = "PC_"
    ),
    default = as.matrix(counts)
  )
  args <- list(cells.1 = c("c1", "c2"), cells.2 = c("c3", "c4"), random.seed = 7L)
  for (class in names(objects)) {
    object <- objects[[class]]
    result <- do.call(dispatch, c(list(object = object), args))
    expect_identical(result$method, class)
    expect_identical(result$object, object)
    expect_identical(result$args, args)
  }
})

test_that("the sampled Seurat all-marker reference matches pairwise v3 and v5 assays", {
  set.seed(20261005)
  counts <- Matrix::Matrix(
    matrix(
      stats::rpois(24L * 32L, lambda = 3),
      nrow = 24L,
      dimnames = list(paste0("g", seq_len(24L)), paste0("c", seq_len(32L)))
    ),
    sparse = TRUE
  )
  groups <- rep(c("A", "B", "C", "D"), each = 8L)
  for (assay_version in c("v3", "v5")) {
    withr::local_options(Seurat.object.assay.version = assay_version)
    object <- SeuratObject::CreateSeuratObject(counts)
    object <- seurat_reference_method("NormalizeData", "Seurat", object, verbose = FALSE)
    SeuratObject::Idents(object) <- groups
    args <- list(
      object = object,
      logfc.threshold = 0,
      min.pct = 0,
      only.pos = FALSE,
      max.cells.per.ident = 5L,
      random.seed = 42L,
      exact = FALSE,
      verbose = FALSE
    )
    expected <- do.call(rbind, lapply(sort(unique(groups)), function(group) {
      markers <- do.call(
        get("FindMarkers.Seurat", asNamespace("Seurat")),
        c(args, list(ident.1 = group))
      )
      markers$cluster <- factor(group, levels = sort(unique(groups)))
      markers$gene <- rownames(markers)
      markers
    }))
    expect_warning(
      actual <- do.call(seurat_reference_find_all_markers, c(args, list(return.thresh = Inf))),
      NA
    )
    expect_gt(nrow(actual), 0L)
    expect_setequal(as.character(unique(actual$cluster)), unique(groups))
    expected <- expected[order(expected$cluster, expected$gene), , drop = FALSE]
    actual <- actual[order(actual$cluster, actual$gene), , drop = FALSE]
    rownames(expected) <- rownames(actual) <- NULL
    expect_equal(actual, expected, tolerance = 1e-12, info = assay_version)
  }
})
