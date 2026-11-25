test_that("feature metadata only drops generated identifier columns", {
  feature_names <- c("g1", "g2", "g3")
  generated <- data.frame(
    features = feature_names,
    symbol = c("A", "B", "C"),
    row.names = feature_names
  )
  external <- data.frame(
    features = c("protein_coding", "lncRNA", "protein_coding"),
    symbol = c("A", "B", "C"),
    row.names = feature_names
  )

  generated_out <- getFromNamespace("prepare_adata_feature_metadata", "scop")(
    generated,
    feature_names
  )
  external_import <- getFromNamespace("prepare_adata_feature_metadata", "scop")(
    external,
    feature_names
  )
  external_export <- getFromNamespace("prepare_adata_feature_metadata", "scop")(
    external,
    feature_names,
    reserve_features = TRUE
  )

  expect_false("features" %in% colnames(generated_out))
  expect_identical(external_import$features, external$features)
  expect_false("features" %in% colnames(external_export))
  expect_identical(external_export$features.metadata, external$features)
})

test_that("AnnData observation metadata retains colliding normalized columns", {
  cells <- c("cell3", "cell1", "cell2")
  obs <- data.frame(
    `a-b` = c("C", "A", "B"),
    `a b` = c("three", "one", "two"),
    a.b.1 = c(30, 10, 20),
    keep_me = c(TRUE, FALSE, TRUE),
    row.names = cells,
    check.names = FALSE
  )
  adata <- structure(
    list(
      X = matrix(1:6, nrow = 3L),
      obs = obs,
      obs_names = cells,
      obs_keys = function() names(obs),
      var_names = c("gene1", "gene2"),
      var_keys = function() character()
    ),
    class = "InMemoryAnnData"
  )
  srt <- adata_to_srt(adata, prepare_env = FALSE, verbose = FALSE)
  metadata <- srt[[]]
  expected_names <- c("a.b", "a.b.2", "a.b.1", "keep_me")
  expect_identical(rownames(metadata), cells)
  expect_identical(names(metadata)[names(metadata) %in% expected_names], expected_names)
  expect_identical(
    unname(as.list(metadata[, expected_names, drop = FALSE])),
    unname(as.list(obs))
  )
})

test_that("srt <-> adata roundtrip keeps var column names unique", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("callr")

  set.seed(1)
  counts <- Matrix::rsparsematrix(20, 12, density = 0.2)
  counts@x <- abs(round(counts@x * 10)) + 1
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("c", seq_len(ncol(counts)))

  srt <- Seurat::CreateSeuratObject(counts)

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  source_tree <- file.exists(file.path(package_path, ".Rbuildignore"))
  result <- callr::r(
    function(srt, package_path, source_tree, libpath) {
      if (!source_tree) {
        libpath <- unique(c(dirname(package_path), libpath))
      }
      .libPaths(libpath)
      if (source_tree) {
        pkgload::load_all(package_path, quiet = TRUE)
      }

      to_adata <- getExportedValue("scop", "srt_to_adata")
      to_srt <- getExportedValue("scop", "adata_to_srt")

      conv_err <- NULL
      out <- tryCatch(
        {
          adata1 <- to_adata(srt, prepare_env = FALSE, verbose = FALSE)
          srt2 <- to_srt(adata1, prepare_env = FALSE, verbose = FALSE)
          adata2 <- to_adata(srt2, prepare_env = FALSE, verbose = FALSE)

          rna_assay <- srt2[["RNA"]]
          rna_meta_cols <- if (inherits(rna_assay, "Assay5")) {
            colnames(rna_assay@meta.data)
          } else {
            colnames(rna_assay@meta.features)
          }

          var_cols <- as.character(reticulate::py_to_r(
            adata2$var$columns$tolist()
          ))
          list(
            available = TRUE,
            var_cols = var_cols,
            features_col = as.character(reticulate::py_to_r(
              adata2$var[["features"]]$tolist()
            )),
            assay_meta_cols = rna_meta_cols
          )
        },
        error = function(e) {
          conv_err <<- conditionMessage(e)
          list(available = FALSE)
        }
      )
      attr(out, "error_message") <- conv_err
      out
    },
    args = list(srt, package_path, source_tree, .libPaths())
  )
  skip_if_not(
    isTRUE(result$available),
    paste0("Python anndata environment unavailable", ifelse(
      !is.null(attr(result, "error_message")),
      paste0(": ", attr(result, "error_message")),
      ""
    ))
  )

  expect_identical(anyDuplicated(result$var_cols), 0L)
  expect_identical(sum(result$var_cols == "features"), 1L)
  expect_identical(result$features_col, rownames(srt))
  expect_false("features" %in% result$assay_meta_cols)
})

test_that("adata_to_srt(prepare_env = FALSE) does not prepare the Scanpy environment", {
  n_prepare <- 0L
  testthat::local_mocked_bindings(
    .package = "scop",
    PrepareEnv = function(...) {
      n_prepare <<- n_prepare + 1L
      invisible(NULL)
    }
  )
  adata <- structure(
    list(X = structure(list(), class = "python.builtin.object")),
    class = "python.builtin.object"
  )
  try(
    adata_to_srt(adata, prepare_env = FALSE, verbose = FALSE),
    silent = TRUE
  )
  expect_identical(n_prepare, 0L)

  try(
    getFromNamespace("adata_matrix_to_r", "scop")(adata$X),
    silent = TRUE
  )
  expect_identical(n_prepare, 0L)
})

test_that("AnnData keys cannot collide in the Seurat component namespace", {
  layer <- matrix(101:112, nrow = 4L)
  embedding <- matrix(1:8, nrow = 4L)
  adata <- structure(
    list(
      X = matrix(1:12, nrow = 4L),
      layers = list(RNA = layer, `signal-A` = layer, signal.A = layer),
      obsm = list(shared = embedding, `graph-key` = embedding),
      obsp = list(shared = diag(4L), graph.key = diag(4L))
    ),
    class = "InMemoryAnnData"
  )
  err <- tryCatch(
    adata_to_srt(adata, prepare_env = FALSE, verbose = FALSE),
    error = identity
  )
  expect_s3_class(err, "error")
  message <- conditionMessage(err)
  expect_match(message, "Conflicting AnnData keys")
  expect_match(message, "RNA.*X.*RNA.*layers.*RNA")
  expect_match(message, "signal[.]A.*layers.*signal-A.*layers.*signal[.]A")
  expect_match(message, "shared.*obsm.*shared.*obsp.*shared")
  expect_match(message, "graph[.]key.*obsm.*graph-key.*obsp.*graph[.]key")
})
