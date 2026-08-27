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

  generated_out <- scop:::prepare_adata_feature_metadata(
    generated,
    feature_names
  )
  external_import <- scop:::prepare_adata_feature_metadata(
    external,
    feature_names
  )
  external_export <- scop:::prepare_adata_feature_metadata(
    external,
    feature_names,
    reserve_features = TRUE
  )

  expect_false("features" %in% colnames(generated_out))
  expect_identical(external_import$features, external$features)
  expect_false("features" %in% colnames(external_export))
  expect_identical(external_export$features.metadata, external$features)
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
          adata1 <- to_adata(srt, verbose = FALSE)
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
