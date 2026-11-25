test_that("AnnData layer alignment uses identifiers rather than shape", {
  align <- getFromNamespace("align_adata_layer", "scop")
  genes <- c("g1", "g2", "g3")
  cells <- c("c1", "c2", "c3")
  mat <- Matrix::Matrix(
    matrix(seq_len(9L), 3L, dimnames = list(genes, cells)),
    sparse = TRUE
  )

  reorder <- function(matrix, features = genes, cells_use = cells) {
    align(matrix, features, cells_use, "RNA", "spliced", "counts")
  }
  expect_identical(reorder(mat), mat)
  expect_identical(reorder(mat[3:1, c(2, 3, 1)]), mat)
  expect_identical(reorder(mat, genes[2], cells[3]), mat[2, 3, drop = FALSE])
  expect_s4_class(reorder(mat[3:1, 3:1]), "sparseMatrix")

  bad <- mat
  rownames(bad)[3] <- "g4"
  colnames(bad)[3] <- "c4"
  expect_error(reorder(bad), "missing 1 requested[[:space:]]+feature ID")
  rownames(bad)[3] <- "g2"
  expect_error(reorder(bad), "non-missing,[[:space:]]+non-empty and unique")
  expect_error(reorder(unname(mat)), "non-missing,[[:space:]]+non-empty and unique")
})

test_that("srt_to_adata converts velocity matrices stored as assay layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("reticulate")
  skip_if_not_installed("callr")

  set.seed(1)
  genes <- paste0("g", seq_len(20))
  cells <- paste0("c", seq_len(12))
  counts <- Matrix::sparseMatrix(
    i = rep(seq_along(genes), times = length(cells)),
    j = rep(seq_along(cells), each = length(genes)),
    x = stats::rpois(length(genes) * length(cells), 5) + 1,
    dimnames = list(genes, cells)
  )
  unspliced <- counts * 0.5

  package_path <- getNamespaceInfo(asNamespace("scop"), "path")
  source_tree <- file.exists(file.path(package_path, ".Rbuildignore"))
  result <- callr::r(
    function(counts, unspliced, package_path, source_tree, libpath) {
      if (!source_tree) {
        libpath <- unique(c(dirname(package_path), libpath))
      }
      .libPaths(libpath)
      if (source_tree) {
        pkgload::load_all(package_path, quiet = TRUE)
      }
      to_adata <- getExportedValue("scop", "srt_to_adata")

      error_message <- NULL
      available <- tryCatch(
        {
          reticulate::import("anndata", convert = FALSE)
          reticulate::import("scipy.sparse", convert = FALSE)
          TRUE
        },
        error = function(e) {
          error_message <<- paste(class(e)[1], conditionMessage(e))
          FALSE
        }
      )
      if (!isTRUE(available)) {
        out <- list(available = FALSE)
        attr(out, "error_message") <- error_message
        return(out)
      }

      layer_keys <- function(adata) {
        keys <- reticulate::iterate(adata$layers$keys())
        vapply(
          keys,
          function(key) as.character(reticulate::py_to_r(key)),
          character(1)
        )
      }
      layer_sum <- function(adata, layer) {
        value <- adata$layers$get(layer)
        if (inherits(value, "python.builtin.NoneType")) {
          return(NA_real_)
        }
        as.numeric(reticulate::py_to_r(value$sum()))
      }

      srt <- Seurat::CreateSeuratObject(counts = counts)

      srt_assays <- srt
      srt_assays[["spliced"]] <- SeuratObject::CreateAssayObject(
        counts = counts[rev(rownames(counts)), , drop = FALSE]
      )
      srt_assays[["unspliced"]] <- SeuratObject::CreateAssayObject(
        counts = unspliced
      )
      adata_assays <- to_adata(
        object = srt_assays,
        prepare_env = FALSE,
        verbose = FALSE
      )

      adata_one <- to_adata(
        object = srt_assays,
        features = rownames(counts)[2],
        prepare_env = FALSE,
        verbose = FALSE
      )

      srt_split <- srt_assays
      srt_split[["RNA"]] <- split(
        srt_split[["RNA"]],
        f = rep(c("even", "odd"), length.out = ncol(counts))
      )
      adata_split <- to_adata(
        object = srt_split,
        prepare_env = FALSE,
        verbose = FALSE
      )
      split_cells <- as.character(reticulate::py_to_r(adata_split$obs_names$tolist()))
      split_genes <- as.character(reticulate::py_to_r(adata_split$var_names$tolist()))
      split_layer <- reticulate::py_to_r(adata_split$layers$get("spliced")$toarray())

      srt_layers <- srt
      srt_layers[["RNA"]]$spliced <- counts
      srt_layers[["RNA"]]$unspliced <- unspliced
      adata_layers <- to_adata(
        object = srt_layers,
        prepare_env = FALSE,
        verbose = FALSE
      )

      warnings_seen <- character(0)
      invisible(withCallingHandlers(
        to_adata(
          object = srt_layers,
          assay_y = c("spliced", "unspliced", "nonexistent"),
          prepare_env = FALSE,
          verbose = TRUE
        ),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ))

      list(
        available = TRUE,
        assays_keys = layer_keys(adata_assays),
        assays_spliced = reticulate::py_to_r(adata_assays$layers$get("spliced")$toarray()),
        one_gene = as.character(reticulate::py_to_r(adata_one$var_names$tolist())),
        one_spliced = reticulate::py_to_r(adata_one$layers$get("spliced")$toarray()),
        split_cells = split_cells,
        split_genes = split_genes,
        split_layer = split_layer,
        layers_keys = layer_keys(adata_layers),
        layers_spliced_sum = layer_sum(adata_layers, "spliced"),
        layers_unspliced_sum = layer_sum(adata_layers, "unspliced"),
        missing_warning = paste(warnings_seen, collapse = " ")
      )
    },
    args = list(counts, unspliced, package_path, source_tree, .libPaths())
  )
  skip_if_not(
    isTRUE(result$available),
    paste0("Python anndata environment unavailable", ifelse(
      !is.null(attr(result, "error_message")),
      paste0(": ", attr(result, "error_message")),
      ""
    ))
  )

  expect_true(all(c("spliced", "unspliced") %in% result$assays_keys))
  expect_equal(result$assays_spliced, unname(t(as.matrix(counts))))
  expect_identical(result$one_gene, genes[2])
  expect_equal(result$one_spliced, unname(t(as.matrix(counts[2, , drop = FALSE]))))
  expect_setequal(result$split_cells, cells)
  expect_identical(result$split_genes, genes)
  expect_equal(
    result$split_layer,
    unname(t(as.matrix(counts[genes, result$split_cells, drop = FALSE])))
  )
  expect_true(all(c("spliced", "unspliced") %in% result$layers_keys))
  expect_equal(result$layers_spliced_sum, sum(counts))
  expect_equal(result$layers_unspliced_sum, sum(unspliced))
  expect_match(result$missing_warning, "nonexistent")
  expect_match(result$missing_warning, "cannot be converted")
})
