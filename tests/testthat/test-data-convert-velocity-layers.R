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
      srt_assays[["spliced"]] <- SeuratObject::CreateAssayObject(counts = counts)
      srt_assays[["unspliced"]] <- SeuratObject::CreateAssayObject(
        counts = unspliced
      )
      adata_assays <- to_adata(
        object = srt_assays,
        prepare_env = FALSE,
        verbose = FALSE
      )

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
  expect_true(all(c("spliced", "unspliced") %in% result$layers_keys))
  expect_equal(result$layers_spliced_sum, sum(counts))
  expect_equal(result$layers_unspliced_sum, sum(unspliced))
  expect_match(result$missing_warning, "nonexistent")
  expect_match(result$missing_warning, "cannot be converted")
})
