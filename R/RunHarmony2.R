#' Run Harmony algorithm
#'
#' This is a modified version of harmony::RunHarmony specifically designed for compatibility with RunSymphonyMap.
#'
#' @md
#' @param object A Seurat object.
#' @param group.by.vars A character vector specifying the batch variable name.
#' @param reduction A character string specifying the reduction to be used. Default is "pca".
#' @param dims.use An integer vector specifying the dimensions to be used. Default is 1:30.
#' @param reduction.name A character string specifying the name of the reduction to be stored in the Seurat object. Default is "Harmony".
#' @param reduction.key A character string specifying the prefix for the column names of the Harmony embeddings. Default is "Harmony_".
#' @param project.dim A logical value indicating whether to project dimension reduction loadings. Default is TRUE.
#' @param verbose A logical value indicating whether to print verbose output. Default is TRUE.
#' @param seed.use An integer specifying the random seed to be used. Default is 11.
#' @param ... Additional arguments to be passed to the [harmony::RunHarmony] function.
#'
#' @rdname RunHarmony2
#' @export
#'
#' @examples
#' panc8_sub <- RunHarmony2(
#'   panc8_sub,
#'   group.by.vars = "tech",
#'   reduction = "pca"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "pca"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "Harmony"
#' )
RunHarmony2 <- function(object, ...) {
  UseMethod(generic = "RunHarmony2", object = object)
}

#' @rdname RunHarmony2
#' @method RunHarmony2 Seurat
#' @export
RunHarmony2.Seurat <- function(
  object,
  group.by.vars,
  reduction = "pca",
  dims.use = 1:30,
  project.dim = TRUE,
  reduction.name = "Harmony",
  reduction.key = "Harmony_",
  verbose = TRUE,
  seed.use = 11L,
  ...
) {
  check_r("harmony@1.1.0")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }

  data.use <- Seurat::Embeddings(object[[reduction]])
  if (max(dims.use) > ncol(data.use)) {
    stop("trying to use more dimensions than computed")
  }

  assay <- SeuratObject::DefaultAssay(object = object[[reduction]])
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = rownames(data.use)
  )

  harmonyObject <- harmony::RunHarmony(
    data_mat = data.use[, dims.use, drop = FALSE],
    meta_data = metavars_df,
    vars_use = group.by.vars,
    verbose = verbose,
    return_object = TRUE,
    ...
  )

  harmonyEmbed <- Matrix::t(Matrix::as.matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(data.use)
  colnames(harmonyEmbed) <- paste0(
    reduction.name,
    "_",
    seq_len(ncol(harmonyEmbed))
  )

  harmonyClusters <- Matrix::t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(data.use)
  colnames(harmonyClusters) <- paste0("R", seq_len(ncol(harmonyClusters)))

  object[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = assay,
    key = reduction.key,
    misc = list(
      R = harmonyClusters,
      reduction_use = reduction,
      reduction_dims = dims.use
    )
  )

  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.name,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}
