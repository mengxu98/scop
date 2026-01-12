#' @title Run Harmony algorithm
#'
#' @description
#' This is a modified version of [harmony::RunHarmony] specifically designed for compatibility with [RunSymphonyMap].
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunUMAP2
#' @param object A Seurat object.
#' @param group.by.vars The batch variable name.
#' @param dims.use The dimensions to be used.
#' Default is `1:30`.
#' @param reduction.name The name of the reduction to be stored in the Seurat object.
#' Default is `"Harmony"`.
#' @param reduction.key The prefix for the column names of the Harmony embeddings.
#' Default is `"Harmony_"`.
#' @param project.dim Whether to project dimension reduction loadings.
#' Default is `TRUE`.
#' @param ... Additional arguments to be passed to [harmony::RunHarmony].
#'
#' @rdname RunHarmony2
#' @export
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- standard_scop(panc8_sub)
#' panc8_sub <- RunHarmony2(
#'   panc8_sub,
#'   group.by.vars = "tech",
#'   reduction = "pca"
#' )
#'
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "pca"
#' )
#'
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "Harmony"
#' )
#'
#' panc8_sub <- standard_scop(
#'   panc8_sub,
#'   prefix = "Harmony",
#'   linear_reduction = "Harmony"
#' )
#'
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "StandardpcaUMAP2D"
#' )
#'
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "HarmonyUMAP2D"
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
    assay = NULL,
    reduction = "pca",
    dims.use = 1:30,
    project.dim = TRUE,
    reduction.name = "Harmony",
    reduction.key = "Harmony_",
    verbose = TRUE,
    seed.use = 11,
    ...) {
  check_r("harmony", verbose = FALSE)
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  if (length(dims.use) == 1) {
    log_message(
      "Only specified one dimension in {.arg dims.use}",
      message_type = "error"
    )
  }
  reduction <- DefaultReduction(
    object,
    pattern = reduction
  )
  data_use <- Seurat::Embeddings(object, reduction = reduction)
  if (max(dims.use) > ncol(data_use)) {
    log_message(
      "Trying to use more dimensions than computed",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = rownames(data_use)
  )

  harmonyObject <- harmony::RunHarmony(
    data_mat = data_use[, dims.use, drop = FALSE],
    meta_data = metavars_df,
    vars_use = group.by.vars,
    verbose = verbose,
    return_object = TRUE,
    ...
  )

  harmonyEmbed <- Matrix::t(as_matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(data_use)
  colnames(harmonyEmbed) <- paste0(
    reduction.name,
    "_",
    seq_len(ncol(harmonyEmbed))
  )

  harmonyClusters <- Matrix::t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(data_use)
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
