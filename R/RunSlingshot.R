#' RunSlingshot
#'
#' Runs the Slingshot algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param group.by The variable to group the cells by.
#' @param reduction The reduction technique to use for dimensionality reduction. Default is NULL, which uses the default reduction for the Seurat object.
#' @param dims The dimensions to use for the Slingshot algorithm. Default is NULL, which uses first two dimensions.
#' @param start The starting group for the Slingshot algorithm. Default is NULL.
#' @param end The ending group for the Slingshot algorithm. Default is NULL.
#' @param prefix The prefix to add to the column names of the resulting pseudotime variable. Default is NULL.
#' @param reverse Logical value indicating whether to reverse the pseudotime variable. Default is FALSE.
#' @param align_start Logical value indicating whether to align the starting pseudotime values at the maximum pseudotime. Default is FALSE.
#' @param show_plot Logical value indicating whether to show the dimensionality plot. Default is TRUE.
#' @param lineage_palette The color palette to use for the lineages in the plot. Default is "Dark2".
#' @param seed The random seed to use for reproducibility. Default is 11.
#' @param ... Additional arguments to be passed to the \code{\link[slingshot]{slingshot}} function.
#'
#' @seealso \code{\link{CellDimPlot}} \code{\link{RunDynamicFeatures}}
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "PCA",
#'   dims = 1:10
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_span = 0.1
#' )
#'
#' # 3D lineage
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "StandardpcaUMAP3D"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP",
#'   lineages = paste0("Lineage", 1:2),
#'   lineages_span = 0.1,
#'   lineages_trim = c(0.05, 0.95)
#' )
RunSlingshot <- function(
    srt,
    group.by,
    reduction = NULL,
    dims = NULL,
    start = NULL,
    end = NULL,
    prefix = NULL,
    reverse = FALSE,
    align_start = FALSE,
    show_plot = TRUE,
    lineage_palette = "Dark2",
    seed = 11,
    ...) {
  if (missing(group.by)) {
    log_message(
      "group.by is missing",
      message_type = "error"
    )
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (is.null(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }

  if (min(table(srt[[group.by]])) < 2) {
    celltypes <- names(which(table(srt[[group.by]]) < 2))
    log_message(
      paste(celltypes, collapse = ", "),
      " have less than 2 cells. Removed from the analysis.",
      message_type = "warning"
    )
    celltypes <- setdiff(names(table(srt[[group.by]])), celltypes)
    srt <- srt[, select_cells(srt, celltypes, group.by)]
  }

  srt_sub <- srt[, !is.na(srt[[group.by, drop = TRUE]])]

  if (is.null(dims)) {
    dims <- 1:2
  }

  set.seed(seed)
  sl <- slingshot::slingshot(
    data = as.data.frame(
      srt_sub[[reduction]]@cell.embeddings[, dims]
    ),
    clusterLabels = as.character(
      srt_sub[[group.by, drop = TRUE]]
    ),
    start.clus = start,
    end.clus = end,
    ...
  )

  srt@tools[[paste("Slingshot", group.by, reduction, sep = "_")]] <- sl
  df <- as.data.frame(slingshot::slingPseudotime(sl))
  colnames(df) <- paste0(prefix, colnames(df))
  if (isTRUE(reverse)) {
    if (isTRUE(align_start)) {
      df <- apply(df, 2, function(x) max(x, na.rm = TRUE) - x)
    } else {
      df <- max(df, na.rm = TRUE) - df
    }
  }
  srt <- Seurat::AddMetaData(srt, metadata = df)
  srt <- Seurat::AddMetaData(
    srt,
    metadata = slingshot::slingBranchID(sl),
    col.name = paste0(prefix, "BranchID")
  )

  if (isTRUE(show_plot)) {
    if (
      ncol(srt[[reduction]]@cell.embeddings) == 2 ||
        ncol(srt[[reduction]]@cell.embeddings) > 3
    ) {
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scop(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 2, type = "lineages", col = "black")
      # plot(srt[[reduction]]@cell.embeddings, col = palette_scop(srt[[group.by, drop = TRUE]], matched = TRUE), asp = 1, pch = 16)
      # lines(slingshot::SlingshotDataSet(sl), lwd = 3, col = 1:length(slingshot::SlingshotDataSet(sl)@lineages))
      p <- CellDimPlot(
        srt,
        group.by = group.by,
        reduction = reduction,
        dims = c(1, 2),
        lineages = colnames(df)
      )
      print(p)
    } else if (ncol(srt[[reduction]]@cell.embeddings) == 3) {
      p <- CellDimPlot3D(
        srt,
        group.by = group.by,
        reduction = reduction,
        lineages = colnames(df)
      )
      print(p)
    }
  }
  return(srt)
}
