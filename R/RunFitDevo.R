#' @title Run FitDevo developmental potential scoring
#'
#' @md
#' @param object A `Seurat` object or expression matrix with genes in rows and
#' cells in columns.
#' @param assay Assay used for Seurat input.
#' @param layer Layer used for Seurat input.
#' @param features Features used for scoring. If `NULL`, the most variable
#' genes are selected.
#' @param nfeatures Number of variable genes selected when `features = NULL`.
#' @param reference.by Optional metadata column containing ordered development
#' labels. Numeric values are used directly; factors use their level order.
#' @param score.name Metadata column for the developmental potential score.
#' @param relative.name Metadata column for the relative rank.
#' @param tool_name Name used in `srt@tools`.
#' @param verbose Whether to print progress messages.
#'
#' @return A modified `Seurat` object or a result bundle for matrix input.
#'
#' @references
#' Zhang F, Yang C, Wang Y, Jiao H, Wang Z, Shen J, Li L. FitDevo:
#' accurate inference of single-cell developmental potential using
#' sample-specific gene weight. Briefings in Bioinformatics, 2022.
#' doi:10.1093/bib/bbac293.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunFitDevo(
#'   pancreas_sub,
#'   nfeatures = 300
#' )
#' FeatureDimPlot(pancreas_sub, features = "FitDevo_Score")
#' FitDevoPlot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
RunFitDevo <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nfeatures = 2000,
  reference.by = NULL,
  score.name = "FitDevo_Score",
  relative.name = "FitDevo_Relative",
  tool_name = "FitDevo",
  verbose = TRUE
) {
  input <- scop_expr_input(object, assay = assay, layer = layer)
  mat <- input$matrix
  features <- resolve_scop_features(mat, features, nfeatures)
  mat <- mat[features, , drop = FALSE]
  target <- NULL
  if (!is.null(reference.by)) {
    if (!inherits(object, "Seurat") || !reference.by %in% colnames(object@meta.data)) {
      log_message("{.arg reference.by} must be a Seurat metadata column.", message_type = "error")
    }
    target <- ordered_numeric(object@meta.data[[reference.by]])
  }
  log_message(
    "Run FitDevo scoring with {.val {length(features)}} features and {.val {ncol(mat)}} cells",
    verbose = verbose
  )
  bundle <- fitdevo_score(mat, target = target)
  bundle$parameters <- list(
    assay = input$assay,
    layer = layer,
    features = features,
    nfeatures = nfeatures,
    reference.by = reference.by
  )

  if (inherits(object, "Seurat")) {
    meta <- data.frame(
      FitDevo_Score = bundle$scores[colnames(object)],
      FitDevo_Relative = bundle$relative[colnames(object)],
      row.names = colnames(object)
    )
    colnames(meta) <- c(score.name, relative.name)
    object <- Seurat::AddMetaData(object, meta)
    object@tools[[tool_name]] <- bundle
    object <- suppressWarnings(Seurat::LogSeuratCommand(object))
    return(object)
  }
  bundle
}

#' @title Plot FitDevo results
#'
#' @description
#' Visualize FitDevo developmental potential scores with `scop` dimensional and
#' grouped summary plots.
#'
#' @md
#' @param srt A `Seurat` object processed by [RunFitDevo()].
#' @param reduction Reduction used by [FeatureDimPlot()] and [CellDimPlot()].
#' @param group.by Optional metadata column used for phenotype and score
#' distribution plots.
#' @param score.name Metadata column containing the FitDevo score.
#' @param relative.name Metadata column containing the FitDevo relative rank.
#' @param combine Whether to combine plots with `patchwork`.
#' @param nrow,ncol,byrow Layout arguments passed to [patchwork::wrap_plots()].
#' @param pt.size,pt.alpha Point size and alpha.
#' @param palette,palcolor Palette arguments for grouped plots.
#' @param theme_use,theme_args Theme arguments passed to `scop` plot helpers.
#' @param ... Additional arguments passed to [FeatureDimPlot()] and
#' [CellDimPlot()].
#'
#' @return A patchwork object or a named list of `ggplot` objects.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunFitDevo(pancreas_sub, verbose = FALSE)
#' FitDevoPlot(pancreas_sub)
FitDevoPlot <- function(
  srt,
  reduction = NULL,
  group.by = NULL,
  score.name = "FitDevo_Score",
  relative.name = "FitDevo_Relative",
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  ...
) {
  required_cols <- c(score.name, relative.name)
  missing_cols <- required_cols[!required_cols %in% colnames(srt@meta.data)]
  if (length(missing_cols) > 0L) {
    log_message("Missing FitDevo result columns: {.val {missing_cols}}.", message_type = "error")
  }
  plist <- list(
    Score = FeatureDimPlot(
      srt,
      features = score.name,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = FALSE,
      ...
    )[[1]] + ggplot2::guides(color = ggplot2::guide_colorbar(title = "FitDevo score")),
    Relative = FeatureDimPlot(
      srt,
      features = relative.name,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = FALSE,
      ...
    )[[1]] + ggplot2::guides(color = ggplot2::guide_colorbar(title = "Relative rank"))
  )
  if (!is.null(group.by)) {
    plist[["Phenotype"]] <- CellDimPlot(
      srt,
      group.by = group.by,
      reduction = reduction,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      palette = palette,
      palcolor = palcolor,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = FALSE,
      ...
    )[[1]]
    plist[["Boxplot"]] <- potency_boxplot(
      srt = srt,
      score.name = score.name,
      group.by = group.by,
      ylab = "FitDevo score",
      pt.alpha = pt.alpha,
      theme_use = theme_use,
      theme_args = theme_args,
      show_potency_axis = FALSE
    )
  }
  if (isTRUE(combine)) {
    if (is.null(nrow) && is.null(ncol) && all(c("Score", "Relative", "Phenotype", "Boxplot") %in% names(plist))) {
      return(
        (plist[["Score"]] | plist[["Relative"]]) /
          (plist[["Phenotype"]] | plist[["Boxplot"]])
      )
    }
    return(patchwork::wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow))
  }
  plist
}
