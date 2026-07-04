#' @title Run FWP feature-weight phenotype scoring
#'
#' @md
#' @inheritParams RunFitDevo
#' @param phenotype.by Metadata column used to train feature weights.
#' @param positive Positive phenotype label. If `NULL`, the last ordered label
#' is used.
#' @param weights Optional named vector of precomputed feature weights.
#' @param score.name Metadata column for the FWP score.
#'
#' @return A modified `Seurat` object or a result bundle for matrix input.
#'
#' @references
#' Wang Q, Song JJ, Zhang F. Feature-weight based measurement of cancerous
#' transcriptome using cohort-wide and sample-specific information. Cellular
#' Oncology, 2024. doi:10.1007/s13402-023-00879-6.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub$IsEndocrine <- pancreas_sub$CellType == "Endocrine"
#' pancreas_sub <- RunFWP(
#'   pancreas_sub,
#'   phenotype.by = "IsEndocrine",
#'   nfeatures = 300
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "FWP_Score"
#' ) + CellDimPlot(
#'   pancreas_sub,
#'   group.by = "IsEndocrine",
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
RunFWP <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  nfeatures = 2000,
  phenotype.by = NULL,
  positive = NULL,
  weights = NULL,
  score.name = "FWP_Score",
  tool_name = "FWP",
  verbose = TRUE
) {
  input <- scop_expr_input(object, assay = assay, layer = layer)
  mat <- input$matrix
  features <- resolve_scop_features(mat, features, nfeatures)
  mat <- mat[features, , drop = FALSE]
  y <- NULL
  if (is.null(weights)) {
    if (
      !inherits(object, "Seurat") ||
        is.null(phenotype.by) ||
        !phenotype.by %in% colnames(object@meta.data)
    ) {
      log_message(
        "Provide {.arg weights} or a valid Seurat {.arg phenotype.by} column.",
        message_type = "error"
      )
    }
    y <- binary_numeric(object@meta.data[[phenotype.by]], positive = positive)
  }
  log_message(
    "Run FWP scoring with {.val {length(features)}} features",
    verbose = verbose
  )
  bundle <- fwp_score(mat, y = y, weights = weights)
  bundle$parameters <- list(
    assay = input$assay,
    layer = layer,
    features = features,
    nfeatures = nfeatures,
    phenotype.by = phenotype.by,
    positive = positive
  )
  if (inherits(object, "Seurat")) {
    object <- Seurat::AddMetaData(
      object,
      data.frame(
        FWP_Score = bundle$score[colnames(object)],
        row.names = colnames(object)
      )
    )
    if (!identical(score.name, "FWP_Score")) {
      object[[score.name]] <- object[["FWP_Score", drop = TRUE]]
    }
    object@tools[[tool_name]] <- bundle
    object <- suppressWarnings(Seurat::LogSeuratCommand(object))
    return(object)
  }
  bundle
}
