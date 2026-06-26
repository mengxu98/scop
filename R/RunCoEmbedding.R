#' @title Co-embed reference and query cells
#'
#' @description
#' Impute reference expression into query cells with `TransferData`, merge
#' reference and query cells, and compute PCA/UMAP on the shared expression
#' space. The current implementation mirrors the Seurat ATAC integration
#' vignette for ATAC query and RNA reference.
#'
#' @md
#' @inheritParams RunLabelTransfer
#' @param genes_use Genes/features used for expression imputation and co-embedding.
#' If `NULL`, reference variable features are used.
#' @param imputed_assay Assay name for imputed RNA expression in ATAC cells.
#' @param modality_col Metadata column storing RNA/ATAC origin after merge.
#' @param coembed_prefix Prefix for co-embedding reductions.
#' @param npcs Number of PCs to compute.
#' @param umap_dims PC dimensions used for UMAP.
#' @param seed Random seed used for UMAP.
#'
#' @return A merged `Seurat` object containing RNA reference and ATAC query cells.
#' @export
#' @examples
#' data(pbmcmultiome_sub)
#' pbmcmultiome_sub <- standard_scop(
#'   pbmcmultiome_sub,
#'   assay = c("RNA", "peaks"),
#'   linear_reduction_dims = 20
#' )
#' coembed <- RunCoEmbedding(
#'   srt = pbmcmultiome_sub,
#'   reference = pbmcmultiome_sub,
#'   assay = "peaks",
#'   reference_assay = "RNA",
#'   gene_activity_assay = "RNA",
#'   reference_reduction = "RNApca",
#'   reference_dims = 1:10,
#'   dims = 2:10,
#'   umap_dims = 1:10
#' )
#'
#' CellDimPlot(
#'   coembed,
#'   group.by = c("modality", "CellType"),
#'   xlab = "CoEmbedUMAP_1",
#'   ylab = "CoEmbedUMAP_2"
#' )
RunCoEmbedding <- function(
  srt,
  reference,
  assay = NULL,
  reference_assay = NULL,
  reference_reduction = "pca",
  reference_dims = 1:30,
  gene_activity_assay = "ACTIVITY",
  weight_reduction = NULL,
  dims = 2:30,
  genes_use = NULL,
  imputed_assay = "RNA",
  modality_col = "modality",
  coembed_prefix = "CoEmbed",
  npcs = 30,
  umap_dims = 1:30,
  k.weight = 100,
  verbose = TRUE,
  seed = 11
) {
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!inherits(srt[[assay]], "ChromatinAssay")) {
    log_message(
      "{.arg assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} is not a {.cls Seurat}",
      message_type = "error"
    )
  }

  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)
  genes_use <- genes_use %||% SeuratObject::VariableFeatures(reference, assay = reference_assay)
  genes_use <- intersect(genes_use, rownames(reference[[reference_assay]]))
  if (length(genes_use) == 0) {
    log_message(
      "No shared reference features available for ATAC co-embedding",
      message_type = "error"
    )
  }

  srt <- atac_add_activity(
    srt = srt,
    assay = assay,
    gene_activity_assay = gene_activity_assay,
    verbose = verbose
  )

  anchor_features <- intersect(genes_use, rownames(srt[[gene_activity_assay]]))
  if (length(anchor_features) == 0) {
    log_message(
      "No shared features between reference and ATAC gene activity assay",
      message_type = "error"
    )
  }
  weight_reduction <- atac_weight_reduction(
    srt = srt,
    prefix = NULL,
    weight_reduction = weight_reduction,
    verbose = verbose
  )
  anchor_params <- atac_anchor_params(
    reference_cells = ncol(reference),
    query_cells = ncol(srt),
    verbose = verbose
  )

  anchor_reduction <- if (!is.null(reference_reduction)) "pcaproject" else "cca"
  log_message("Finding RNA-to-ATAC transfer anchors...", verbose = verbose)
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = srt,
    reference.assay = reference_assay,
    query.assay = gene_activity_assay,
    reduction = anchor_reduction,
    reference.reduction = reference_reduction,
    dims = reference_dims,
    features = anchor_features,
    k.anchor = anchor_params$k.anchor,
    k.filter = anchor_params$k.filter,
    k.score = anchor_params$k.score,
    verbose = FALSE
  )
  k.weight_use <- atac_k_weight(
    anchors = anchors,
    query_cells = ncol(srt),
    k.weight = k.weight,
    verbose = verbose
  )

  ref_data <- GetAssayData5(
    reference,
    assay = reference_assay,
    layer = "data"
  )
  ref_data <- ref_data[anchor_features, , drop = FALSE]
  reference[[reference_assay]] <- SeuratObject::CreateAssay5Object(data = ref_data)

  log_message("Imputing RNA expression into ATAC cells...", verbose = verbose)
  imputed <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_data,
    weight.reduction = srt[[weight_reduction]],
    dims = dims,
    k.weight = k.weight_use,
    verbose = FALSE
  )
  if (inherits(imputed, "Assay") || inherits(imputed, "StdAssay")) {
    if (imputed_assay %in% names(srt@assays)) {
      SeuratObject::DefaultAssay(srt) <- assay
      srt[[imputed_assay]] <- NULL
    }
    imputed_data <- SeuratObject::GetAssayData(imputed, layer = "data")
    srt[[imputed_assay]] <- SeuratObject::CreateAssay5Object(data = imputed_data)
  } else {
    if (imputed_assay %in% names(srt@assays)) {
      SeuratObject::DefaultAssay(srt) <- assay
      srt[[imputed_assay]] <- NULL
    }
    srt[[imputed_assay]] <- SeuratObject::CreateAssay5Object(data = imputed)
  }

  reference[[modality_col]] <- "RNA"
  srt[[modality_col]] <- "ATAC"
  SeuratObject::DefaultAssay(reference) <- reference_assay
  SeuratObject::DefaultAssay(srt) <- imputed_assay

  coembed <- merge(reference, srt)
  SeuratObject::DefaultAssay(coembed) <- reference_assay
  coembed <- Seurat::ScaleData(
    coembed,
    features = anchor_features,
    verbose = FALSE
  )
  coembed <- Seurat::RunPCA(
    coembed,
    features = anchor_features,
    npcs = npcs,
    reduction.name = paste0(coembed_prefix, "PCA"),
    reduction.key = paste0(coembed_prefix, "PCA_"),
    verbose = FALSE
  )
  coembed <- Seurat::RunUMAP(
    coembed,
    reduction = paste0(coembed_prefix, "PCA"),
    dims = umap_dims,
    reduction.name = paste0(coembed_prefix, "UMAP"),
    reduction.key = paste0(coembed_prefix, "UMAP_"),
    verbose = FALSE,
    seed.use = seed
  )
  coembed@misc[["scop_atac_coembedding"]] <- list(
    query_assay = assay,
    gene_activity_assay = gene_activity_assay,
    imputed_assay = imputed_assay,
    features = anchor_features
  )
  coembed
}
