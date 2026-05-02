#' @title Run BayesSpace spatial clustering
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param q Number of BayesSpace clusters.
#' @param platform Spatial sequencing platform.
#' @param image Name of the Seurat spatial image used to recover spot
#' coordinates when they are not already present in metadata.
#' @param use_reduction Optional Seurat reduction to pass to BayesSpace as PCA.
#' @param dims Dimensions from `use_reduction` to use.
#' @param preprocess Whether to run `BayesSpace::spatialPreprocess()`.
#' @param n.PCs,n.HVGs Parameters passed to `spatialPreprocess()`.
#' @param skip.PCA Whether to skip PCA inside `spatialPreprocess()`.
#' @param spatial_preprocess_params Additional parameters passed to
#' `BayesSpace::spatialPreprocess()`.
#' @param spatial_cluster_params Additional parameters passed to
#' `BayesSpace::spatialCluster()`.
#' @param cluster_colname Metadata column used for BayesSpace clusters.
#' @param init_colname Metadata column used for BayesSpace initial clusters.
#' @param store_sce Whether to store the BayesSpace `SingleCellExperiment`
#' in `srt@tools`.
#'
#' @return A `Seurat` object with BayesSpace clusters in metadata and raw
#' results in `srt@tools[["BayesSpace"]]`.
#' @export
RunBayesSpace <- function(
  srt,
  q,
  assay = NULL,
  platform = c("Visium", "VisiumHD", "ST"),
  image = NULL,
  use_reduction = NULL,
  dims = 1:15,
  preprocess = TRUE,
  n.PCs = 15,
  n.HVGs = 2000,
  skip.PCA = !is.null(use_reduction),
  spatial_preprocess_params = list(),
  spatial_cluster_params = list(),
  cluster_colname = "BayesSpace_cluster",
  init_colname = "BayesSpace_init",
  store_sce = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (missing(q) || length(q) != 1L || !is.numeric(q) || is.na(q) || q < 2) {
    log_message(
      "{.arg q} must be a single number >= 2",
      message_type = "error"
    )
  }
  platform <- match.arg(platform)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  check_r(
    c("BayesSpace", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors"),
    verbose = FALSE
  )

  log_message(
    "Convert {.cls Seurat} to {.cls SingleCellExperiment} for {.pkg BayesSpace}",
    verbose = verbose
  )
  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  sce <- bayesspace_add_spatial_coords(
    srt = srt,
    sce = sce,
    image = image
  )

  use_dimred <- "PCA"
  d_use <- min(length(dims), n.PCs)
  if (!is.null(use_reduction)) {
    if (!use_reduction %in% SeuratObject::Reductions(srt)) {
      log_message(
        "{.arg use_reduction} {.val {use_reduction}} is not present in {.cls Seurat}",
        message_type = "error"
      )
    }
    emb <- SeuratObject::Embeddings(srt, reduction = use_reduction)
    dims_use <- dims[dims <= ncol(emb)]
    if (length(dims_use) == 0L) {
      log_message(
        "No valid {.arg dims} are available in reduction {.val {use_reduction}}",
        message_type = "error"
      )
    }
    emb <- emb[colnames(sce), dims_use, drop = FALSE]
    SingleCellExperiment::reducedDim(sce, use_dimred) <- emb
    d_use <- length(dims_use)
  }

  if (isTRUE(preprocess)) {
    preprocess_args <- c(
      list(
        sce = sce,
        platform = platform,
        n.PCs = n.PCs,
        n.HVGs = n.HVGs,
        skip.PCA = skip.PCA
      ),
      spatial_preprocess_params
    )
    sce <- do.call(BayesSpace::spatialPreprocess, preprocess_args)
  }

  cluster_args <- c(
    list(
      sce = sce,
      q = as.integer(q),
      use.dimred = use_dimred,
      d = d_use,
      platform = platform
    ),
    spatial_cluster_params
  )
  log_message(
    "Run {.pkg BayesSpace} spatial clustering with {.arg q = {q}}",
    verbose = verbose
  )
  sce <- do.call(BayesSpace::spatialCluster, cluster_args)

  cdata <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!"spatial.cluster" %in% colnames(cdata)) {
    log_message(
      "{.pkg BayesSpace} did not return {.val spatial.cluster}",
      message_type = "error"
    )
  }
  cluster_df <- data.frame(
    BayesSpace_cluster = as.character(cdata[["spatial.cluster"]]),
    row.names = colnames(sce),
    stringsAsFactors = FALSE
  )
  colnames(cluster_df) <- cluster_colname
  srt <- Seurat::AddMetaData(srt, metadata = cluster_df)

  if ("cluster.init" %in% colnames(cdata) && !is.null(init_colname)) {
    init_df <- data.frame(
      BayesSpace_init = as.character(cdata[["cluster.init"]]),
      row.names = colnames(sce),
      stringsAsFactors = FALSE
    )
    colnames(init_df) <- init_colname
    srt <- Seurat::AddMetaData(srt, metadata = init_df)
  }

  tool <- list(
    colData = cdata,
    parameters = list(
      assay = assay,
      q = q,
      platform = platform,
      image = image,
      use_reduction = use_reduction,
      dims = dims,
      preprocess = preprocess,
      n.PCs = n.PCs,
      n.HVGs = n.HVGs,
      skip.PCA = skip.PCA,
      spatial_preprocess_params = spatial_preprocess_params,
      spatial_cluster_params = spatial_cluster_params,
      cluster_colname = cluster_colname,
      init_colname = init_colname
    )
  )
  if (isTRUE(store_sce)) {
    tool$sce <- sce
  }
  srt@tools[["BayesSpace"]] <- tool

  log_message(
    "{.pkg BayesSpace} clusters stored in metadata column {.val {cluster_colname}}",
    verbose = verbose
  )
  srt
}

bayesspace_add_spatial_coords <- function(srt, sce, image = NULL) {
  cdata <- as.data.frame(SummarizedExperiment::colData(sce))
  if (all(c("array_row", "array_col") %in% colnames(cdata))) {
    return(sce)
  }

  coords <- bayesspace_get_seurat_coords(srt, image = image)
  if (is.null(coords)) {
    log_message(
      "BayesSpace requires spatial coordinates. Provide a Seurat object with image coordinates or metadata columns {.val array_row}/{.val array_col}.",
      message_type = "error"
    )
  }

  coords <- bayesspace_normalize_coords(coords)
  common <- intersect(colnames(sce), rownames(coords))
  if (length(common) == 0L) {
    log_message(
      "Unable to match spatial coordinates to {.cls SingleCellExperiment} cells",
      message_type = "error"
    )
  }
  cdata[common, colnames(coords)] <- coords[common, , drop = FALSE]
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cdata)
  sce
}

bayesspace_get_seurat_coords <- function(srt, image = NULL) {
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  if (length(images) == 0L) {
    return(NULL)
  }
  image <- image %||% images[1]
  if (!image %in% images) {
    log_message(
      "{.arg image} {.val {image}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  coords <- tryCatch(
    {
      as.data.frame(SeuratObject::GetTissueCoordinates(srt[[image]]))
    },
    error = function(e) {
      tryCatch(
        as.data.frame(slot(srt@images[[image]], "coordinates")),
        error = function(e2) NULL
      )
    }
  )
  coords
}

bayesspace_normalize_coords <- function(coords) {
  coords <- as.data.frame(coords)
  row_col <- bayesspace_pick_col(coords, c("array_row", "row", "arrayrow"))
  col_col <- bayesspace_pick_col(coords, c("array_col", "col", "arraycol"))
  pxl_row_col <- bayesspace_pick_col(
    coords,
    c("pxl_row_in_fullres", "imagerow", "image_row", "y")
  )
  pxl_col_col <- bayesspace_pick_col(
    coords,
    c("pxl_col_in_fullres", "imagecol", "image_col", "x")
  )
  if (is.null(row_col) || is.null(col_col)) {
    log_message(
      "BayesSpace coordinate table must contain array row/column coordinates",
      message_type = "error"
    )
  }
  out <- data.frame(
    array_row = coords[[row_col]],
    array_col = coords[[col_col]],
    row = coords[[row_col]],
    col = coords[[col_col]],
    row.names = rownames(coords)
  )
  if (!is.null(pxl_row_col)) {
    out$pxl_row_in_fullres <- coords[[pxl_row_col]]
  }
  if (!is.null(pxl_col_col)) {
    out$pxl_col_in_fullres <- coords[[pxl_col_col]]
  }
  out
}

bayesspace_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}
