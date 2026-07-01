#' @title Run BayesSpace spatial clustering
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param q Number of BayesSpace clusters.
#' @param platform Spatial sequencing platform.
#' @param image Name of the Seurat spatial image used to recover spot
#' coordinates when they are not already present in metadata. For regular
#' Visium data with only pixel `x`/`y` coordinates, BayesSpace array
#' coordinates are inferred from the spatial grid.
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
#'
#' @export
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial$BayesSpace_cluster <- factor(
#'   paste0("domain_", (seq_len(ncol(spatial)) - 1) %% 3 + 1)
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "BayesSpace_cluster",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   requireNamespace("BayesSpace", quietly = TRUE)
#' ) {
#' spatial <- RunBayesSpace(
#'   spatial,
#'   q = 3,
#'   n.PCs = 5,
#'   n.HVGs = 200,
#'   store_sce = FALSE,
#'   spatial_cluster_params = list(
#'     nrep = 200,
#'     burn.in = 50,
#'     thin = 10,
#'     save.chain = FALSE
#'   )
#' )
#' table(spatial$BayesSpace_cluster)
#' }
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
    image = image,
    platform = platform
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

bayesspace_add_spatial_coords <- function(
  srt,
  sce,
  image = NULL,
  platform = "Visium"
) {
  cdata <- as.data.frame(SummarizedExperiment::colData(sce))
  if (all(c("array_row", "array_col") %in% colnames(cdata))) {
    return(sce)
  }

  coords <- bayesspace_get_seurat_coords(srt, image = image)
  if (is.null(coords)) {
    coords <- cdata
  }

  coords <- bayesspace_normalize_coords(coords, platform = platform)
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

bayesspace_normalize_coords <- function(coords, platform = "Visium") {
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
    coords <- bayesspace_infer_array_coords_from_pixels(
      coords = coords,
      platform = platform,
      pxl_row_col = pxl_row_col,
      pxl_col_col = pxl_col_col
    )
    row_col <- "array_row"
    col_col <- "array_col"
    pxl_row_col <- bayesspace_pick_col(
      coords,
      c("pxl_row_in_fullres", "imagerow", "image_row", "y")
    )
    pxl_col_col <- bayesspace_pick_col(
      coords,
      c("pxl_col_in_fullres", "imagecol", "image_col", "x")
    )
  }
  if (is.null(row_col) || is.null(col_col)) {
    log_message(
      "BayesSpace coordinate table must contain array row/column coordinates or regular pixel coordinates",
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

bayesspace_infer_array_coords_from_pixels <- function(
  coords,
  platform = "Visium",
  pxl_row_col = NULL,
  pxl_col_col = NULL
) {
  if (is.null(pxl_row_col) || is.null(pxl_col_col)) {
    log_message(
      "BayesSpace coordinate table must contain array row/column coordinates or pixel coordinate columns such as {.val x}/{.val y}",
      message_type = "error"
    )
  }
  pxl_row <- suppressWarnings(as.numeric(coords[[pxl_row_col]]))
  pxl_col <- suppressWarnings(as.numeric(coords[[pxl_col_col]]))
  keep <- is.finite(pxl_row) & is.finite(pxl_col)
  if (sum(keep) < 3L) {
    log_message(
      "At least three finite spatial coordinate pairs are required to infer BayesSpace array coordinates",
      message_type = "error"
    )
  }
  platform <- match.arg(platform, c("Visium", "VisiumHD", "ST"))
  if (platform == "Visium") {
    array_row <- bayesspace_infer_axis_index(pxl_row)
    same_row_diffs <- unlist(tapply(pxl_col, array_row, function(x) {
      diff(sort(unique(x[is.finite(x)])))
    }))
    col_pair_step <- bayesspace_nearest_regular_step(same_row_diffs)
    if (!is.finite(col_pair_step) || col_pair_step <= 0) {
      log_message(
        "Unable to infer Visium array columns from pixel coordinates",
        message_type = "error"
      )
    }
    array_col <- round((pxl_col - min(pxl_col, na.rm = TRUE)) /
      (col_pair_step / 2))
  } else {
    array_row <- bayesspace_infer_axis_index(pxl_row)
    array_col <- bayesspace_infer_axis_index(pxl_col)
  }
  if (any(!is.finite(array_row[keep])) || any(!is.finite(array_col[keep]))) {
    log_message(
      "Unable to infer BayesSpace array coordinates from pixel coordinates",
      message_type = "error"
    )
  }
  out <- data.frame(
    array_row = as.integer(array_row),
    array_col = as.integer(array_col),
    row = as.integer(array_row),
    col = as.integer(array_col),
    pxl_row_in_fullres = pxl_row,
    pxl_col_in_fullres = pxl_col,
    row.names = rownames(coords),
    stringsAsFactors = FALSE
  )
  out
}

bayesspace_infer_axis_index <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  step <- bayesspace_regular_step(diff(sort(unique(values[is.finite(values)]))))
  if (!is.finite(step) || step <= 0) {
    log_message(
      "Unable to infer regular spatial grid spacing from pixel coordinates",
      message_type = "error"
    )
  }
  round((values - min(values, na.rm = TRUE)) / step)
}

bayesspace_regular_step <- function(diffs) {
  diffs <- as.numeric(diffs)
  diffs <- diffs[is.finite(diffs) & diffs > 0]
  if (length(diffs) == 0L) {
    return(NA_real_)
  }
  baseline <- stats::median(diffs, na.rm = TRUE)
  large <- diffs[diffs > baseline * 3]
  if (length(large) == 0L) {
    large <- diffs[diffs > max(1, baseline)]
  }
  if (length(large) == 0L) {
    large <- diffs
  }
  stats::median(large, na.rm = TRUE)
}

bayesspace_nearest_regular_step <- function(diffs) {
  diffs <- as.numeric(diffs)
  diffs <- diffs[is.finite(diffs) & diffs > 0]
  if (length(diffs) == 0L) {
    return(NA_real_)
  }
  low <- diffs[diffs <= stats::quantile(diffs, 0.25, na.rm = TRUE)]
  if (length(low) == 0L) {
    low <- diffs
  }
  stats::median(low, na.rm = TRUE)
}

bayesspace_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}
