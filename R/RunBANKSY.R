#' @title Run BANKSY spatial clustering
#'
#' @description
#' Build neighborhood-augmented BANKSY features from a spatial `Seurat` object
#' and store spatial domain or microenvironment clusters in metadata.
#'
#' @md
#' @inheritParams RunBayesSpace
#' @param assay Expression assay used by BANKSY. When `NULL`, the assay linked
#' to the selected image is used first, followed by `Spatial`, `RNA`, and the
#' object's default assay.
#' @param layer Assay layer used as BANKSY input.
#' @param image Spatial image used to recover coordinates. A single available
#' image is selected automatically; choose one explicitly when several exist,
#' or provide a named sample-to-image map with `sample.by`.
#' @param features Optional features to use. If `NULL`, all assay features are
#' used after zero-count filtering.
#' @param coord.cols Metadata coordinate columns used when no image coordinate
#' source is available. The default is resolved from common coordinate names.
#' @param lambda BANKSY spatial weighting parameter.
#' @param k_geom Unitless number of spatial neighbors used by BANKSY.
#' @param M Highest azimuthal Fourier harmonic passed to BANKSY.
#' @param npcs Number of principal components to compute.
#' @param use_agf Whether to use azimuthal Gabor filters.
#' @param algo Clustering algorithm passed to `Banksy::clusterBanksy()`.
#' @param k_neighbors Number of neighbors for graph clustering.
#' @param resolution Graph clustering resolution.
#' @param group Optional metadata column passed to BANKSY for multi-sample
#' scaling. It is an algorithm parameter and does not split sample fits.
#' @param sample.by Optional metadata column for independent per-sample fits.
#' Combined labels are prefixed with the sample name; domains are not aligned
#' across samples.
#' @param seed Optional seed for PCA and clustering.
#' @param compute_banksy_params Additional parameters passed to
#' `Banksy::computeBanksy()`.
#' @param run_pca_params Additional parameters passed to
#' `Banksy::runBanksyPCA()`.
#' @param cluster_banksy_params Additional parameters passed to
#' `Banksy::clusterBanksy()`.
#' @param cluster_source Optional BANKSY `colData` column to copy. If `NULL`,
#' the first cluster name reported by `Banksy::clusterNames()` is used when
#' available.
#' @param cluster_colname Metadata column used for BANKSY clusters.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param store_results Whether to store detailed BANKSY results in
#' `object@tools`; cluster assignments are still written to metadata when
#' `FALSE`; a lightweight index retains the result column and analyzed cells.
#' @param coordinate_space Coordinate space used for BANKSY spatial input.
#' The default is raw acquisition coordinates, so geometry and distance
#' weighting use raw coordinate units. Use `"legacy_display"` explicitly to
#' reproduce the display-scaled coordinates used before scop 0.9.0.
#'
#' @return A `Seurat` object with BANKSY clusters in metadata. When
#' `store_results = TRUE`, detailed results are stored in
#' `srt@tools[[tool_name]]`; when `FALSE`, clusters remain in metadata and only
#' their lightweight retrieval index is retained in the tool entry.
#' @seealso [GetSpatialResult()]
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' keep_spots <- unique(round(seq(1, ncol(visium_human_pancreas_sub), length.out = 400)))
#' spatial <- visium_human_pancreas_sub[, keep_spots]
#' \dontrun{
#' spatial <- RunBANKSY(
#'   spatial,
#'   layer = "counts",
#'   features = rownames(spatial)[1:200],
#'   resolution = 0.6,
#'   verbose = FALSE
#' )
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "BANKSY_cluster",
#'   pt.size = 1.5
#' )
#' }
RunBANKSY <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  lambda = 0.2,
  k_geom = 15,
  M = 1,
  npcs = 20,
  use_agf = FALSE,
  algo = "leiden",
  k_neighbors = 50,
  resolution = 0.6,
  group = NULL,
  seed = 1,
  compute_banksy_params = list(),
  run_pca_params = list(),
  cluster_banksy_params = list(),
  cluster_source = NULL,
  cluster_colname = "BANKSY_cluster",
  tool_name = "BANKSY",
  store_results = TRUE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  srt = NULL,
  sample.by = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  validate_scalar_string(cluster_colname, "cluster_colname", require_character = FALSE)
  validate_scalar_string(tool_name, "tool_name", require_character = FALSE)
  validate_scalar_flag(store_results, "store_results")
  validate_named_param_list(compute_banksy_params, "compute_banksy_params", require_list = TRUE)
  validate_named_param_list(run_pca_params, "run_pca_params", require_list = TRUE)
  validate_named_param_list(cluster_banksy_params, "cluster_banksy_params", require_list = TRUE)
  coordinate_space <- match.arg(coordinate_space)

  if (!is.null(sample.by)) {
    return(banksy_run_by_sample(
      srt = srt,
      assay = assay,
      layer = layer,
      features = features,
      image = image,
      coord.cols = coord.cols,
      lambda = lambda,
      k_geom = k_geom,
      M = M,
      npcs = npcs,
      use_agf = use_agf,
      algo = algo,
      k_neighbors = k_neighbors,
      resolution = resolution,
      group = group,
      sample.by = sample.by,
      seed = seed,
      compute_banksy_params = compute_banksy_params,
      run_pca_params = run_pca_params,
      cluster_banksy_params = cluster_banksy_params,
      cluster_source = cluster_source,
      cluster_colname = cluster_colname,
      tool_name = tool_name,
      store_results = store_results,
      verbose = verbose,
      coordinate_space = coordinate_space
    ))
  }

  image_use <- spatial_image_resolve(
    srt = srt,
    image = image,
    image_policy = "strict"
  )$image
  assay_info <- banksy_resolve_assay(srt, assay = assay, image = image_use)
  assay <- assay_info$assay
  banksy_require_layer(srt, assay = assay, layer = layer)
  features_use <- features %||% rownames(srt[[assay]])
  features_use <- intersect(features_use, rownames(srt[[assay]]))
  if (length(features_use) == 0L) {
    log_message(
      "No features are available for {.fn RunBANKSY}",
      message_type = "error"
    )
  }
  expr <- banksy_get_matrix(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features_use
  )
  keep_features <- Matrix::rowSums(expr) > 0
  keep_spots <- Matrix::colSums(expr) > 0
  expr <- expr[keep_features, keep_spots, drop = FALSE]
  if (nrow(expr) == 0L || ncol(expr) == 0L) {
    log_message(
      "No non-zero features or spots remain for {.fn RunBANKSY}",
      message_type = "error"
    )
  }
  coords <- resolve_spatial_spot_coords(
    srt = srt,
    spot_ids = colnames(expr),
    image = image_use,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )
  coordinate_source <- attr(coords, "spatial_source", exact = TRUE)
  coord_cols_use <- coordinate_source$coord.cols %||% coord.cols

  if (is.null(assay_info$requested)) {
    log_message(
      "Using assay: {.val {assay}} ({assay_info$source})",
      level = 2,
      timestamp = FALSE,
      verbose = verbose
    )
  }
  if (is.null(image) && !is.null(image_use)) {
    log_message(
      "Using image: {.val {image_use}} (auto-selected)",
      level = 2,
      timestamp = FALSE,
      verbose = verbose
    )
  }
  coordinate_label <- if (length(coordinate_source$coord.cols) > 0L) {
    paste(coordinate_source$coord.cols, collapse = ", ")
  } else {
    paste(coord_cols_use, collapse = ", ")
  }
  log_message(
    "Using coordinates: {.val {coordinate_label}} ({coordinate_space} space)",
    level = 2,
    timestamp = FALSE,
    verbose = verbose
  )
  coldata <- banksy_coldata(
    srt = srt,
    spot_ids = colnames(expr),
    coords = coords,
    group = group
  )

  log_message(
    "Run {.pkg Banksy} with {.val {nrow(expr)}} features and {.val {ncol(expr)}} spatial spots",
    verbose = verbose
  )
  backend <- banksy_run_backend(
    expr = expr,
    coords = coords,
    coldata = coldata,
    assay_name = "scop_input",
    lambda = lambda,
    k_geom = k_geom,
    M = M,
    npcs = npcs,
    use_agf = use_agf,
    algo = algo,
    k_neighbors = k_neighbors,
    resolution = resolution,
    group = group,
    seed = seed,
    compute_banksy_params = compute_banksy_params,
    run_pca_params = run_pca_params,
    cluster_banksy_params = cluster_banksy_params
  )

  cluster_source <- banksy_resolve_cluster_source(
    se = backend$se,
    before_cols = backend$before_cols,
    cluster_source = cluster_source
  )
  cdata <- as.data.frame(SummarizedExperiment::colData(backend$se))
  clusters <- as.character(cdata[colnames(expr), cluster_source, drop = TRUE])
  cluster_df <- data.frame(
    BANKSY_cluster = clusters,
    row.names = colnames(expr),
    stringsAsFactors = FALSE
  )
  colnames(cluster_df) <- cluster_colname
  # Replace the whole metadata column so a rerun cannot retain labels for
  # cells excluded by this run (for example, zero-count spots).
  metadata_clusters <- stats::setNames(rep(NA_character_, ncol(srt)), colnames(srt))
  metadata_clusters[rownames(cluster_df)] <- cluster_df[[cluster_colname]]
  srt <- Seurat::AddMetaData(srt, metadata = metadata_clusters, col.name = cluster_colname)
  domain_summary <- spatial_domain_summary(cluster_df[[cluster_colname]])
  n_spots <- nrow(cluster_df)
  n_domains <- nrow(domain_summary)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      clusters = cluster_df,
      cluster_source = cluster_source,
      colData = cdata,
      coords = coords,
      features = rownames(expr),
      se = backend$se,
      summary = list(
        n_spots = n_spots,
        domains = domain_summary
      ),
      parameters = list(
        assay = assay,
        layer = layer,
        image = image_use,
        coord.cols = coord_cols_use,
        coordinate_space = coordinate_space,
        lambda = lambda,
        k_geom = k_geom,
        M = M,
        npcs = npcs,
        use_agf = use_agf,
        algo = algo,
        k_neighbors = k_neighbors,
        resolution = resolution,
        group = group,
        seed = seed,
        compute_banksy_params = compute_banksy_params,
        run_pca_params = run_pca_params,
        cluster_banksy_params = cluster_banksy_params,
        cluster_source = cluster_source,
        cluster_colname = cluster_colname,
        tool_name = tool_name
      )
    )
    srt@tools[[tool_name]] <- spatial_tag_coordinate_contract(srt@tools[[tool_name]])
  } else {
    srt@tools[[tool_name]] <- list()
  }

  srt@tools[[tool_name]]$result_index <- list(
    method = "BANKSY", cluster_colname = cluster_colname, cells = colnames(expr)
  )

  image_use <- coordinate_source$image
  if (length(image_use) != 1L || is.na(image_use) || !nzchar(image_use)) {
    image_use <- NULL
  }
  has_image <- !is.null(image_use)
  display_scale <- if (has_image && isTRUE(thisutils::get_verbose(verbose))) {
    spatial_run_receipt_display_scale(srt, image_use)
  } else {
    NULL
  }
  plot_args <- paste0(
    "group.by = ",
    spatial_run_receipt_quote(cluster_colname, "cluster_colname")
  )
  if (has_image && !is.null(display_scale)) {
    plot_args <- c(plot_args, paste0("image = ", spatial_run_receipt_quote(image_use, "image")))
    if (identical(display_scale, "hires")) {
      plot_args <- c(
        plot_args,
        paste0("image.scale = ", spatial_run_receipt_quote(display_scale, "image.scale"))
      )
    }
  } else if (!has_image) {
    plot_args <- c(
      plot_args,
      paste0(
        "coord.cols = ",
        deparse1(unname(as.character(coord_cols_use)), width.cutoff = 500L)
      )
    )
  }
  plot_call <- if (!has_image || !is.null(display_scale)) {
    paste0("SpatialSpotPlot(<returned_object>, ", paste(plot_args, collapse = ", "), ")")
  } else {
    NULL
  }
  inspect_call <- if (has_image && is.null(display_scale)) {
    paste0(
      "GetSpatialResult(<returned_object>, ",
      spatial_run_receipt_quote(tool_name, "tool_name"),
      ")"
    )
  } else {
    NULL
  }
  saved <- if (isTRUE(store_results)) {
    paste0(
      "metadata column {.var ", cluster_colname,
      "} and returned object's {.code @tools} entry {.var ", tool_name,
      "}; use GetSpatialResult(<returned_object>, ",
      spatial_run_receipt_quote(tool_name, "tool_name"),
      ") for clusters, parameters, and summary"
    )
  } else {
    paste0(
      "metadata column {.var ", cluster_colname,
      "}; detailed results were not stored, and GetSpatialResult() returns metadata clusters with domain counts"
    )
  }
  scope <- paste0(
    "assay {.val ", assay, "}, ",
    if (has_image) paste0("image {.val ", image_use, "}, ") else "",
    "coordinates {.val ", coordinate_label, "} ({coordinate_space} space)"
  )
  spatial_run_receipt(
    done = paste0(
      "{.pkg BANKSY} completed ({.val ", n_domains, "} domains detected; ",
      "{.val ", n_spots, "} spots)"
    ),
    scope = scope,
    saved = saved,
    plot = plot_call,
    inspect = inspect_call,
    verbose = verbose,
    .envir = environment()
  )
  srt
}

banksy_resolve_assay <- function(srt, assay = NULL, image = NULL) {
  available <- SeuratObject::Assays(srt)
  requested <- assay
  if (!is.null(assay)) {
    if (!is.character(assay) || length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
      log_message("{.arg assay} must be one non-empty assay name", message_type = "error")
    }
    if (!assay %in% available) {
      log_message(
        "Assay {.val {assay}} is not present; available assays: {.val {available}}",
        message_type = "error"
      )
    }
    return(list(assay = assay, requested = requested, source = "explicitly selected"))
  }

  if (!is.null(image)) {
    if (!image %in% names(srt@images)) {
      log_message(
        "Image {.val {image}} is not present in {.arg srt}; available images: {.val {names(srt@images)}}",
        message_type = "error"
      )
    }
    image_assay <- tryCatch(
      methods::slot(srt[[image]], "assay"),
      error = function(e) NULL
    )
    if (is.character(image_assay) && length(image_assay) == 1L && !is.na(image_assay) && nzchar(image_assay)) {
      if (!image_assay %in% available) {
        log_message(
          "Image {.val {image}} refers to assay {.val {image_assay}}, which is not present in {.arg srt}",
          message_type = "error"
        )
      }
      return(list(
        assay = image_assay,
        requested = requested,
        source = paste0("from image ", image)
      ))
    }
  }
  if ("Spatial" %in% available) {
    return(list(assay = "Spatial", requested = requested, source = "Spatial assay"))
  }
  if ("RNA" %in% available) {
    return(list(assay = "RNA", requested = requested, source = "RNA assay"))
  }
  default_assay <- SeuratObject::DefaultAssay(srt)
  if (!default_assay %in% available) {
    log_message("The default assay is not present in {.arg srt}", message_type = "error")
  }
  list(assay = default_assay, requested = requested, source = "DefaultAssay")
}

banksy_require_layer <- function(srt, assay, layer) {
  if (!is.character(layer) || length(layer) != 1L || is.na(layer) || !nzchar(layer)) {
    log_message("{.arg layer} must be one non-empty layer name", message_type = "error")
  }
  assay_object <- srt[[assay]]
  available <- tryCatch(SeuratObject::Layers(assay_object), error = function(e) character())
  if (length(available) == 0L) {
    slots <- intersect(c("counts", "data", "scale.data"), methods::slotNames(assay_object))
    available <- slots[vapply(slots, function(slot_name) {
      value <- methods::slot(assay_object, slot_name)
      !is.null(value) && length(dim(value)) == 2L && all(dim(value) > 0L)
    }, logical(1))]
  }
  # A layer may be split into sample-specific layers such as "data.sample1".
  exact_match <- available == layer
  prefix_match <- startsWith(available, paste0(layer, "."))
  matches <- exact_match | prefix_match

  if (!any(matches)) {
    available_label <- if (length(available) == 0L) "none" else paste(available, collapse = ", ")
    log_message(
      paste0(
        "Layer ", shQuote(layer), " is not present in assay ", shQuote(assay),
        "; available layers: ", available_label,
        ". Supply an existing {.arg layer}; RunBANKSY does not normalize data or switch layers automatically."
      ),
      message_type = "error"
    )
  }
  invisible(available)
}

banksy_run_by_sample <- function(
  srt,
  assay,
  layer,
  features,
  image,
  coord.cols,
  lambda,
  k_geom,
  M,
  npcs,
  use_agf,
  algo,
  k_neighbors,
  resolution,
  group,
  sample.by,
  seed,
  compute_banksy_params,
  run_pca_params,
  cluster_banksy_params,
  cluster_source,
  cluster_colname,
  tool_name,
  store_results,
  verbose,
  coordinate_space
) {
  if (!is.character(sample.by) || length(sample.by) != 1L || is.na(sample.by) ||
    !nzchar(sample.by) || !sample.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg sample.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  sample_values <- as.character(srt@meta.data[[sample.by]])
  if (length(sample_values) != ncol(srt) || anyNA(sample_values) || any(!nzchar(sample_values))) {
    log_message("{.arg sample.by} must identify every cell or spot", message_type = "error")
  }
  names(sample_values) <- rownames(srt@meta.data)
  samples <- unique(sample_values)
  image_map <- spatial_resolve_sample_images(srt, sample.by, image = image)

  combined <- stats::setNames(rep(NA_character_, ncol(srt)), colnames(srt))
  analyzed_samples <- stats::setNames(character(), character())
  sample_results <- stats::setNames(vector("list", length(samples)), samples)
  sample_summaries <- stats::setNames(vector("list", length(samples)), samples)
  assay_by_sample <- stats::setNames(character(length(samples)), samples)
  coordinate_sources <- stats::setNames(vector("list", length(samples)), samples)

  for (sample_name in samples) {
    cells <- colnames(srt)[sample_values[colnames(srt)] == sample_name]
    sample_srt <- srt[, cells]
    image_use <- image_map[[sample_name]]
    if (is.na(image_use)) image_use <- NULL
    sample_output <- tryCatch({
      assay_info <- banksy_resolve_assay(sample_srt, assay = assay, image = image_use)
      coords <- resolve_spatial_spot_coords(
        srt = sample_srt,
        spot_ids = cells,
        image = image_use,
        coord.cols = coord.cols,
        coordinate_space = coordinate_space
      )
      result <- RunBANKSY(
        object = sample_srt,
        assay = assay,
        layer = layer,
        features = features,
        image = image_use,
        coord.cols = coord.cols,
        lambda = lambda,
        k_geom = k_geom,
        M = M,
        npcs = npcs,
        use_agf = use_agf,
        algo = algo,
        k_neighbors = k_neighbors,
        resolution = resolution,
        group = group,
        seed = seed,
        compute_banksy_params = compute_banksy_params,
        run_pca_params = run_pca_params,
        cluster_banksy_params = cluster_banksy_params,
        cluster_source = cluster_source,
        cluster_colname = cluster_colname,
        tool_name = tool_name,
        store_results = store_results,
        verbose = FALSE,
        coordinate_space = coordinate_space
      )
      list(
        result = result,
        assay = assay_info$assay,
        coordinate_source = attr(coords, "spatial_source", exact = TRUE)
      )
    },
      error = function(e) {
        log_message(
          paste0("BANKSY failed for sample ", shQuote(sample_name), ": ", conditionMessage(e)),
          message_type = "error"
        )
      }
    )
    sample_result <- sample_output$result
    assay_by_sample[[sample_name]] <- sample_output$assay
    coordinate_sources[[sample_name]] <- sample_output$coordinate_source
    analyzed_cells <- sample_result@tools[[tool_name]]$result_index$cells
    analyzed_samples[analyzed_cells] <- sample_name
    sample_clusters <- as.character(sample_result@meta.data[cells, cluster_colname, drop = TRUE])
    names(sample_clusters) <- cells
    assigned <- !is.na(sample_clusters) & nzchar(sample_clusters)
    combined[cells[assigned]] <- paste(sample_name, sample_clusters[assigned], sep = "_")
    sample_summaries[[sample_name]] <- list(
      n_spots = sum(assigned),
      domains = spatial_domain_summary(sample_clusters)
    )
    if (isTRUE(store_results)) {
      sample_results[[sample_name]] <- sample_result@tools[[tool_name]]
    }
  }

  combined_df <- data.frame(
    stats::setNames(list(unname(combined)), cluster_colname),
    row.names = names(combined),
    stringsAsFactors = FALSE
  )
  srt <- Seurat::AddMetaData(srt, metadata = combined_df)
  domains <- spatial_domain_summary(combined)
  n_spots <- sum(!is.na(combined) & nzchar(combined))
  n_domains <- nrow(domains)

  if (isTRUE(store_results)) {
    parameters <- list(
      assay = assay,
      assay_by_sample = assay_by_sample,
      layer = layer,
      image = image_map,
      coord.cols = lapply(coordinate_sources, function(source) source$coord.cols),
      coordinate_space = coordinate_space,
      lambda = lambda,
      k_geom = k_geom,
      M = M,
      npcs = npcs,
      use_agf = use_agf,
      algo = algo,
      k_neighbors = k_neighbors,
      resolution = resolution,
      group = group,
      sample.by = sample.by,
      seed = seed,
      compute_banksy_params = compute_banksy_params,
      run_pca_params = run_pca_params,
      cluster_banksy_params = cluster_banksy_params,
      cluster_source = cluster_source,
      cluster_colname = cluster_colname,
      tool_name = tool_name
    )
    srt@tools[[tool_name]] <- spatial_tag_coordinate_contract(list(
      clusters = combined_df,
      per_sample = sample_results,
      summary = list(n_spots = n_spots, domains = domains),
      parameters = parameters
    ))
  } else {
    srt@tools[[tool_name]] <- list()
  }

  srt@tools[[tool_name]]$result_index <- list(
    method = "BANKSY", cluster_colname = cluster_colname,
    samples = analyzed_samples
  )

  if (isTRUE(thisutils::get_verbose(verbose))) {
    for (sample_name in samples) {
      sample_summary <- sample_summaries[[sample_name]]
      log_message(
        paste0(
          sample_name, ": ", nrow(sample_summary$domains), " domains (",
          sample_summary$n_spots, " spots)"
        ),
        message_type = "success",
        verbose = TRUE
      )
    }
  }

  assay_scope <- paste(
    paste0(samples, "=", vapply(assay_by_sample, function(x) {
      spatial_run_receipt_quote(x, "assay")
    }, character(1))),
    collapse = ", "
  )
  image_labels <- ifelse(is.na(image_map), "metadata coordinates", unname(image_map))
  image_scope <- paste(paste0(samples, "=", image_labels), collapse = ", ")
  coordinate_labels <- vapply(coordinate_sources, function(source) {
    paste(source$coord.cols, collapse = ", ")
  }, character(1))
  coordinate_scope <- paste(
    paste0(samples, "=", coordinate_labels),
    collapse = ", "
  )
  scope <- c(
    paste0(
      "sample.by ", spatial_run_receipt_quote(sample.by, "sample.by"), ": ",
      paste(samples, collapse = ", ")
    ),
    paste0("assay by sample: ", assay_scope),
    paste0("image by sample: ", image_scope),
    paste0("coordinates by sample: ", coordinate_scope, " (", coordinate_space, " space)")
  )
  saved <- if (isTRUE(store_results)) {
    paste0(
      "metadata column {.var ", cluster_colname,
      "} and per-sample tool results; for example, call GetSpatialResult(<returned_object>, ",
      spatial_run_receipt_quote(tool_name, "tool_name"),
      ", sample = ", spatial_run_receipt_quote(samples[[1L]], "sample"), ")"
    )
  } else {
    paste0(
      "metadata column {.var ", cluster_colname,
      "}; detailed results were not stored"
    )
  }

  images_use <- unique(unname(image_map[!is.na(image_map)]))
  one_shared_image <- length(images_use) == 1L && all(!is.na(image_map))
  display_scale <- if (one_shared_image && isTRUE(thisutils::get_verbose(verbose))) {
    spatial_run_receipt_display_scale(srt, images_use[[1L]])
  } else {
    NULL
  }
  plot_call <- if (one_shared_image && !is.null(display_scale)) {
    plot_args <- c(
      paste0("group.by = ", spatial_run_receipt_quote(cluster_colname, "cluster_colname")),
      paste0("image = ", spatial_run_receipt_quote(images_use[[1L]], "image"))
    )
    if (identical(display_scale, "hires")) {
      plot_args <- c(
        plot_args,
        paste0("image.scale = ", spatial_run_receipt_quote(display_scale, "image.scale"))
      )
    }
    paste0("SpatialSpotPlot(<returned_object>, ", paste(plot_args, collapse = ", "), ")")
  } else {
    NULL
  }
  inspect_call <- if (is.null(plot_call)) {
    paste0(
      "GetSpatialResult(<returned_object>, ",
      spatial_run_receipt_quote(tool_name, "tool_name"),
      ")"
    )
  } else {
    NULL
  }
  spatial_run_receipt(
    done = paste0(
      "{.pkg BANKSY} completed across {.val ", length(samples), "} samples ({.val ",
      n_domains, "} domains; {.val ", n_spots, "} spots)"
    ),
    scope = scope,
    saved = saved,
    plot = plot_call,
    inspect = inspect_call,
    verbose = verbose,
    .envir = environment()
  )
  srt
}

banksy_run_backend <- function(
  expr,
  coords,
  coldata,
  assay_name,
  lambda,
  k_geom,
  M,
  npcs,
  use_agf,
  algo,
  k_neighbors,
  resolution,
  group,
  seed,
  compute_banksy_params,
  run_pca_params,
  cluster_banksy_params
) {
  check_r(
    c("Banksy", "SpatialExperiment", "SummarizedExperiment", "S4Vectors"),
    verbose = FALSE
  )
  spatial_experiment <- get_namespace_fun("SpatialExperiment", "SpatialExperiment")
  se <- spatial_experiment(
    assays = stats::setNames(list(expr), assay_name),
    colData = S4Vectors::DataFrame(coldata),
    spatialCoords = as.matrix(coords)
  )
  before_cols <- colnames(as.data.frame(SummarizedExperiment::colData(se)))
  compute_fun <- get_namespace_fun("Banksy", "computeBanksy")
  pca_fun <- get_namespace_fun("Banksy", "runBanksyPCA")
  cluster_fun <- get_namespace_fun("Banksy", "clusterBanksy")

  compute_args <- c(
    list(
      assay_name = assay_name,
      coord_names = c("x", "y"),
      compute_agf = use_agf,
      M = M,
      k_geom = k_geom
    ),
    compute_banksy_params
  )
  se <- banksy_do_call(compute_fun, se, compute_args)

  pca_args <- c(
    list(
      assay_name = assay_name,
      M = M,
      lambda = lambda,
      npcs = npcs,
      use_agf = use_agf,
      group = group,
      seed = seed
    ),
    run_pca_params
  )
  se <- banksy_do_call(pca_fun, se, pca_args)

  cluster_args <- c(
    list(
      assay_name = assay_name,
      M = M,
      lambda = lambda,
      use_agf = use_agf,
      npcs = npcs,
      algo = algo,
      k_neighbors = k_neighbors,
      resolution = resolution,
      group = group,
      seed = seed
    ),
    cluster_banksy_params
  )
  se <- banksy_do_call(cluster_fun, se, cluster_args)
  list(se = se, before_cols = before_cols)
}

banksy_get_matrix <- function(srt, assay, layer, features) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  mat <- mat[features, , drop = FALSE]
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(
      if (is.data.frame(mat)) as.matrix(mat) else mat,
      sparse = TRUE
    )
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  mat@x[!is.finite(mat@x)] <- 0
  Matrix::drop0(mat)
}

banksy_coldata <- function(srt, spot_ids, coords, group = NULL) {
  coldata <- data.frame(
    x = coords$x,
    y = coords$y,
    row.names = spot_ids,
    stringsAsFactors = FALSE
  )
  if (!is.null(group)) {
    if (length(group) != 1L || !group %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg group} must be a single metadata column in {.arg srt}",
        message_type = "error"
      )
    }
    coldata[[group]] <- srt@meta.data[spot_ids, group, drop = TRUE]
  }
  coldata
}

banksy_resolve_cluster_source <- function(se, before_cols, cluster_source = NULL) {
  cdata <- as.data.frame(SummarizedExperiment::colData(se))
  if (!is.null(cluster_source)) {
    if (length(cluster_source) != 1L || !cluster_source %in% colnames(cdata)) {
      log_message(
        "{.arg cluster_source} must be a single BANKSY colData column",
        message_type = "error"
      )
    }
    return(cluster_source)
  }
  cluster_names <- tryCatch(
    get_namespace_fun("Banksy", "clusterNames")(se),
    error = function(e) character()
  )
  cluster_names <- cluster_names[cluster_names %in% colnames(cdata)]
  if (length(cluster_names) > 0L) {
    return(cluster_names[1L])
  }
  new_cols <- setdiff(colnames(cdata), before_cols)
  candidate <- new_cols[vapply(cdata[, new_cols, drop = FALSE], function(x) {
    is.factor(x) || is.character(x) || is.numeric(x)
  }, logical(1))]
  if (length(candidate) == 0L) {
    log_message(
      "{.pkg Banksy} did not add a detectable cluster column",
      message_type = "error"
    )
  }
  candidate[1L]
}

banksy_do_call <- function(fun, se, args) {
  args <- args[!vapply(args, is.null, logical(1))]
  fmls <- names(formals(fun))
  if (!"..." %in% fmls) {
    args <- args[names(args) %in% fmls]
  }
  do.call(fun, c(list(se), args))
}
