#' @title Run SpotSweeper spatial quality control
#'
#' @description
#' Run `SpotSweeper` spatially aware spot-level quality control on a spatial
#' `Seurat` object. The wrapper computes standard spot QC metrics, runs local
#' outlier detection for each metric, optionally runs regional artifact
#' detection per sample, and writes scop-style pass/fail metadata that can be
#' visualized with [SpatialSpotPlot()].
#'
#' @md
#' @param srt A `Seurat` object.
#' @param assay Assay used for expression. If `NULL`, the default assay is used.
#' @param layer Assay layer used for expression values.
#' @param coord.cols Metadata coordinate columns used when no Seurat image is
#' available.
#' @param image Name of the Seurat spatial image. Required when multiple images
#' are present; a single image is selected automatically when `NULL`.
#' @param coordinate_space Coordinate system used for spatial neighborhoods.
#' @param sample.by Optional metadata column identifying samples or images.
#' If `NULL`, all spots are treated as one sample.
#' @param metrics QC metrics used by `SpotSweeper::localOutliers()`. If `NULL`,
#' `nCount_<assay>`, `nFeature_<assay>`, and `percent.mito` are used.
#' @param directions Outlier direction for each metric. If `NULL`, count and
#' feature metrics use `"lower"` and mitochondrial metrics use `"higher"`.
#' @param n_neighbors Number of nearest spatial neighbors for local outlier
#' detection.
#' @param cutoff Modified z-score cutoff passed to local outlier detection.
#' @param log Whether SpotSweeper should log1p-transform local outlier metrics.
#' @param run_artifact Whether to run `SpotSweeper::findArtifacts()` per sample.
#' @param mito_pattern Regex prefixes used to identify mitochondrial genes.
#' @param mito_gene Optional explicit mitochondrial gene vector. When provided,
#' `mito_pattern` is ignored.
#' @param mito_percent Metadata column used as mitochondrial percent for
#' artifact detection. If `NULL`, `percent.mito` is used.
#' @param mito_sum Metadata column used as mitochondrial counts for artifact
#' detection. If `NULL`, mitochondrial counts are computed.
#' @param n_order,shape Parameters passed to `SpotSweeper::findArtifacts()`.
#' @param prefix Prefix used for metadata columns.
#' @param tool_name Name used to store detailed results in `srt@tools`.
#' @param return_filtered Whether to return only spots passing SpotSweeper QC.
#' @param store_results Whether to store detailed results in `srt@tools`.
#' @param workers Number of workers passed to SpotSweeper local outlier
#' detection.
#' @param verbose Whether to print progress messages.
#' @param ... Additional named arguments passed to matching SpotSweeper backend
#' functions when those arguments are supported by the installed version.
#'
#' @return A `Seurat` object with SpotSweeper QC metadata. When
#' `store_results = TRUE`, detailed results are stored in
#' `srt@tools[[tool_name]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#' spatial$SpotSweeper_QC <- factor(
#'   ifelse(seq_len(ncol(spatial)) %% 9 == 0, "Fail", "Pass"),
#'   levels = c("Pass", "Fail")
#' )
#' spatial$SpotSweeper_nCount_Spatial_z <- as.numeric(scale(spatial$nCount_Spatial))
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "SpotSweeper_QC",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "SpotSweeper_nCount_Spatial_z",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#'   spatial <- RunSpotSweeper(
#'     spatial,
#'     assay = "Spatial",
#'     coord.cols = c("x", "y"),
#'     n_neighbors = 12,
#'     run_artifact = FALSE,
#'     verbose = FALSE
#'   )
#'
#'   SpatialSpotPlot(
#'     spatial,
#'     group.by = "SpotSweeper_QC",
#'     overlay_image = FALSE,
#'     coord.cols = c("x", "y")
#'   )
RunSpotSweeper <- function(
  srt,
  assay = NULL,
  layer = "counts",
  coord.cols = c("col", "row"),
  image = NULL,
  sample.by = NULL,
  metrics = NULL,
  directions = NULL,
  n_neighbors = 36,
  cutoff = 3,
  log = TRUE,
  run_artifact = TRUE,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  mito_percent = NULL,
  mito_sum = NULL,
  n_order = 5,
  shape = c("hexagonal", "square"),
  prefix = "SpotSweeper",
  tool_name = "SpotSweeper",
  return_filtered = FALSE,
  store_results = TRUE,
  workers = 1,
  verbose = TRUE,
  coordinate_space = c("legacy_display", "raw"),
  ...
) {
  coordinate_space <- match.arg(coordinate_space)
  log_message(
    "Running SpotSweeper spatial quality control",
    message_type = "running",
    verbose = verbose
  )
  spot_sweeper_validate_srt(srt)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  shape <- match.arg(shape)
  spot_sweeper_assert_string(prefix, "prefix")
  spot_sweeper_assert_string(tool_name, "tool_name")
  spot_sweeper_assert_flag(log, "log")
  spot_sweeper_assert_flag(run_artifact, "run_artifact")
  spot_sweeper_assert_flag(return_filtered, "return_filtered")
  spot_sweeper_assert_flag(store_results, "store_results")
  n_neighbors <- spot_sweeper_assert_positive_integer(n_neighbors, "n_neighbors")
  n_order <- spot_sweeper_assert_positive_integer(n_order, "n_order")
  workers <- spot_sweeper_assert_positive_integer(workers, "workers")
  cutoff <- spot_sweeper_assert_number(cutoff, "cutoff")
  extra_args <- list(...)
  spot_sweeper_validate_named_list(extra_args, "...")

  check_r(
    c("MicTott/SpotSweeper", "SpatialExperiment", "SummarizedExperiment", "S4Vectors"),
    verbose = FALSE
  )

  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  inputs <- spot_sweeper_prepare_inputs(
    srt = srt,
    expr = expr,
    coords = coords,
    assay = assay,
    sample.by = sample.by,
    metrics = metrics,
    directions = directions,
    mito_pattern = mito_pattern,
    mito_gene = mito_gene,
    mito_percent = mito_percent,
    mito_sum = mito_sum,
    verbose = verbose
  )

  spe <- spot_sweeper_make_spe(
    expr = expr[inputs$features, inputs$spots, drop = FALSE],
    coldata = inputs$coldata,
    coords = inputs$coords
  )

  local_out <- spot_sweeper_run_local_outliers(
    spe = spe,
    metrics = inputs$metrics,
    directions = inputs$directions,
    sample_col = inputs$sample_col,
    n_neighbors = n_neighbors,
    cutoff = cutoff,
    log = log,
    workers = workers,
    extra_args = extra_args,
    verbose = verbose
  )
  spe <- local_out$spe

  artifact_out <- NULL
  if (isTRUE(run_artifact)) {
    artifact_out <- spot_sweeper_run_artifacts(
      spe = spe,
      mito_percent = inputs$mito_percent,
      mito_sum = inputs$mito_sum,
      sample_col = inputs$sample_col,
      n_order = n_order,
      shape = shape,
      log = log,
      extra_args = extra_args,
      verbose = verbose
    )
    spe <- artifact_out$spe
  }

  finalized <- spot_sweeper_finalize_metadata(
    srt = srt,
    spe = spe,
    prefix = prefix,
    metrics = inputs$metrics,
    metric_columns = inputs$metric_columns,
    run_artifact = run_artifact
  )
  srt <- Seurat::AddMetaData(srt, metadata = finalized$metadata)

  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      local_outliers = finalized$local_outliers,
      artifacts = finalized$artifacts,
      colData = finalized$coldata,
      coords = inputs$coords,
      metrics = inputs$metrics,
      directions = inputs$directions,
      metric_columns = inputs$metric_columns,
      tested_spots = inputs$spots,
      spe = spe,
      parameters = list(
        assay = assay,
        layer = layer,
        coord.cols = coord.cols,
        image = image,
        coordinate_space = coordinate_space,
        sample.by = sample.by,
        n_neighbors = n_neighbors,
        cutoff = cutoff,
        log = log,
        run_artifact = run_artifact,
        mito_pattern = mito_pattern,
        mito_gene = mito_gene,
        mito_percent = inputs$mito_percent,
        mito_sum = inputs$mito_sum,
        n_order = n_order,
        shape = shape,
        prefix = prefix,
        tool_name = tool_name,
        return_filtered = return_filtered,
        workers = workers,
        backend_args = extra_args
      ),
      backend = list(
        local_outliers = local_out$result,
        artifacts = if (is.null(artifact_out)) NULL else artifact_out$result
      )
    )
    srt@tools[[tool_name]] <- spatial_result_build(
      bundle = srt@tools[[tool_name]],
      method = "SpotSweeper",
      result_type = "quality_control",
      provenance = list(producer = "RunSpotSweeper", backend_id = "spotsweeper")
    )
  }

  failed <- sum(srt[[paste0(prefix, "_QC"), drop = TRUE]] == "Fail", na.rm = TRUE)
  log_message(
    "{.pkg SpotSweeper} flagged {.val {failed}} spots with prefix {.val {prefix}}",
    message_type = "success",
    verbose = verbose
  )

  if (isTRUE(return_filtered)) {
    srt <- srt[, srt[[paste0(prefix, "_QC"), drop = TRUE]] == "Pass"]
  }
  srt
}

spot_sweeper_prepare_inputs <- function(
  srt,
  expr,
  coords,
  assay,
  sample.by = NULL,
  metrics = NULL,
  directions = NULL,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  mito_percent = NULL,
  mito_sum = NULL,
  verbose = TRUE
) {
  common <- intersect(colnames(srt), rownames(coords))
  if (length(common) == 0L) {
    log_message(
      "No spatial coordinates match cells or spots in {.arg srt}",
      message_type = "error"
    )
  }
  coords <- coords[common, , drop = FALSE]
  keep_coords <- is.finite(coords$x) & is.finite(coords$y)
  if (!all(keep_coords)) {
    log_message(
      "Drop {.val {sum(!keep_coords)}} spots with missing or non-finite spatial coordinates",
      verbose = verbose
    )
  }
  coords <- coords[keep_coords, , drop = FALSE]
  spots <- rownames(coords)
  if (length(spots) < 3L) {
    log_message(
      "At least three spots with finite coordinates are required",
      message_type = "error"
    )
  }
  expr <- expr[, spots, drop = FALSE]
  features <- rownames(expr)
  if (length(features) == 0L) {
    log_message(
      "Selected assay layer has no features",
      message_type = "error"
    )
  }
  coldata <- as.data.frame(srt@meta.data[spots, , drop = FALSE], check.names = FALSE)
  sample_col <- spot_sweeper_resolve_sample_col(coldata, sample.by)
  if (is.null(sample.by)) {
    coldata[[sample_col]] <- "sample1"
  }

  metric_data <- spot_sweeper_metric_data(
    expr = expr,
    coldata = coldata,
    assay = assay,
    metrics = metrics,
    mito_pattern = mito_pattern,
    mito_gene = mito_gene,
    mito_percent = mito_percent,
    mito_sum = mito_sum
  )
  coldata <- metric_data$coldata
  metrics <- metric_data$metrics
  directions <- spot_sweeper_resolve_directions(metrics, directions)
  list(
    spots = spots,
    features = features,
    coords = coords,
    coldata = coldata,
    sample_col = sample_col,
    metrics = metrics,
    directions = directions,
    metric_columns = metric_data$metric_columns,
    mito_percent = metric_data$mito_percent,
    mito_sum = metric_data$mito_sum
  )
}

spot_sweeper_metric_data <- function(
  expr,
  coldata,
  assay,
  metrics = NULL,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  mito_percent = NULL,
  mito_sum = NULL
) {
  ncount_col <- paste0("nCount_", assay)
  nfeature_col <- paste0("nFeature_", assay)
  coldata[[ncount_col]] <- Matrix::colSums(expr)
  coldata[[nfeature_col]] <- Matrix::colSums(expr > 0)

  mito_features <- spot_qc_mito_features(
    features = rownames(expr),
    mito_pattern = mito_pattern,
    mito_gene = mito_gene
  )
  if (!is.null(mito_sum)) {
    spot_sweeper_assert_string(mito_sum, "mito_sum")
    if (!mito_sum %in% colnames(coldata)) {
      log_message(
        "{.arg mito_sum} {.val {mito_sum}} is not present in metadata",
        message_type = "error"
      )
    }
    mito_sum_col <- mito_sum
  } else {
    mito_sum_col <- "SpotSweeper_mito_sum"
    mito_counts <- rep(0, ncol(expr))
    names(mito_counts) <- colnames(expr)
    if (length(mito_features) > 0L) {
      mito_counts <- Matrix::colSums(expr[mito_features, , drop = FALSE])
    }
    coldata[[mito_sum_col]] <- mito_counts
  }

  if (!is.null(mito_percent)) {
    spot_sweeper_assert_string(mito_percent, "mito_percent")
    if (!mito_percent %in% colnames(coldata)) {
      log_message(
        "{.arg mito_percent} {.val {mito_percent}} is not present in metadata",
        message_type = "error"
      )
    }
    mito_percent_col <- mito_percent
  } else {
    mito_percent_col <- "percent.mito"
    ncount <- coldata[[ncount_col]]
    coldata[[mito_percent_col]] <- ifelse(ncount > 0, coldata[[mito_sum_col]] / ncount * 100, 0)
  }

  if (is.null(metrics)) {
    metrics <- c(ncount_col, nfeature_col, mito_percent_col)
  }
  metrics <- unique(as.character(metrics))
  missing_metrics <- setdiff(metrics, colnames(coldata))
  if (length(missing_metrics) > 0L) {
    log_message(
      "{.arg metrics} column{?s} not found in metadata: {.val {missing_metrics}}",
      message_type = "error"
    )
  }
  non_numeric <- metrics[!vapply(coldata[, metrics, drop = FALSE], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    log_message(
      "{.arg metrics} must be numeric. Non-numeric column{?s}: {.val {non_numeric}}",
      message_type = "error"
    )
  }
  list(
    coldata = coldata,
    metrics = metrics,
    metric_columns = c(nCount = ncount_col, nFeature = nfeature_col, mito_percent = mito_percent_col),
    mito_percent = mito_percent_col,
    mito_sum = mito_sum_col
  )
}

spot_sweeper_make_spe <- function(expr, coldata, coords) {
  check_r("SpatialExperiment", verbose = FALSE)
  spatial_experiment <- get_namespace_fun("SpatialExperiment", "SpatialExperiment")
  spatial_experiment(
    assays = list(counts = expr),
    colData = S4Vectors::DataFrame(coldata),
    spatialCoords = as.matrix(coords[, c("x", "y"), drop = FALSE])
  )
}

spot_sweeper_run_local_outliers <- function(
  spe,
  metrics,
  directions,
  sample_col,
  n_neighbors,
  cutoff,
  log,
  workers,
  extra_args,
  verbose = TRUE
) {
  cdata <- spot_sweeper_coldata(spe)
  sample_sizes <- table(cdata[[sample_col]])
  if (min(sample_sizes) < 2L) {
    log_message(
      "Each sample must contain at least two spots for SpotSweeper local outlier detection",
      message_type = "error"
    )
  }
  n_neighbors_use <- min(n_neighbors, min(sample_sizes) - 1L)
  if (n_neighbors_use < n_neighbors) {
    log_message(
      "Use {.arg n_neighbors = {n_neighbors_use}} for the smallest sample",
      verbose = verbose
    )
  }
  local_fun <- get_namespace_fun("SpotSweeper", "localOutliers")
  result <- list()
  for (metric in metrics) {
    log_message(
      "Run {.pkg SpotSweeper} local outlier detection for {.val {metric}}",
      verbose = verbose
    )
    args <- c(
      list(
        spe = spe,
        metric = metric,
        direction = directions[[metric]],
        n_neighbors = n_neighbors_use,
        samples = sample_col,
        log = log,
        cutoff = cutoff,
        workers = workers
      ),
      spot_sweeper_filter_args(local_fun, extra_args)
    )
    spe <- do.call(local_fun, args)
    result[[metric]] <- spot_sweeper_coldata(spe)
  }
  list(spe = spe, result = result)
}

spot_sweeper_run_artifacts <- function(
  spe,
  mito_percent,
  mito_sum,
  sample_col,
  n_order,
  shape,
  log,
  extra_args,
  verbose = TRUE
) {
  artifact_fun <- get_namespace_fun("SpotSweeper", "findArtifacts")
  cdata <- spot_sweeper_coldata(spe)
  samples <- unique(as.character(cdata[[sample_col]]))
  combined <- cdata
  combined$artifact <- FALSE
  combined$.SpotSweeper_artifact_skip_reason <- NA_character_
  result <- list()
  for (sample in samples) {
    idx <- which(as.character(cdata[[sample_col]]) == sample)
    spe_sample <- spe[, idx]
    sample_coldata <- spot_sweeper_coldata(spe_sample)
    artifact_ready <- spot_sweeper_artifact_ready(
      sample_coldata = sample_coldata,
      mito_percent = mito_percent,
      mito_sum = mito_sum
    )
    if (!isTRUE(artifact_ready$ready)) {
      log_message(
        "Skip {.pkg SpotSweeper} artifact detection for sample {.val {sample}}: {artifact_ready$reason}",
        verbose = verbose
      )
      sample_coldata$artifact <- FALSE
      sample_coldata$.SpotSweeper_artifact_skip_reason <- artifact_ready$reason
      combined <- spot_sweeper_merge_coldata(combined, sample_coldata)
      result[[sample]] <- sample_coldata
      next
    }
    log_message(
      "Run {.pkg SpotSweeper} artifact detection for sample {.val {sample}}",
      verbose = verbose
    )
    args <- c(
      list(
        spe = spe_sample,
        mito_percent = mito_percent,
        mito_sum = mito_sum,
        samples = sample_col,
        n_order = n_order,
        shape = shape,
        log = log,
        name = "artifact"
      ),
      spot_sweeper_filter_args(artifact_fun, extra_args)
    )
    spe_sample <- do.call(artifact_fun, args)
    sample_coldata <- spot_sweeper_coldata(spe_sample)
    sample_coldata$.SpotSweeper_artifact_skip_reason <- NA_character_
    combined <- spot_sweeper_merge_coldata(combined, sample_coldata)
    result[[sample]] <- sample_coldata
  }
  SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(combined[colnames(spe), , drop = FALSE])
  list(spe = spe, result = result)
}

spot_sweeper_finalize_metadata <- function(
  srt,
  spe,
  prefix,
  metrics,
  metric_columns,
  run_artifact = TRUE
) {
  cdata <- spot_sweeper_coldata(spe)
  all_spots <- colnames(srt)
  metadata <- srt@meta.data[, 0, drop = FALSE]
  local_table <- list()
  local_fail <- rep(FALSE, length(all_spots))
  names(local_fail) <- all_spots

  for (metric in metrics) {
    outlier_raw <- paste0(metric, "_outliers")
    z_raw <- paste0(metric, "_z")
    outlier_col <- spot_sweeper_metadata_col(prefix, metric, "outlier")
    z_col <- spot_sweeper_metadata_col(prefix, metric, "z")
    metadata[[outlier_col]] <- FALSE
    metadata[[z_col]] <- NA_real_
    outlier <- if (outlier_raw %in% colnames(cdata)) {
      as.logical(cdata[[outlier_raw]])
    } else {
      rep(FALSE, nrow(cdata))
    }
    z <- if (z_raw %in% colnames(cdata)) {
      suppressWarnings(as.numeric(cdata[[z_raw]]))
    } else {
      rep(NA_real_, nrow(cdata))
    }
    metadata[rownames(cdata), outlier_col] <- outlier
    metadata[rownames(cdata), z_col] <- z
    local_fail[rownames(cdata)] <- local_fail[rownames(cdata)] | outlier
    local_table[[metric]] <- data.frame(
      spot = rownames(cdata),
      metric = metric,
      outlier = outlier,
      z = z,
      stringsAsFactors = FALSE
    )
  }
  local_table <- do.call(rbind, local_table)
  rownames(local_table) <- NULL

  artifact_fail <- rep(FALSE, length(all_spots))
  names(artifact_fail) <- all_spots
  artifact_table <- NULL
  if (isTRUE(run_artifact) && "artifact" %in% colnames(cdata)) {
    artifact <- as.logical(cdata[["artifact"]])
    metadata[[paste0(prefix, "_artifact")]] <- FALSE
    metadata[rownames(cdata), paste0(prefix, "_artifact")] <- artifact
    artifact_fail[rownames(cdata)] <- artifact
    artifact_table <- data.frame(
      spot = rownames(cdata),
      artifact = artifact,
      stringsAsFactors = FALSE
    )
    var_cols <- grep("^k[0-9]+$", colnames(cdata), value = TRUE)
    if (length(var_cols) > 0L) {
      artifact_table <- cbind(artifact_table, cdata[, var_cols, drop = FALSE])
    }
    skip_col <- ".SpotSweeper_artifact_skip_reason"
    if (skip_col %in% colnames(cdata)) {
      artifact_table$skip_reason <- cdata[[skip_col]]
    }
  }

  metadata[[paste0(prefix, "_local_outlier_qc")]] <- spot_sweeper_pass_fail(local_fail)
  metadata[[paste0(prefix, "_artifact_qc")]] <- spot_sweeper_pass_fail(artifact_fail)
  metadata[[paste0(prefix, "_QC")]] <- spot_sweeper_pass_fail(local_fail | artifact_fail)

  list(
    metadata = metadata,
    local_outliers = local_table,
    artifacts = artifact_table,
    coldata = cdata,
    metric_columns = metric_columns
  )
}

spot_sweeper_merge_coldata <- function(combined, sample_coldata) {
  missing_combined <- setdiff(colnames(sample_coldata), colnames(combined))
  for (col in missing_combined) {
    combined[[col]] <- NA
  }
  missing_sample <- setdiff(colnames(combined), colnames(sample_coldata))
  for (col in missing_sample) {
    sample_coldata[[col]] <- NA
  }
  sample_coldata <- sample_coldata[, colnames(combined), drop = FALSE]
  combined[rownames(sample_coldata), ] <- sample_coldata
  combined
}

spot_sweeper_artifact_ready <- function(sample_coldata, mito_percent, mito_sum) {
  missing_cols <- setdiff(c(mito_percent, mito_sum), colnames(sample_coldata))
  if (length(missing_cols) > 0L) {
    return(list(
      ready = FALSE,
      reason = paste0("missing mitochondrial column(s): ", paste(missing_cols, collapse = ", "))
    ))
  }
  for (col in c(mito_percent, mito_sum)) {
    values <- suppressWarnings(as.numeric(sample_coldata[[col]]))
    values <- values[is.finite(values)]
    if (length(unique(values)) < 2L) {
      return(list(
        ready = FALSE,
        reason = paste0("column ", col, " has fewer than two finite values")
      ))
    }
  }
  list(ready = TRUE, reason = NA_character_)
}

spot_sweeper_resolve_sample_col <- function(coldata, sample.by = NULL) {
  if (is.null(sample.by)) {
    return(".SpotSweeper_sample")
  }
  spot_sweeper_assert_string(sample.by, "sample.by")
  if (!sample.by %in% colnames(coldata)) {
    log_message(
      "{.arg sample.by} {.val {sample.by}} is not present in metadata",
      message_type = "error"
    )
  }
  sample.by
}

spot_sweeper_resolve_directions <- function(metrics, directions = NULL) {
  if (is.null(directions)) {
    out <- ifelse(grepl("mito|percent", metrics, ignore.case = TRUE), "higher", "lower")
    names(out) <- metrics
    return(out)
  }
  directions <- as.character(directions)
  if (!is.null(names(directions)) && any(nzchar(names(directions)))) {
    missing <- setdiff(metrics, names(directions))
    if (length(missing) > 0L) {
      log_message(
        "{.arg directions} is missing metric{?s}: {.val {missing}}",
        message_type = "error"
      )
    }
    directions <- directions[metrics]
  } else if (length(directions) == 1L) {
    directions <- rep(directions, length(metrics))
  } else if (length(directions) != length(metrics)) {
    log_message(
      "{.arg directions} must have length 1, one value per metric, or be named by metric",
      message_type = "error"
    )
  }
  invalid <- setdiff(directions, c("lower", "higher", "both"))
  if (length(invalid) > 0L) {
    log_message(
      "{.arg directions} must contain only {.val lower}, {.val higher}, or {.val both}",
      message_type = "error"
    )
  }
  names(directions) <- metrics
  directions
}

spot_sweeper_filter_args <- function(fun, args) {
  if (length(args) == 0L) {
    return(args)
  }
  keep <- intersect(names(args), names(formals(fun)))
  args[keep]
}

spot_sweeper_coldata <- function(spe) {
  out <- as.data.frame(SummarizedExperiment::colData(spe), stringsAsFactors = FALSE)
  rownames(out) <- colnames(spe)
  out
}

spot_sweeper_metadata_col <- function(prefix, metric, suffix) {
  make.names(paste(prefix, metric, suffix, sep = "_"))
}

spot_sweeper_pass_fail <- function(fail) {
  factor(ifelse(fail, "Fail", "Pass"), levels = c("Pass", "Fail"))
}

spot_sweeper_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spot_sweeper_assert_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spot_sweeper_assert_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be TRUE or FALSE",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spot_sweeper_assert_positive_integer <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x)) {
    log_message(
      "{.arg {arg}} must be a positive integer",
      message_type = "error"
    )
  }
  as.integer(x)
}

spot_sweeper_assert_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    log_message(
      "{.arg {arg}} must be a single number",
      message_type = "error"
    )
  }
  as.numeric(x)
}

spot_sweeper_validate_named_list <- function(x, arg) {
  if (!is.list(x)) {
    log_message(
      "{.arg {arg}} must be a list",
      message_type = "error"
    )
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    log_message(
      "{.arg {arg}} entries must be named",
      message_type = "error"
    )
  }
  invisible(TRUE)
}
