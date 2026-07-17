#' @title Run CytoSPACE spatial assignment
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param reference Reference `Seurat` object containing annotated single cells.
#' @param reference_label Metadata column in `reference` with cell type labels.
#' @param reference_assay Assay used in `reference`.
#' @param layer,reference_layer Assay layers used for spatial and reference
#' expression.
#' @param features Features used for assignment. If `NULL`, shared features are
#' used.
#' @param cell_fractions Optional cell-type fractions. Provide a named numeric
#' vector, one-row matrix/data.frame, or a spot-by-cell-type matrix/data.frame.
#' Spot-level rows are aggregated to the global composition used by the default
#' CytoSPACE assignment workflow.
#' @param n_cells_per_spot Optional number of cells assigned to each spatial
#' spot. If `NULL`, counts are estimated from spatial RNA reads with
#' `mean_cell_numbers`.
#' @param mean_cell_numbers Mean number of cells per spot. Default `5`, matching
#' the CytoSPACE Visium default.
#' @param scRNA_max_transcripts_per_cell Maximum reference transcripts per cell
#' before assignment. Default `1500`, matching CytoSPACE.
#' @param sampling_method Sampling method. Only `"duplicates"` is supported in
#' the package runtime.
#' @param seed Random seed used for deterministic reference downsampling and
#' duplicate sampling.
#' @param prefix Prefix for metadata columns.
#' @param store_results Whether to store detailed assignment results in
#' `srt@tools`.
#' @param image Optional Seurat image used for spatial coordinates.
#' @param coord.cols Metadata coordinate columns used when no image is selected.
#' @param coordinate_space Coordinate space used for assignment locations. The
#'   default preserves the historical coordinate behavior.
#'
#' @return A `Seurat` object with CytoSPACE metadata columns and detailed
#' results stored in `srt@tools[["CytoSPACE"]]`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' data(pancreas_sub)
#' features_use <- intersect(
#'   rownames(visium_human_pancreas_sub),
#'   rownames(pancreas_sub)
#' )
#' spatial <- RunCytoSPACE(
#'   visium_human_pancreas_sub,
#'   reference = pancreas_sub,
#'   reference_label = "CellType",
#'   features = features_use,
#'   mean_cell_numbers = 1
#' )
#'
#' SpatialSpotPlot(
#'   visium_human_pancreas_sub,
#'   group.by = "coda_label",
#'   theme_use = "theme_scop"
#' )
#'
#' SpatialSpotPlot(
#'   spatial,
#'   group.by = "CytoSPACE_dominant_type",
#'   theme_use = "theme_scop"
#' )
RunCytoSPACE <- function(
  srt,
  reference,
  reference_label,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  cell_fractions = NULL,
  n_cells_per_spot = NULL,
  mean_cell_numbers = 5,
  scRNA_max_transcripts_per_cell = 1500,
  sampling_method = "duplicates",
  seed = 1,
  prefix = "CytoSPACE",
  store_results = TRUE,
  verbose = TRUE,
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("legacy_display", "raw")
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!inherits(reference, "Seurat")) {
    log_message(
      "{.arg reference} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  if (!identical(sampling_method, "duplicates")) {
    log_message(
      "{.arg sampling_method} only supports {.val duplicates} in the native package runtime",
      message_type = "error"
    )
  }
  coordinate_space <- match.arg(coordinate_space)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(reference)

  labels <- cytospace_resolve_reference_labels(reference, reference_label)
  keep_ref <- !is.na(labels)
  if (!all(keep_ref)) {
    log_message(
      "Drop {.val {sum(!keep_ref)}} reference cells with missing {.arg reference_label}",
      verbose = verbose
    )
    reference <- reference[, keep_ref]
    labels <- labels[keep_ref]
  }
  cell_types <- unique(as.character(labels))
  labels <- factor(labels, levels = cell_types)
  if (nlevels(labels) < 1L) {
    log_message(
      "{.arg reference_label} must contain at least one non-missing class",
      message_type = "error"
    )
  }

  features_use <- cytospace_common_features(
    srt,
    reference,
    assay,
    reference_assay,
    features
  )
  if (length(features_use) == 0L) {
    log_message(
      "No shared features are available for {.fn RunCytoSPACE}",
      message_type = "error"
    )
  }

  st_expr <- cytospace_get_expr_matrix(srt, assay, layer, features_use)
  ref_expr <- cytospace_get_expr_matrix(
    reference,
    reference_assay,
    reference_layer,
    features_use
  )
  st_expr <- st_expr[features_use, , drop = FALSE]
  ref_expr <- ref_expr[features_use, , drop = FALSE]

  spot_ids <- colnames(st_expr)
  n_cells <- cytospace_resolve_n_cells_per_spot(
    n_cells_per_spot = n_cells_per_spot,
    spot_ids = spot_ids,
    st_expr = st_expr,
    mean_cell_numbers = mean_cell_numbers
  )
  fractions <- cytospace_resolve_cell_fractions(
    cell_fractions = cell_fractions,
    st_expr = st_expr,
    ref_expr = ref_expr,
    labels = labels,
    cell_types = cell_types,
    n_cells_per_spot = n_cells
  )
  target_cell_type_counts <- cytospace_fraction_to_cell_type_counts(
    fractions,
    sum(n_cells)
  )
  log_message(
    "Estimated CytoSPACE cell type fractions for {.val {length(cell_types)}} reference labels",
    verbose = verbose
  )
  ref_sample <- cytospace_sample_reference_cells(
    ref_expr = ref_expr,
    labels = labels,
    cell_types = cell_types,
    target_cell_type_counts = target_cell_type_counts,
    scRNA_max_transcripts_per_cell = scRNA_max_transcripts_per_cell,
    seed = seed
  )

  coords <- cytospace_get_spatial_coords(
    srt,
    spot_ids,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )

  log_message(
    "Run CytoSPACE with {.val {length(features_use)}} features, {.val {ncol(ref_sample$expr)}} sampled reference cells, and {.val {ncol(st_expr)}} spatial spots",
    verbose = verbose
  )

  result <- cytospace_assign(
    sc_expr = ref_sample$expr,
    st_expr = st_expr,
    spot_capacities = as.integer(n_cells),
    seed = seed,
    verbose = verbose
  )
  log_message(
    "Build CytoSPACE assignment summaries",
    verbose = verbose
  )

  assignments <- cytospace_build_assignment_table(
    result = result,
    sampled_cells = ref_sample$cell_ids,
    sampled_labels = ref_sample$labels,
    spot_ids = spot_ids,
    coords = coords
  )
  spot_summary <- cytospace_build_spot_summary(
    assignments,
    spot_ids,
    cell_types
  )
  assigned_expression <- cytospace_build_assigned_expression(
    sampled_expr = ref_sample$expr,
    assignments = assignments
  )
  srt <- cytospace_add_metadata(srt, spot_summary, prefix)

  if (isTRUE(store_results)) {
    srt@tools[["CytoSPACE"]] <- list(
      assigned_locations = assignments,
      cell_type_assignments_by_spot = spot_summary$counts,
      fractional_abundances_by_spot = spot_summary$fractions,
      assigned_expression = assigned_expression,
      spatial_coords = coords,
      n_cells_per_spot = n_cells,
      cell_fractions = fractions,
      target_cell_type_counts = target_cell_type_counts,
      features = features_use,
      cell_types = cell_types,
      parameters = list(
        assay = assay,
        reference_assay = reference_assay,
        layer = layer,
        reference_layer = reference_layer,
          reference_label = reference_label,
          image = image,
          coord.cols = coord.cols,
          coordinate_space = coordinate_space,
        mean_cell_numbers = mean_cell_numbers,
        scRNA_max_transcripts_per_cell = scRNA_max_transcripts_per_cell,
        sampling_method = sampling_method,
        distance_metric = "Pearson_correlation",
        seed = seed,
        prefix = prefix
      )
    )
    srt@tools[["CytoSPACE"]] <- spatial_result_build(
      bundle = srt@tools[["CytoSPACE"]],
      method = "CytoSPACE",
      result_type = "mapping",
      source = c(
        attr(coords, "spatial_source") %||% list(),
        list(transform = attr(coords, "spatial_transform"))
      ),
      provenance = list(producer = "RunCytoSPACE", backend_id = "core")
    )
  }

  log_message(
    "{.pkg CytoSPACE} assignments stored in {.code srt@tools[['CytoSPACE']]} and metadata prefix {.val {prefix}}",
    verbose = verbose
  )
  srt
}

cytospace_resolve_reference_labels <- function(reference, reference_label) {
  if (
    missing(reference_label) ||
      is.null(reference_label) ||
      length(reference_label) != 1L
  ) {
    log_message(
      "{.arg reference_label} must be a single reference metadata column",
      message_type = "error"
    )
  }
  if (!reference_label %in% colnames(reference[[]])) {
    log_message(
      "{.arg reference_label} {.val {reference_label}} is not present in {.arg reference}",
      message_type = "error"
    )
  }
  reference[[reference_label, drop = TRUE]]
}

cytospace_common_features <- function(
  srt,
  reference,
  assay,
  reference_assay,
  features = NULL
) {
  common <- intersect(
    rownames(srt[[assay]]),
    rownames(reference[[reference_assay]])
  )
  if (!is.null(features)) {
    common <- intersect(features, common)
  }
  common
}

cytospace_get_expr_matrix <- function(srt, assay, layer, features) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  mat <- mat[features, , drop = FALSE]
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}

cytospace_resolve_n_cells_per_spot <- function(
  n_cells_per_spot,
  spot_ids,
  st_expr,
  mean_cell_numbers = 5
) {
  if (!is.null(n_cells_per_spot)) {
    if (!is.numeric(n_cells_per_spot)) {
      log_message(
        "{.arg n_cells_per_spot} must be numeric",
        message_type = "error"
      )
    }
    if (!is.null(names(n_cells_per_spot))) {
      missing_spots <- setdiff(spot_ids, names(n_cells_per_spot))
      if (length(missing_spots) > 0L) {
        log_message(
          "{.arg n_cells_per_spot} is missing {.val {length(missing_spots)}} spatial spots",
          message_type = "error"
        )
      }
      n_cells_per_spot <- n_cells_per_spot[spot_ids]
    }
    if (length(n_cells_per_spot) != length(spot_ids)) {
      log_message(
        "{.arg n_cells_per_spot} must have one value per spatial spot",
        message_type = "error"
      )
    }
    n_cells <- as.integer(n_cells_per_spot)
    if (any(is.na(n_cells) | n_cells < 0L)) {
      log_message(
        "{.arg n_cells_per_spot} cannot contain missing or negative values",
        message_type = "error"
      )
    }
    names(n_cells) <- spot_ids
    return(n_cells)
  }

  if (
    !is.numeric(mean_cell_numbers) ||
      length(mean_cell_numbers) != 1L ||
      is.na(mean_cell_numbers) ||
      mean_cell_numbers < 1
  ) {
    log_message(
      "{.arg mean_cell_numbers} must be a single number >= 1",
      message_type = "error"
    )
  }
  norm_st <- cytospace_log_cpm_r(st_expr)
  rna_reads <- colSums(norm_st)
  mean_reads <- mean(rna_reads)
  min_reads <- min(rna_reads)
  min_cells <- if (min_reads > 0) 1 else 0
  if (
    !is.finite(mean_reads) || abs(mean_reads - min_reads) < .Machine$double.eps
  ) {
    n_cells <- rep(as.integer(mean_cell_numbers), length(spot_ids))
  } else {
    slope <- (mean_cell_numbers - min_cells) / (mean_reads - min_reads)
    intercept <- min_cells - slope * min_reads
    n_cells <- as.integer(slope * rna_reads + intercept)
  }
  n_cells[is.na(n_cells) | n_cells < 0L] <- 0L
  names(n_cells) <- spot_ids
  n_cells
}

cytospace_resolve_cell_fractions <- function(
  cell_fractions,
  st_expr,
  ref_expr,
  labels,
  cell_types,
  n_cells_per_spot
) {
  if (!is.null(cell_fractions)) {
    fractions <- cytospace_parse_cell_fractions(
      cell_fractions,
      cell_types,
      colnames(st_expr),
      n_cells_per_spot
    )
  } else {
    fractions <- cytospace_estimate_fractions(
      st_expr,
      ref_expr,
      labels,
      cell_types,
      n_cells_per_spot
    )
  }
  fractions[!is.finite(fractions) | fractions < 0] <- 0
  if (sum(fractions) <= 0) {
    fractions[] <- 1 / length(fractions)
  }
  fractions <- fractions / sum(fractions)
  names(fractions) <- cell_types
  fractions
}

cytospace_parse_cell_fractions <- function(
  cell_fractions,
  cell_types,
  spot_ids,
  n_cells_per_spot
) {
  if (is.numeric(cell_fractions) && is.null(dim(cell_fractions))) {
    fractions <- cell_fractions
    if (is.null(names(fractions))) {
      if (length(fractions) != length(cell_types)) {
        log_message(
          "{.arg cell_fractions} must be named or have one value per reference label",
          message_type = "error"
        )
      }
      names(fractions) <- cell_types
    }
    missing_types <- setdiff(cell_types, names(fractions))
    if (length(missing_types) > 0L) {
      log_message(
        "{.arg cell_fractions} must contain all reference labels",
        message_type = "error"
      )
    }
    return(as.numeric(fractions[cell_types]))
  }

  mat <- as.matrix(cell_fractions)
  storage.mode(mat) <- "double"
  if (is.null(colnames(mat))) {
    if (ncol(mat) != length(cell_types)) {
      log_message(
        "{.arg cell_fractions} must have one column per reference label",
        message_type = "error"
      )
    }
    colnames(mat) <- cell_types
  }
  missing_types <- setdiff(cell_types, colnames(mat))
  if (length(missing_types) > 0L) {
    log_message(
      "{.arg cell_fractions} must contain all reference labels",
      message_type = "error"
    )
  }
  mat <- mat[, cell_types, drop = FALSE]
  if (nrow(mat) == 1L) {
    return(as.numeric(mat[1, ]))
  }
  if (is.null(rownames(mat))) {
    if (nrow(mat) != length(spot_ids)) {
      log_message(
        "{.arg cell_fractions} must have one row or one row per spatial spot",
        message_type = "error"
      )
    }
    rownames(mat) <- spot_ids
  }
  missing_spots <- setdiff(spot_ids, rownames(mat))
  if (length(missing_spots) > 0L) {
    log_message(
      "{.arg cell_fractions} must contain all spatial spots",
      message_type = "error"
    )
  }
  mat <- mat[spot_ids, , drop = FALSE]
  colSums(mat * n_cells_per_spot) / sum(n_cells_per_spot)
}

cytospace_estimate_fractions <- function(
  st_expr,
  ref_expr,
  labels,
  cell_types,
  n_cells_per_spot
) {
  norm_st <- cytospace_log_cpm_r(st_expr)
  norm_ref <- cytospace_log_cpm_r(ref_expr)
  profiles <- vapply(
    cell_types,
    function(ct) {
      idx <- which(labels == ct)
      if (length(idx) == 1L) {
        norm_ref[, idx]
      } else {
        rowMeans(norm_ref[, idx, drop = FALSE])
      }
    },
    numeric(nrow(norm_ref))
  )
  rownames(profiles) <- rownames(norm_ref)
  scores <- suppressWarnings(stats::cor(
    x = norm_st,
    y = profiles,
    method = "pearson"
  ))
  scores <- as.matrix(scores)
  scores[!is.finite(scores) | scores < 0] <- 0
  rownames(scores) <- colnames(norm_st)
  colnames(scores) <- cell_types
  row_sums <- rowSums(scores)
  scores[row_sums <= 0, ] <- 1 / ncol(scores)
  row_sums <- rowSums(scores)
  scores <- sweep(scores, 1, row_sums, "/")
  colSums(scores * n_cells_per_spot) / sum(n_cells_per_spot)
}

cytospace_log_cpm_r <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  totals <- colSums(x)
  out <- x
  keep <- is.finite(totals) & totals > 0
  out[, keep] <- log2(
    sweep(out[, keep, drop = FALSE], 2, totals[keep], "/") * 1e6 + 1
  )
  out[, !keep] <- 0
  out[!is.finite(out)] <- 0
  out
}

cytospace_fraction_to_cell_type_counts <- function(fractions, total_cells) {
  raw <- fractions * total_cells
  counts <- as.integer(raw)
  remain <- total_cells - sum(counts)
  if (remain > 0L) {
    counts[1L] <- counts[1L] + remain
  }
  names(counts) <- names(fractions)
  counts
}

cytospace_downsample_counts <- function(x, target_count, seed = 1) {
  if (is.null(target_count) || is.na(target_count) || target_count <= 0) {
    return(x)
  }
  set.seed(seed)
  out <- x
  for (j in seq_len(ncol(out))) {
    total <- sum(out[, j])
    if (!is.finite(total) || total <= target_count) {
      next
    }
    prob <- out[, j] / total
    out[, j] <- as.numeric(stats::rmultinom(
      1L,
      size = target_count,
      prob = prob
    ))
  }
  out
}

cytospace_sample_reference_cells <- function(
  ref_expr,
  labels,
  cell_types,
  target_cell_type_counts,
  scRNA_max_transcripts_per_cell = 1500,
  seed = 1
) {
  ref_use <- cytospace_downsample_counts(
    ref_expr,
    scRNA_max_transcripts_per_cell,
    seed
  )
  set.seed(seed)
  sampled_index <- integer()
  for (ct in cell_types) {
    idx <- which(labels == ct)
    desired <- target_cell_type_counts[ct]
    if (length(idx) == 0L && desired > 0L) {
      log_message(
        "Cell type {.val {ct}} is absent from the reference object",
        message_type = "error"
      )
    }
    if (desired <= 0L) {
      next
    }
    if (desired > length(idx)) {
      chosen <- c(idx, sample(idx, desired - length(idx), replace = TRUE))
    } else {
      chosen <- sample(idx, desired, replace = FALSE)
    }
    sampled_index <- c(sampled_index, chosen)
  }
  expr <- ref_use[, sampled_index, drop = FALSE]
  colnames(expr) <- make.unique(colnames(ref_expr)[sampled_index], sep = "_dup")
  list(
    expr = expr,
    source_index = sampled_index,
    cell_ids = colnames(ref_expr)[sampled_index],
    labels = as.character(labels[sampled_index])
  )
}

cytospace_get_spatial_coords <- function(
  srt,
  spot_ids,
  image = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("legacy_display", "raw")
) {
  coordinate_space <- match.arg(coordinate_space)
  resolved <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space,
    image_policy = "strict"
  )
  matched <- match(spot_ids, resolved$data$cell_id)
  if (anyNA(matched)) {
    log_message("Spatial coordinates are missing for one or more requested spots", message_type = "error")
  }
  coords <- resolved$data[matched, c("x", "y"), drop = FALSE]
  rownames(coords) <- spot_ids
  attr(coords, "spatial_source") <- resolved$source
  attr(coords, "spatial_transform") <- resolved$transform
  coords
}

cytospace_build_assignment_table <- function(
  result,
  sampled_cells,
  sampled_labels,
  spot_ids,
  coords
) {
  cell_index <- result[["cell_index"]]
  spot_index <- result[["spot_index"]]
  assigned <- data.frame(
    UniqueCID = paste0(
      "UCID",
      formatC(
        seq_along(cell_index) - 1L,
        width = nchar(length(cell_index)),
        flag = "0"
      )
    ),
    OriginalCID = sampled_cells[cell_index],
    CellType = sampled_labels[cell_index],
    SpotID = spot_ids[spot_index],
    score = result[["score"]],
    stringsAsFactors = FALSE
  )
  coords_use <- coords[assigned$SpotID, , drop = FALSE]
  rownames(coords_use) <- NULL
  cbind(assigned, coords_use)
}

cytospace_build_spot_summary <- function(assignments, spot_ids, cell_types) {
  counts <- table(
    factor(assignments$SpotID, levels = spot_ids),
    factor(assignments$CellType, levels = cell_types)
  )
  counts <- as.data.frame.matrix(counts)
  counts[["Total cells"]] <- rowSums(counts)
  fractions <- counts[, cell_types, drop = FALSE]
  totals <- counts[["Total cells"]]
  fractions[] <- lapply(fractions, function(x) {
    ifelse(totals > 0, x / totals, 0)
  })
  list(counts = counts, fractions = fractions)
}

cytospace_build_assigned_expression <- function(sampled_expr, assignments) {
  expr <- sampled_expr[, seq_len(nrow(assignments)), drop = FALSE]
  colnames(expr) <- assignments$UniqueCID
  expr
}

cytospace_add_metadata <- function(srt, spot_summary, prefix = "CytoSPACE") {
  counts <- spot_summary$counts
  fractions <- spot_summary$fractions
  meta <- data.frame(
    total_cells = counts[["Total cells"]],
    row.names = rownames(counts),
    check.names = FALSE
  )
  colnames(meta) <- paste0(prefix, "_total_cells")
  for (ct in colnames(fractions)) {
    suffix <- make.names(ct)
    meta[[paste0(prefix, "_count_", suffix)]] <- counts[[ct]]
    meta[[paste0(prefix, "_frac_", suffix)]] <- fractions[[ct]]
  }
  if (ncol(fractions) > 0L) {
    meta[[paste0(prefix, "_dominant_type")]] <- colnames(fractions)[
      max.col(as.matrix(fractions), ties.method = "first")
    ]
  }
  Seurat::AddMetaData(srt, metadata = meta)
}
