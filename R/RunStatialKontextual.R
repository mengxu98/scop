#' @title Run Statial Kontextual spatial relationships
#'
#' @description
#' Run `Statial::Kontextual()` on a spatial `Seurat` object to quantify
#' pairwise cell or spot label relationships relative to a parent context.
#' Results are stored as a compact SCOP bundle with raw Statial output,
#' standardized summary, and parameters. `Statial` is an optional Bioconductor
#' dependency installable with `BiocManager::install("Statial")`.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param group.by Metadata column containing cell or spot labels.
#' @param r Numeric radius or radii used by `Statial::Kontextual()`, expressed
#' in the selected coordinate units.
#' @param from,to,parent Cell or spot labels passed to `Statial::Kontextual()`.
#' Ignored when `parent_df` is supplied.
#' @param parent_df Optional data frame from `Statial::parentCombinations()`.
#' @param image Name of the Seurat spatial image. Required when multiple images
#' are present; a single image is selected automatically when `NULL`.
#' @param sample.by Optional metadata column used as Statial `imageID`. If
#' `NULL`, all cells or spots are treated as one image.
#' @param images Optional Statial image filter passed to `Kontextual(image = )`.
#' @param coord.cols Metadata coordinate columns used when no Seurat image
#' coordinates are available.
#' @param coordinate_space Coordinate system used for spatial relationships.
#' The default is raw acquisition coordinates; `"legacy_display"` remains an
#' explicit compatibility option.
#' @param inhom Whether Statial should account for inhomogeneity.
#' @param edge_correct Whether Statial should perform edge correction.
#' @param window,window.length Window arguments passed to
#' `Statial::Kontextual()`. Numeric window lengths use the selected coordinate
#' units.
#' @param include_original Whether to include original L-function values.
#' @param cores Number of cores passed to `Statial::Kontextual()`.
#' @param tool_name Name used to store results in `srt@tools`.
#' @param store_results Whether to store results in `srt@tools`.
#' @param store_input Whether to store the backend input cell table in
#' `srt@tools`.
#' @param ... Additional named arguments passed to `Statial::Kontextual()`.
#'
#' @return A `Seurat` object with Statial results stored in
#' `srt@tools[[tool_name]]` when `store_results = TRUE`.
#' @concept spatial-producer
#' @export
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- visium_human_pancreas_sub
#'
#' labels <- unique(as.character(spatial$coda_label))
#' if (length(labels) >= 2) {
#'   spatial <- RunStatialKontextual(
#'     spatial,
#'     group.by = "coda_label",
#'     r = 50,
#'     from = labels[1],
#'     to = labels[2],
#'     parent = labels[1:2],
#'     coord.cols = c("x", "y"),
#'     verbose = FALSE
#'   )
#'   spatial@tools$StatialKontextual$summary
#' }
RunStatialKontextual <- function(
  srt,
  group.by,
  r,
  from = NULL,
  to = NULL,
  parent = NULL,
  parent_df = NULL,
  image = NULL,
  sample.by = NULL,
  images = NULL,
  coord.cols = c("col", "row"),
  inhom = FALSE,
  edge_correct = TRUE,
  window = c("convex", "square", "concave"),
  window.length = NA_real_,
  include_original = TRUE,
  cores = 1,
  tool_name = "StatialKontextual",
  store_results = TRUE,
  store_input = FALSE,
  verbose = TRUE,
  coordinate_space = c("raw", "legacy_display"),
  ...
) {
  coordinate_space <- match.arg(coordinate_space)
  log_message(
    "Running Statial Kontextual spatial relationships",
    message_type = "running",
    verbose = verbose
  )
  statial_validate_srt(srt)
  statial_assert_string(group.by, "group.by")
  statial_assert_string(tool_name, "tool_name")
  statial_assert_flag(inhom, "inhom")
  statial_assert_flag(edge_correct, "edge_correct")
  statial_assert_flag(include_original, "include_original")
  statial_assert_flag(store_results, "store_results")
  statial_assert_flag(store_input, "store_input")
  window <- match.arg(window)
  r <- statial_validate_radius(r)
  cores <- statial_validate_positive_integer(cores, "cores")
  statial_validate_kontextual_relationship(
    parent_df = parent_df,
    from = from,
    to = to,
    parent = parent
  )
  extra_args <- list(...)
  statial_validate_named_list(extra_args, "...")

  cells <- statial_prepare_cells(
    srt = srt,
    group.by = group.by,
    image = image,
    sample.by = sample.by,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space,
    images = images,
    verbose = verbose
  )

  check_r("Statial", verbose = FALSE)
  kontextual <- get_namespace_fun("Statial", "Kontextual")
  args <- list(
    cells = cells,
    r = r,
    parentDf = parent_df,
    from = from,
    to = to,
    parent = parent,
    image = images,
    inhom = inhom,
    edgeCorrect = edge_correct,
    window = window,
    window.length = window.length,
    includeOriginal = include_original,
    spatialCoords = c("x", "y"),
    cellType = "cellType",
    imageID = "imageID",
    cores = cores
  )
  args <- utils::modifyList(args, extra_args)
  raw <- do.call(kontextual, args)
  table <- statial_standardize_kontextual(raw)
  summary <- statial_kontextual_summary(table)

  if (isTRUE(store_results)) {
    bundle <- list(
      table = table,
      raw = raw,
      summary = summary,
      parameters = list(
        group.by = group.by,
        r = r,
        from = from,
        to = to,
        parent = parent,
        parent_df = parent_df,
        image = image,
        sample.by = sample.by,
        coordinate_space = coordinate_space,
        images = images,
        coord.cols = coord.cols,
        inhom = inhom,
        edge_correct = edge_correct,
        window = window,
        window.length = window.length,
        include_original = include_original,
        cores = cores,
        tool_name = tool_name,
        store_results = store_results,
        store_input = store_input,
        backend_args = extra_args
      )
    )
    if (isTRUE(store_input)) {
      bundle$input <- cells
    }
    bundle <- spatial_result_build(
      bundle = bundle,
      method = "StatialKontextual",
      result_type = "neighborhood",
      provenance = list(producer = "RunStatialKontextual", backend_id = "statial")
    )
    srt@tools[[tool_name]] <- bundle
  }

  log_message(
    "{.pkg Statial} Kontextual completed {.val {nrow(table)}} relationship record{?s}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

statial_prepare_cells <- function(
  srt,
  group.by,
  image = NULL,
  sample.by = NULL,
  coord.cols = c("col", "row"),
  coordinate_space = c("raw", "legacy_display"),
  images = NULL,
  verbose = TRUE
) {
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} {.val {group.by}} is not present in metadata",
      message_type = "error"
    )
  }
  if (!is.null(sample.by)) {
    statial_assert_string(sample.by, "sample.by")
    if (!sample.by %in% colnames(srt@meta.data)) {
      log_message(
        "{.arg sample.by} {.val {sample.by}} is not present in metadata",
        message_type = "error"
      )
    }
  }
  coords <- spatial_analysis_coords(
    srt = srt,
    image = image,
    coord.cols = coord.cols,
    coordinate_space = coordinate_space
  )$data
  common <- colnames(srt)[colnames(srt) %in% rownames(coords)]
  if (length(common) == 0L) {
    log_message(
      "No cells or spots match spatial coordinates",
      message_type = "error"
    )
  }
  coords <- coords[common, , drop = FALSE]
  meta <- srt@meta.data[common, , drop = FALSE]
  sample <- if (is.null(sample.by)) {
    rep("sample1", length(common))
  } else {
    as.character(meta[[sample.by]])
  }
  cells <- data.frame(
    imageID = sample,
    cellType = as.character(meta[[group.by]]),
    x = as.numeric(coords$x),
    y = as.numeric(coords$y),
    row.names = common,
    stringsAsFactors = FALSE
  )
  keep <- !is.na(cells$imageID) &
    nzchar(cells$imageID) &
    !is.na(cells$cellType) &
    nzchar(cells$cellType) &
    is.finite(cells$x) &
    is.finite(cells$y)
  if (!all(keep)) {
    log_message(
      "Drop {.val {sum(!keep)}} cells or spots with missing Statial labels or coordinates",
      verbose = verbose
    )
    cells <- cells[keep, , drop = FALSE]
  }
  if (!is.null(images)) {
    cells <- cells[cells$imageID %in% images, , drop = FALSE]
  }
  if (nrow(cells) < 3L) {
    log_message(
      "At least three cells or spots are required for {.fn RunStatialKontextual}",
      message_type = "error"
    )
  }
  if (length(unique(cells$cellType)) < 2L) {
    log_message(
      "At least two {.arg group.by} labels are required for {.fn RunStatialKontextual}",
      message_type = "error"
    )
  }
  cells
}

statial_standardize_kontextual <- function(raw) {
  table <- as.data.frame(raw, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("imageID", "test", "kontextual", "r")
  missing <- setdiff(required, colnames(table))
  if (length(missing) > 0L) {
    log_message(
      "{.pkg Statial} Kontextual result is missing required column{?s}: {.val {missing}}",
      message_type = "error"
    )
  }
  table$imageID <- as.character(table$imageID)
  table$test <- as.character(table$test)
  table$kontextual <- suppressWarnings(as.numeric(table$kontextual))
  table$r <- suppressWarnings(as.numeric(table$r))
  if ("original" %in% colnames(table)) {
    table$original <- suppressWarnings(as.numeric(table$original))
  }
  if ("inhomL" %in% colnames(table)) {
    table$inhomL <- as.logical(table$inhomL)
  }
  rownames(table) <- NULL
  table
}

statial_kontextual_summary <- function(table, top_n = 10L) {
  valid <- table[is.finite(table$kontextual), , drop = FALSE]
  top <- valid[order(abs(valid$kontextual), decreasing = TRUE), , drop = FALSE]
  top <- utils::head(top, top_n)
  list(
    n_records = nrow(table),
    n_images = length(unique(table$imageID)),
    n_tests = length(unique(table$test)),
    radii = sort(unique(table$r)),
    n_localized = sum(valid$kontextual > 0, na.rm = TRUE),
    n_dispersed = sum(valid$kontextual < 0, na.rm = TRUE),
    top_relationships = top
  )
}

#' @title Plot stored Statial Kontextual results
#'
#' @description Plot contextual relationship scores across radii from a result
#' produced by [RunStatialKontextual()] without rerunning Statial.
#'
#' @param object Optional `Seurat` object containing the result.
#' @param res Optional result list, usually
#'   `object@tools$StatialKontextual`.
#' @param tests Optional relationship names to retain.
#' @param images Optional image identifiers to retain.
#' @return A `ggplot` object.
#' @export
StatialKontextualPlot <- function(object = NULL, res = NULL, tests = NULL, images = NULL) {
  if (is.null(res)) {
    if (is.null(object) || !inherits(object, "Seurat")) {
      log_message("Provide a {.cls Seurat} {.arg object} or a Statial {.arg res}", message_type = "error")
    }
    res <- object@tools[["StatialKontextual"]]
  }
  tab <- res$table %||% NULL
  if (!is.data.frame(tab) || !all(c("imageID", "test", "r", "kontextual") %in% colnames(tab))) {
    log_message("{.arg res} is not a plottable StatialKontextual result", message_type = "error")
  }
  if (!is.null(tests)) tab <- tab[tab$test %in% tests, , drop = FALSE]
  if (!is.null(images)) tab <- tab[tab$imageID %in% images, , drop = FALSE]
  if (nrow(tab) == 0L) log_message("No Statial records remain after filtering", message_type = "error")
  p <- ggplot2::ggplot(
    tab,
    ggplot2::aes(x = .data$r, y = .data$kontextual, color = .data$test, group = interaction(.data$imageID, .data$test))
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Radius", y = "Kontextual score", color = "Relationship") +
    theme_scop()
  if (length(unique(tab$imageID)) > 1L) p <- p + ggplot2::facet_wrap(~imageID)
  p
}

statial_validate_kontextual_relationship <- function(parent_df = NULL, from = NULL, to = NULL, parent = NULL) {
  if (!is.null(parent_df)) {
    parent_df <- as.data.frame(parent_df)
    missing <- setdiff(c("from", "to", "parent"), colnames(parent_df))
    if (length(missing) > 0L) {
      log_message(
        "{.arg parent_df} must contain column{?s}: {.val {missing}}",
        message_type = "error"
      )
    }
    return(invisible(TRUE))
  }
  missing_args <- c(
    if (is.null(from)) "from",
    if (is.null(to)) "to",
    if (is.null(parent)) "parent"
  )
  if (length(missing_args) > 0L) {
    log_message(
      "Provide {.arg parent_df} or all of {.arg from}, {.arg to}, and {.arg parent}. Missing: {.val {missing_args}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

statial_validate_radius <- function(r) {
  r <- suppressWarnings(as.numeric(r))
  if (length(r) == 0L || any(!is.finite(r) | r <= 0)) {
    log_message("{.arg r} must contain positive finite radii", message_type = "error")
  }
  unique(r)
}

statial_validate_positive_integer <- function(x, arg) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 1) {
    log_message("{.arg {arg}} must be a positive integer", message_type = "error")
  }
  as.integer(x)
}

statial_validate_srt <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  invisible(TRUE)
}

statial_assert_string <- function(x, arg) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    log_message("{.arg {arg}} must be a single non-empty string", message_type = "error")
  }
  invisible(TRUE)
}

statial_assert_flag <- function(x, arg) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    log_message("{.arg {arg}} must be TRUE or FALSE", message_type = "error")
  }
  invisible(TRUE)
}

statial_validate_named_list <- function(x, arg) {
  if (length(x) == 0L) {
    return(invisible(TRUE))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    log_message("{.arg {arg}} must contain named arguments only", message_type = "error")
  }
  invisible(TRUE)
}
